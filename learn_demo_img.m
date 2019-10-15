function [ model ] = learn_demo_img( dataset, modelPath, conf )
% learning function

%load(getDatasetPath(conf), 'classes', 'dataset');
%classid = find(strcmp(classes,cls));
load(modelPath, 'model');
% ------------------------------------------------------------------
%                                                    Run SVM struct
% ------------------------------------------------------------------

svm = model.svm;
sizeFeatV = size(svm.wV);
sizeFeatE = size(svm.wE);
nDimAttr = [ numel(svm.rho_bins)  numel(svm.theta_bins) ];

switch lower(conf.learningType)
    case 'harg' % full learning
        parm.dimension = prod(sizeFeatV) + prod(sizeFeatE);
        parm.featureFn = @featureCB;
        parm.constraintFn  = @constraintCB;
    case 'dw' % indivisual weight learning
        parm.dimension = prod(sizeFeatV) + 2*size(svm.wE,2); % distance & angle
        parm.featureFn = @featureCB_dw;
        parm.constraintFn  = @constraintCB_dw;
    case 'sw' % shared weight learning
        parm.dimension = sizeFeatV(1) + 2;
        parm.featureFn = @featureCB_sw;
        parm.constraintFn  = @constraintCB_sw;
    case 'non'
        return;
end  

% initialize training and test examples for the given class
exmpAll = makeExampleSet_img(model.classid, dataset);
exmpTrain = exmpAll( find([ exmpAll.type ] == 0) );
%exmpTest = exmpAll( find([ exmpAll.type ] == 1) );

switch lower(conf.learningType)
    case { 'harg', 'dw', 'sw' } % SSVM optimization
        patterns = {} ;
        labels = {} ;
        for i = 1:numel(exmpTrain)
            patterns{i} = exmpTrain(i) ;       % x

            y.cls = conf.class;
            y.seq = exmpTrain(i).seqGT;
            y.pt = exmpTrain(i).ptGT;
            y.frame = exmpTrain(i).frameGT;
            y.distThres = conf.distRatioTrue*sqrt(sum(exmpTrain(i).object_size.^2));
            %y.distThres = conf.distRatioTrue*sqrt(sum(exmpTrain(i).size.^2));            
            y.modelPath = modelPath;

            y.sizeFeatV = sizeFeatV;
            y.sizeFeatE = sizeFeatE;
            y.nDimAttr = nDimAttr;
            labels{i} = y;   % y
        end

        parm.patterns = patterns ;
        parm.labels = labels ;
        parm.lossFn = @lossCB ;
        parm.verbose = conf.verbose;
        output = svm_struct_learn(' -c 1.0 -o 2 -p 1 -v 3 -e 0.02 -# 30', parm) ;
                
end
        

%fprintf('positive: %.2f (%d elements)\n', sum(output.w(output.w>eps)), nnz(output.w>eps) );
%fprintf('negative: %.2f (%d elements)\n', sum(output.w(output.w<=eps)), nnz(output.w<=eps) );
%output.w(output.w<=eps) = 0;
%pause;

switch lower(conf.learningType)
    case 'harg' % full learning
        %output.w = output.w - min(output.w);
        svm.wV = single(reshape(output.w(1:prod(sizeFeatV)), sizeFeatV));
        svm.wE = single(reshape(output.w(prod(sizeFeatV)+1:end), sizeFeatE));
    case 'dw' % indivisual weight learning
        nfV = prod(sizeFeatV);
        svm.wV = svm.wV.*reshape(output.w(1:nfV),sizeFeatV);
        svm.wE = svm.wE.* [ repmat( output.w( nfV+1:nfV+sizeFeatE(2) )', [ nDimAttr(1) 1 ]);...
            repmat(output.w(nfV+sizeFeatE(2)+1:end)', [ nDimAttr(2) 1 ]) ];
    case 'sw' % shared weight learning
        svm.wV = svm.wV.*repmat(output.w(1:sizeFeatV(1)),[1 sizeFeatV(2)]);
        svm.wE = svm.wE.*[ output.w(sizeFeatV(1)+1)*ones([ nDimAttr(1) sizeFeatE(2) ]);...
            output.w(sizeFeatV(1)+2)*ones([ nDimAttr(2) sizeFeatE(2) ]) ];
end 

model.svm = svm;

end


% ------------------------------------------------------------------
%                                               SVM struct callbacks
% ------------------------------------------------------------------

function delta = lossCB(param, y, ybar)
  %delta = double(y ~= ybar) ;
  % a normalized Hamming loss for matching
  %delta = 1 - dot(y,ybar)/sum(y); 

  nPt = size(y.pt,2);
  sqDist = sum((y.pt - ybar.pt).^2,1);
  delta = nnz( sqDist > y.distThres^2 )/nPt; 
  %delta = sum(min(1.0, sqDist./y.distThres^2))/nPt; % smooth loss
  if param.verbose
    %fprintf('loss(%3d,%3d) = %.2f\n', numel(y_seq), numel(ybar_seq), delta) ;
    %fprintf('loss: %.2f\n', delta) ;
  end
end

function psi = featureCB(param, x, y)
  % joint feature map of test graph (here, x) and assignment y
  view =  extract_localfeatures_wrapper( x.fname, x.tname );
  fV = view.desc(:,y.seq);
  % this part could be speeded up
  edgeAttr_t = computeEdgeAttr( view.frame(:,y.seq) );
  fE = makeFeatBin(edgeAttr_t);  

  psi = sparse(double( [fV(:); fE(:)]));
%   if param.verbose
%     fprintf('w = psi([%8.3f,%8.3f], %3d) = [%8.3f, %8.3f]\n', ...
%             x, y, full(psi(1)), full(psi(2))) ;
%   end
end


function psi = featureCB_dw(param, x, y)
  % joint feature map of test graph (here, x) and assignment y
  view =  extract_localfeatures_wrapper( x.fname, x.tname );
  fV = view.desc(:,y.seq);
  % this part could be speeded up
  edgeAttr_t = computeEdgeAttr( view.frame(:,y.seq) );
  fE = makeFeatBin(edgeAttr_t);
  
  load(y.modelPath, 'model');
  
  wV = model.svm.wV.*fV;  
  smatE = model.svm.wE.*fE;
  rho_wE = double(sum(smatE(1:y.nDimAttr(1),:),1));
  theta_wE = double(sum(smatE(y.nDimAttr(1)+1:sum(y.nDimAttr(1:2)),:),1));  
  psi = sparse(double([ wV(:); rho_wE(:); theta_wE(:)]));  
end

function psi = featureCB_sw(param, x, y)
  % joint feature map of test graph (here, x) and assignment y
  view =  extract_localfeatures_wrapper( x.fname, x.tname );
  fV = view.desc(:,y.seq);
  % this part could be speeded up
  edgeAttr_t = computeEdgeAttr( view.frame(:,y.seq) );
  fE = makeFeatBin(edgeAttr_t);
  
  load(y.modelPath, 'model');
  
  wV = sum(model.svm.wV.*fV,2);  
  smatE = model.svm.wE.*fE;
  rho_wE = double(sum(sum(smatE(1:y.nDimAttr(1),:))));
  theta_wE = double(sum(sum(smatE(y.nDimAttr(1)+1:sum(y.nDimAttr(1:2)),:))));  
  psi = sparse(double([wV; rho_wE; theta_wE]));  
end


function yhat = constraintCB(param, model_lib, x, y)
% function to find the most violating constraint which is compatible with w

% margin rescaling: argmax_y delta(yi, y) + <psi(x,y), w>
    load(y.modelPath, 'model');  
    svm = model.svm;
    svm.wV = single(reshape(model_lib.w(1:prod(y.sizeFeatV)), y.sizeFeatV));
    svm.wE = single(reshape(model_lib.w(prod(y.sizeFeatV)+1:end), y.sizeFeatE));
    [ res ] = immatching_FMVC_wrapper( x, model.view_m, svm, y );    
    yhat = y;
    if nnz(res.X) == numel(y.seq)
        yhat.seq = res.seq;
        yhat.pt = res.frame(1:2,res.seq);
        yhat.frame = res.frame(:,res.seq);
    else
        rand_y_seg = randperm(size(res.pt,2));
        yhat.seq = rand_y_seg(1:numel(y.seq));
        yhat.pt = res.frame(1:2,yhat.seq);
        yhat.frame = res.frame(:,yhat.seq);
    end
  
end

function yhat = constraintCB_dw(param, model_lib, x, y)
% function to find the most violating constraint which is compatible with w

% margin rescaling: argmax_y delta(yi, y) + <psi(x,y), w>
  load(y.modelPath, 'model');  
    svm = model.svm;
    nfV = prod(y.sizeFeatV);
    svm.wV = svm.wV.*reshape(model_lib.w(1:nfV),y.sizeFeatV);
    svm.wE = svm.wE.* [ repmat( model_lib.w( nfV+1:nfV+y.sizeFeatE(2) )', [ y.nDimAttr(1) 1 ]);...
        repmat(model_lib.w(nfV+y.sizeFeatE(2)+1:end)', [ y.nDimAttr(2) 1 ]) ]; 
    
    [ res ] = immatching_FMVC_wrapper( x, model.view_m, svm, y );    
    yhat = y;
    if nnz(res.X) == numel(y.seq)
        yhat.seq = res.seq;
        yhat.pt = res.frame(1:2,res.seq);
        yhat.frame = res.frame(:,res.seq);
    else
        rand_y_seg = randperm(size(res.pt,2));
        yhat.seq = rand_y_seg(1:numel(y.seq));
        yhat.pt = res.frame(1:2,yhat.seq);
        yhat.frame = res.frame(:,yhat.seq);
    end 

end

function yhat = constraintCB_sw(param, model_lib, x, y)
% function to find the most violating constraint which is compatible with w

% margin rescaling: argmax_y delta(yi, y) + <psi(x,y), w>
  load(y.modelPath, 'model');  
    svm = model.svm;
    svm.wV = svm.wV.*repmat(model_lib.w(1:y.sizeFeatV(1)),[1 y.sizeFeatV(2)]);
    svm.wE = svm.wE.*[ model_lib.w(y.sizeFeatV(1)+1)*ones([ y.nDimAttr(1) y.sizeFeatE(2) ]);...
    model_lib.w(y.sizeFeatV(1)+2)*ones([ y.nDimAttr(2) y.sizeFeatE(2) ]) ];
    [ res ] = immatching_FMVC_wrapper( x, model.view_m, svm, y );    
    yhat = y;
    if nnz(res.X) == numel(y.seq)
        yhat.seq = res.seq;
        yhat.pt = res.frame(1:2,res.seq);
        yhat.frame = res.frame(:,res.seq);
    else
        rand_y_seg = randperm(size(res.pt,2));
        yhat.seq = rand_y_seg(1:numel(y.seq));
        yhat.pt = res.frame(1:2,yhat.seq);
        yhat.frame = res.frame(:,yhat.seq);
    end 

end