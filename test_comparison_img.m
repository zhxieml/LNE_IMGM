function test_comparison_img( iExp, varargin )
% test function for matching between two images

set_conf_img;
set_alg;

cls = conf.class;

datatype = 1; % test as a default
verbose = true;

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'conf'
      conf = arg;   
    case 'datatype'
      datatype = arg; 
    case 'verbose'
      verbose = arg;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

set_conf_img;
% load dataset
load(getDatasetPath(conf), 'classes', 'dataset');
classid = find(strcmp(classes,cls));


model_i = init_model_img(dataset{iExp},conf); % get the initial model
model_i = init_svm_img( model_i ); % a standard inital model

if ~exist(getModelPath(conf,iExp))
    error('model file does not exist! %s!',getModelPath(conf,iExp));
else
    load(getModelPath(conf, iExp), 'model'); % after learning
end

svm = model.svm;
svm_i = model_i.svm;

% initialize training and test examples for the given class
exmpAll = makeExampleSet_img(classid, dataset{iExp});
if datatype == 1
    exmp = exmpAll( find([ exmpAll.type ] == 1) ); % test
else
    exmp = exmpAll( find([ exmpAll.type ] == 0) ); % train
end

bShowImg = true;
%% do matching
if bShowImg 
    hFigRes=figure; hFigRes_i=figure;
end

%disThreshold = 20;
endpointErr = zeros(size(exmp(1).ptGT,2),numel(exmp));
accuracy = zeros(1,numel(exmp));
endpointErr_i = zeros(size(exmp(1).ptGT,2),numel(exmp));
accuracy_i = zeros(1,numel(exmp));
for kk = 1:numel(exmp)
    x = exmp(kk);
    % for leared model
    [ res cdata ] = immatching_wrapper( x, model.view_m, svm );
    endpointErr(:,kk) = sqrt(sum((x.ptGT - res.pt).^2,1))/sqrt(sum(x.object_size.^2));
    accuracy(kk) = nnz(endpointErr(:,kk) < conf.distRatioTrue)/numel(endpointErr(:,kk)); 
    
    % for a standard initial model
    [ res_i cdata_i ] = immatching_wrapper( x, model.view_m, svm_i );
    endpointErr_i(:,kk) = sqrt(sum((x.ptGT - res_i.pt).^2,1))/sqrt(sum(x.object_size.^2));
    accuracy_i(kk) = nnz(endpointErr_i(:,kk) < conf.distRatioTrue)/numel(endpointErr_i(:,kk)); 
    
    % make the indicator function for visualization
    res.indicator = endpointErr(:,kk) < conf.distRatioTrue;
    res_i.indicator = endpointErr_i(:,kk) < conf.distRatioTrue;
    
    if nnz(res.X) ~= size(x.ptGT,2)
        err('# of matches: %d, but # of GT : %d', nnz(res.X), size(x.ptGT,2));
    end
    if bShowImg
        figure(hFigRes);
        set(hFigRes, 'Name', sprintf('HARG learning -- node: %7.4f  edge: %7.4f  score: %7.4f   accuracy: %.2f', res.sV, res.sE, res.score_GM, accuracy(kk)));
        clf(hFigRes);
        showMatchingResult( cdata, res, model.colorCode );
        
        figure(hFigRes_i);
        set(hFigRes_i, 'Name', sprintf('w/o learning -- node: %7.4f  edge: %7.4f  score: %7.4f   accuracy: %.2f', res_i.sV, res_i.sE, res_i.score_GM, accuracy_i(kk)));
        clf(hFigRes_i);
        showMatchingResult( cdata_i, res_i, model.colorCode );
        pause;
    end
      
end

fprintf('** average on %d test examples\n',numel(exmp));
fprintf('- learned:  accuracy (dist ratio:%.4f): %.4f    endpoint error: %.4f\n', conf.distRatioTrue, mean(accuracy), mean(endpointErr(:)) );
fprintf('- standard: accuracy (dist ratio:%.4f): %.4f    endpoint error: %.4f\n', conf.distRatioTrue, mean(accuracy_i), mean(endpointErr_i(:)) );






