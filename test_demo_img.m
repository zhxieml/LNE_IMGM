function [ accuracy error_EP score_GM ] = test_demo_img( iExp, varargin )

set_conf_img;
set_alg;

datatype = -1; % all as a default
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

% load dataset
load(getDatasetPath(conf), 'classes', 'dataset');

if strcmpi(conf.learningType, 'non')
    model = init_model_img(dataset{iExp},conf); % get the initial model
    model = init_svm_img( model ); 
elseif ~exist(getModelPath(conf, iExp))
    error('model file does not exist! %s!',getModelPath(conf, iExp));
else
    load(getModelPath(conf,iExp), 'model');
end

classid = find(strcmp(classes,model.class));

% initialize training and test examples for the given class
exmpAll = makeExampleSet_img(classid, dataset{iExp});
if datatype == 1
    exmp = exmpAll( find([ exmpAll.type ] == 1) ); % test
else
    exmp = exmpAll( find([ exmpAll.type ] == 0) ); % train
end

svm = model.svm;

%bShowImg = true;
bShowImg = false;
%% do matching
if bShowImg 
    hFigRes = figure;
end


error_EP = zeros(size(exmp(1).ptGT,2),numel(exmp));
accuracy = zeros(1,numel(exmp));
score_GM = zeros(1,numel(exmp));
for kk = 1:numel(exmp)
    x = exmp(kk);
    [ res cdata ] = immatching_wrapper( x, model.view_m, svm );
    error_EP(:,kk) = sqrt(sum((x.ptGT - res.pt).^2,1))/sqrt(sum(x.object_size.^2));
    accuracy(kk) = nnz(error_EP(:,kk) < conf.distRatioTrue)/numel(error_EP(:,kk)); 
    score_GM(kk) = res.score_GM;
    
    if nnz(res.X) ~= size(x.ptGT,2)
        err('# of matches: %d, but # of GT : %d', nnz(res.X), size(x.ptGT,2));
    end
    if bShowImg
        figure(hFigRes);
        set(hFigRes, 'Name', sprintf('node: %7.4f  edge: %7.4f  score: %7.4f   accuracy: %.2f', res.sV, res.sE, res.score_GM, accuracy(kk)));
        clf(hFigRes);
        showMatchingResult( cdata, res );
        pause;
    end
      
end

accuracy = mean(accuracy);
error_EP = mean(error_EP(:));
score_GM = mean(score_GM);

fprintf('** average on %d test examples\n',numel(exmp));
fprintf('- accuracy (dist ratio:%.4f): %.4f    endpoint error: %.4f\n', conf.distRatioTrue, accuracy, error_EP );





