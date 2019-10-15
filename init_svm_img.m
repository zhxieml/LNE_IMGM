function [ model ] = init_svm_img( model, varargin )
% initialize a svm based on a given model
% learning options
svm = struct('lambda',0.0,'eta0',0.0,'mu0',1.0,'tstart',0.0,...
    'wV',[],'wE',[],'wDivisor',1.0,'wBias',0,...
    'aV',[],'aE',[],'aDivisor',1.0,'wFraction',0.0,'aBias',0,'t',0);

svm.graph_type = 'full';

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'verbose'
      verbose = arg;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

% Weight masks on the extracted features
svm.init_maskV = 1.0;%0.2; % for SIFT %0.3
svm.init_maskE = [ 1.0 1.0 ];%[ 0.02 0.02 ]; % for rho, theta

if ~isempty(model)
    svm.wV = zeros(size(model.view_m.desc),'single');
    for j=1:size(model.view_m.frame,2)
        desc1 = model.view_m.desc(:,j);
        desc1 = desc1 / sqrt(sum(desc1.^2));    
        svm.wV(:,j) = svm.init_maskV * desc1;
    end

    [ featBin rho_bins theta_bins ] = makeFeatBin(model.graph_m.attrE, 5, 5);
    svm.wE = svm.init_maskE(1) * featBin;
    svm.rho_bins = rho_bins;
    svm.theta_bins = theta_bins;
else
    svm.wV = [];
    svm.wE = [];
    svm.rho_bins = 0;
    svm.theta_bins = 0;
end

model.svm = svm;

