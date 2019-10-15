function [ model ] = show_modelA_img( iExp, varargin )

%if nargin < 2 
typeModel = 1; %learned model
%end

set_conf_img;
set_param_GM;

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'conf'
      conf = arg;  
    %case 'datatype'
    %  datatype = arg; 
    case 'verbose'
      verbose = arg;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

load(getModelPath(conf,iExp), 'model');
strTitle = 'model graph after learning';
if typeModel == 0 
    model = init_svm_img(model);
    strTitle = 'model graph before learning';
end

% speficy graph type
frame_m = model.view_m.frame;
color_m = model.colorCode;
adjmat_m = model.graph_m.adjmat;

%listOfModelDataFile = dir(sprintf(VOCopts.imgpath, [ cls '/*' ]));
%fileName = listOfModelDataFile(1).name(1:end-4);
%filePathName = sprintf(VOCopts.imgpath, [ cls '/' listOfModelDataFile(1).name(1:end-4)] );

%anno_filePathName = sprintf(VOCopts.annopath,[ cls '/annotation_' listOfModelDataFile(1).name(end-7:end-4) ]);
%load(anno_filePathName, 'box_coord', 'obj_contour');%% load the annotated data

%svm = renormA(svm);
wV = sqrt(sum(model.svm.wV.^2,1)); 
wE = sqrt(sum(model.svm.wE.^2,1)); 
%wV = sum(abs(svm.wV),1); 
%wE = sum(abs(svm.wE),1); 

if 0
    % show node shapes (ellipses)
    for i = 1:size(frame_m,2)
        fprintf('%f - node  %d (feat %d) +:%d 0:%d -:%d\n',...
            wV(i),i,model.view_m.type(i), nnz(model.svm.wV(:,i)>0),...
            nnz(model.svm.wV(:,i)==0), nnz(model.svm.wV(:,i)<0));
    end
    [ r, c ] = find(adjmat_m);
    [ ie ] = find(adjmat_m(:));
    for k = 1:numel(r)
        %svm.aE(:,ie(k))
        fprintf('%f - edge (%d,%d)\n',wE(ie(k)),r(k),c(k));
        %pause;
    end
    %fprintf('\n%f - bias \n',svm.aBias);
end

% show the image
figure('Name',strTitle, 'NumberTitle', 'off'); 
imgToShow = rgb2gray(model.view_m.img);
imshow(appendimages( imgToShow, imgToShow )); hold on;

% show GT bbox
%box_handle = rectangle('position', [box_coord(3), box_coord(1), box_coord(4)-box_coord(3), box_coord(2)-box_coord(1)]);
%set(box_handle, 'edgecolor','y', 'linewidth',5);    
visualizeSIFT2((model.svm.wV), frame_m, 'descscale', fparam.descScale);

% visualize nodes with their weights
drawWeightedGraph( frame_m, [ size(model.view_m.img,2) 0 ], color_m, adjmat_m, wV, wE(:) );

hold off
end




