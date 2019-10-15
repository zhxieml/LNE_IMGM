function [ model ] = init_model_img( dataset, conf )
% initialize a model graph

%load(getDatasetPath(conf), 'classes', 'dataset');
%classid = find(strcmp(classes,cls));
%disp( [ '== initialize an object model of ' cls ] );
set_param_GM;
%% construct graphical models
% nodes 
%id_img = find( [ dataset.classid ] == classid && [ dataset.type ] == 0 , 1,'first' ); % first training image
id_cnd = find( [ dataset.type ] == 0 ); % training images
tmp = 1;%randperm(numel(id_cnd)); 
id_img = id_cnd(tmp(1));% randomly select one example
classid = dataset(id_img).classid;
view_m = extract_localfeatures_wrapper( dataset(id_img).imgpath, dataset(id_img).tmppath, fparam );

idxSel = [];
pt = dataset(id_img).ptGT;
for i=1:size(pt,2)
    idx1 = nearestneighbour(pt(:,i), view_m.frame(1:2,:), 'NumberOfNeighbours', 1);
    idxSel = [idxSel idx1];
end

view_m.frame = view_m.frame(:, idxSel);
view_m.type = view_m.type(:, idxSel);
view_m.desc = view_m.desc(:, idxSel);

graph_m.nodes = (1:size(view_m.frame,2))';
% adjacency matrix 
graph_m.adjmat = makeAdjMatOfCompleteGraph(graph_m.nodes);
graph_m.attrE = computeEdgeAttr(view_m.frame);
%graph_m.x = sourceset(classid).x;

% save the initial model
model.class = conf.class;
model.classid = classid;
model.graph_m = graph_m;
model.view_m = view_m;
%model.sourceset = sourceset;
% color code for visualization
model.colorCode = makeColorCode(size(view_m.frame,2));
drawnow;
