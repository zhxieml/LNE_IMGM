function do_makeGT( )
% function to make GTs 

close all; clc;
modeGTguide = true;

set_conf_img;
%% Params for test algorithms
set_param_GM;
set_alg;
cls = conf.class;

disp('== initialize an object model...');
listOfModelDataFile = dir(fullfile(conf.imgDir, cls, '*.png'));

%% load image1
idxTest1 = 1;
filePathName1 = fullfile(conf.imgDir, cls, listOfModelDataFile(idxTest1).name );
anno_filePathName1 = fullfile(conf.annoDir, cls, [ listOfModelDataFile(idxTest1).name(1:end-4) '.mat' ]);
view_m = extract_localfeatures_mcho( filePathName1, fparam, 'verbose', true );
pts_coord = []; 
if exist(anno_filePathName1)
    load(anno_filePathName1, 'pts_coord');
end
GT_M = pts_coord;
if modeGTguide && size(GT_M,2) > 0
    % use graph matching to guide GT
    idxSel = [];
    for i=1:size(GT_M,2)
        idx1 = nearestneighbour(GT_M(:,i), view_m.frame(1:2,:), 'NumberOfNeighbours', 1);
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
    model.classid = 1;
    model.graph_m = graph_m;
    model.view_m = view_m;
    %model.sourceset = sourceset;
    % color code for visualization
    model.colorCode = makeColorCode(size(view_m.frame,2));
    model = init_svm_img( model ); 
end


%% load image2
idxTest2 = [ 2:100 ];
bCont = 1;  i = 1;
%% main loop start
while bCont
        filePathName2 = fullfile(conf.imgDir, cls, listOfModelDataFile(idxTest2(i)).name );
        anno_filePathName2 = fullfile(conf.annoDir, cls, [ listOfModelDataFile(idxTest2(i)).name(1:end-4) '.mat' ]);
        %load(anno_filePathName2, 'box_coord', 'obj_contour');%% load the annotated data
        %bbox_in2   = [ box_coord(3), box_coord(1), box_coord(4), box_coord(2)];
        % extract features
        %viewInfo2 = extract_localfeatures_mcho( filePathName2, fparam, 'bbox', bbox_in2, 'verbose', true );
        viewInfo2 = extract_localfeatures_mcho( filePathName2, fparam, 'verbose', true );

        %% make it into cdata
        cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;
        cdata.view(1) = view_m;  cdata.view(2) = viewInfo2;

        % visualize it
        if 0
            imgInput = appendimages( cdata.view(1).img, cdata.view(2).img );
            imgInput = double(imgInput)./255;
            figure('Name','model initialization', 'NumberTitle', 'off'); 
            imshow(rgb2gray(imgInput)); hold on;

            offset_x = size(cdata.view(1).img,2);
            visualizeFeatures(cdata.view(1).frame, 'style', 'point', 'colorcode', 'g');
            visualizeFeatures(cdata.view(2).frame, 'style', 'point', 'colorcode', 'g', 'offset', [offset_x 0] );

            box_handle1 = rectangle('position', bbox2rect(bbox_in1) );
            set(box_handle1, 'edgecolor','y', 'linewidth',5);
            bbox_in2 = bbox_in2 + [ 1 0 1 0 ]*size(cdata.view(1).img, 2);
            box_handle2 = rectangle('position', bbox2rect(bbox_in2) );
            set(box_handle2, 'edgecolor','y', 'linewidth',5);
        end
        
        %% perform GT labeling
        disp ([ filePathName2 ' file loaded for labeing GT.']);
        
        if i == 0
            visualizeSIFT2(cdata.view(1).desc, cdata.view(1).frame, 'descscale', fparam.descScale);
        end
    
        pts_coord = []; 
        if exist(anno_filePathName1)
            load(anno_filePathName1, 'pts_coord');
        end
        GT_M = pts_coord;
        pts_coord = []; 
        if exist(anno_filePathName2)
            load(anno_filePathName2, 'pts_coord');
        end
        GT_T = pts_coord;
        % GT guise mode using graph matching
        if modeGTguide && size(GT_M,2) > 0 && size(GT_T,2) == 0
            x.fname = filePathName2;
            x.tname = [];
            [ res ] = immatching_wrapper( x, model.view_m, model.svm );
            GT_T = res.pt;
        end
        % make ground truth for evaluation
        [ modeEnd GT_M GT_T ] = make_groundTruth_feat(cdata, GT_M, GT_T);
        %A = cdata.view(1).feat(GT(:,1),1:5); B = cdata.view(2).feat(GT(:,2),1:5);
%         save(['./data_matches/GT_' listOfModelDataFile.name(1:end-4) '+' testDataFile(1:end-4) '.mat'], 'A', 'B');
        if modeEnd > 0 
            if size(GT_M,2) > 0
                pts_coord = GT_M;
                if exist(anno_filePathName1)
                    save(anno_filePathName1, 'pts_coord','-append');
                else
                    save(anno_filePathName1, 'pts_coord');
                end
            end
            if size(GT_T,2) > 0
                pts_coord = GT_T;
                if exist(anno_filePathName2)
                    save(anno_filePathName2, 'pts_coord','-append');
                else
                    save(anno_filePathName2, 'pts_coord');
                end
                disp ([ filePathName2 ' file saved for ground truth information.']);
            end
            
        end
        if modeEnd == 1 % next
            i = i + 1;
        elseif modeEnd == 2 % previous
            i = i - 1; 
        elseif modeEnd == 3 % terminate
            bCont = 0;
        end
       
end

