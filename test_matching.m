function test_matching
% test function for matching between two images

close all;
set_conf_img;
%% Params for test algorithms
set_param_GM;
set_alg;

cls = conf.class;

disp('== initialize an object model...');
listOfModelDataFile = dir(fullfile(conf.imgDir, cls, '*.png'));


hAlg1 = figure('Name',sprintf('%s', alg(1).strName), 'NumberTitle', 'off'); 
%hAlg2 = figure('Name',sprintf('%s', alg(2).strName), 'NumberTitle', 'off'); 
hFeat = figure('Name','model initialization', 'NumberTitle', 'off'); 

%% load image1
idxTest1 = 1;
filePathName1 = fullfile(conf.imgDir, cls, listOfModelDataFile(idxTest1).name );
if 0
    anno_filePathName1 = fullfile(conf.annoDir, cls, [ '/annotation_' listOfModelDataFile(idxTest1).name(end-7:end-4) '.mat' ]);
    load(anno_filePathName1, 'box_coord', 'obj_contour','pts_coord');%% load the annotated data    
    bbox_in1   = [ box_coord(3), box_coord(1), box_coord(4), box_coord(2)]; % x1 y1 x2 y2
else
    anno_filePathName1 = fullfile(conf.annoDir, cls, [ listOfModelDataFile(idxTest1).name(1:end-4) '.mat' ]);
    load(anno_filePathName1, 'pts_coord');%% load the annotated data    
    
end
% extract features
%viewInfo1 = extract_localfeatures_mcho( filePathName1, fparam, 'bbox', bbox_in1, 'minarearatio', 0.03, 'verbose', true );
viewInfo1 = extract_localfeatures_mcho( filePathName1, fparam, 'verbose', true );
if 1
    idx_sel1 = [];
    for i=1:size(pts_coord,2)
        pt1 = pts_coord(1:2,i);
        %idx1 = nearestneighbour(pt1, viewInfo1.frame(1:2,:), 'Radius', distClose);
        idx1 = nearestneighbour(pt1, viewInfo1.frame(1:2,:), 'NumberOfNeighbours', 1);
        idx_sel1 = [ idx_sel1 idx1 ];
    end
    %[ x_sel1 idx_sel1 ] = sel_pts_by_meanshift( viewInfo1.frame(1:2,:), 25);
    %idx_sel1
    viewInfo1.type = viewInfo1.type( idx_sel1 );
    viewInfo1.frame = viewInfo1.frame( :, idx_sel1);
    viewInfo1.desc = viewInfo1.desc( :, idx_sel1);
end

%% load image2
idxTest2 = 1:10; 

for i=1:numel(idxTest2)

    filePathName2 = fullfile(conf.imgDir, cls, listOfModelDataFile(idxTest2(i)).name );
    

    if 0
        anno_filePathName2 = fullfile(conf.annoDir, cls, [ '/annotation_' listOfModelDataFile(idxTest2(i)).name(end-7:end-4) '.mat' ]);
        load(anno_filePathName2, 'box_coord', 'obj_contour');%% load the annotated data
        bbox_in2   = [ box_coord(3), box_coord(1), box_coord(4), box_coord(2)];
    else
        anno_filePathName2 = fullfile(conf.annoDir, cls, [ listOfModelDataFile(idxTest2(i)).name(1:end-4) '.mat' ]);
        load(anno_filePathName2, 'pts_coord');%% load the annotated data
    end
    % extract features
    %viewInfo2 = extract_localfeatures_mcho( filePathName2, fparam, 'bbox', bbox_in2, 'verbose', true );
    viewInfo2 = extract_localfeatures_mcho( filePathName2, fparam, 'verbose', true );

    %% make it into cdata
    cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;
    cdata.view(1) = viewInfo1;  cdata.view(2) = viewInfo2;

    % visualize it
    imgInput = appendimages( cdata.view(1).img, cdata.view(2).img );
    imgInput = double(imgInput)./255;
    offset_x = size(cdata.view(1).img,2);
    set(0, 'CurrentFigure', hFeat); clf;
    %figure('Name','model initialization', 'NumberTitle', 'off'); 
    imshow(rgb2gray(imgInput)); hold on;
    visualizeFeatures(cdata.view(1).frame, 'style', 'point', 'colorcode', 'g');
    visualizeFeatures(cdata.view(2).frame, 'style', 'point', 'colorcode', 'g', 'offset', [offset_x 0] );

    %% perform matching
    res = wrapper_matchingBasic( alg(1), cdata );
%     fprintf('%10s - score:%.5f (%.5f), %.3fsec \n', alg(1).strName, res.score_raw, res.score_GM, res.time);   
    fprintf('%10s - score:%.5f (%.5f), xxxsec \n', alg(1).strName, res.score_raw, res.score_GM);      

    
    set(0, 'CurrentFigure', hAlg1); clf;
    %figure('Name',sprintf('%s', alg(1).strName), 'NumberTitle', 'off'); 
    imshow(rgb2gray(imgInput)); hold on;
    colorCode = makeColorCode(nnz(res.X));

    visualizeMatches( cdata.view(1).frame, cdata.view(2).frame, res.cand_matchlist(:,find(res.X)), ...
        'offset', [size(cdata.view(1).img,2) 0], 'style', 'frame', 'colorcode', colorCode, 'weight', res.X_raw(find(res.X)) ); 

    %if i == 0
    %    visualizeSIFT2(cdata.view(1).desc, cdata.view(1).frame, 'descscale', fparam.descScale);
    %end
    pause;
end



