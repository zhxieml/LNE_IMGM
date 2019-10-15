function [ classes dataset ] = init_dataset_img( conf )
% make img dataset

cls = conf.class;
classes{1} = cls;

warning off;
ims_all = dir(fullfile(conf.imgDir, cls, '*.png'))' ;
%ims = vl_colsubset(ims, conf.numTrain + conf.numTest) ;
%idendify files with match annotations
bValid = zeros(1,length(ims_all));
bFirst = true;
for ck = 1:length(ims_all)
    annopath = fullfile(conf.annoDir, cls, [ ims_all(ck).name(1:end-4) '.mat' ] );
    if exist(annopath)
        pts_coord = [];
        load(annopath, 'pts_coord');%% load the annotated data    
        if ~isempty(pts_coord)
            if bFirst
                nPt = size(pts_coord,2);
                bValid(ck) = true;
                bFirst = false;
            elseif nPt == size(pts_coord,2) && ~any(isnan(pts_coord(:)))
                bValid(ck) = true;
            end
        end            
    end
end
ims = ims_all(find(bValid));

% random splits for multiple trials
for iExp = 1:conf.nExp
    
    perm = randperm(length(ims));
    nSample = min(length(ims), conf.numTrain+conf.numTest+conf.numValid);
    sel  = perm(1:nSample);
    selType = zeros(1,nSample); % image type 0:train
    selType(conf.numTrain+1:min(conf.numTrain+conf.numTest, nSample)) = 1; % 1:test
    selType(conf.numTrain+conf.numTest+1:nSample) = 2; % 2:validation
    [ sel, ii_sel] = sort(sel);
    selType = selType(ii_sel);

    fprintf('Class %20s: %d train, %d test, %d validation images\n',conf.class,...
      nnz(selType==0), nnz(selType==1), nnz(selType==2) );
  
    dataset1 = [];
    for cj = 1:nSample 
      dataset1(cj).imgpath = fullfile(conf.imgDir, cls, ims(sel(cj)).name);
      dataset1(cj).annopath = fullfile(conf.annoDir, cls, [ ims(sel(cj)).name(1:end-4) '.mat' ] );
      dataset1(cj).tmppath   = fullfile(conf.tmpDir, [ cls '_' ims(sel(cj)).name(1:end-4) ] );
      
      if exist(dataset1(cj).annopath)
        %load(dataset(cntImage).annopath, 'box_coord', 'obj_contour','pts_coord');%% load the annotated data    
        load(dataset1(cj).annopath, 'pts_coord');%% load the annotated data    
        bbox   = [];%[ box_coord(3), box_coord(1), box_coord(4), box_coord(2)]; % x1 y1 x2 y2
        sizeBB = 0;%(bbox(3)-bbox(1)+1)*(bbox(4)-bbox(2)+1);
      else
        fprintf(' - no annotation of %s\n',dataset1(cj).imgpath);  
        bbox   = [ ];
        sizeBB = 0;
      end
      
      dataset1(cj).type  = selType(cj); % 0:train 1:test 2:validation
      dataset1(cj).classid  = 1;         % class id
      dataset1(cj).bbox   = bbox;
      dataset1(cj).flip    = false;
      dataset1(cj).trunc   = 0;%rec.objects(j).truncated;
      dataset1(cj).dataid = cj;
      dataset1(cj).size   = sizeBB;
      dataset1(cj).ptGT = pts_coord;        
    end
    
    dataset{iExp} = dataset1;  
end

warning on;
