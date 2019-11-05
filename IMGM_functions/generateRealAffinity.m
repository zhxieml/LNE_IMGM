function [affinity, rawMat] = generateRealAffinity
    % generate affinity matrix
    global target;
    set_param_GM;
    affinityDir = target.config.affinityDir;
    nInlier = target.config.nInlier;
    nOutlier = target.config.nOutlier;
    featDir = target.config.featDir;
    gtDir = target.config.gtDir;
    cls = target.config.class;
    nodeCnt = nOutlier + nInlier;
    graphCnt = target.config.graphCnt;
    affinityRatio = target.config.affinityRatio;
    
<<<<<<< HEAD
    init_path;
    set_conf_img;
    set_param_GM;
    cls = conf.class;
    
    savePath = fullfile(conf.affinityDir, conf.class, sprintf('in%02d_out%02d.mat', conf.numInlier, conf.numOutlier));
    if isfile(savePath)
        load(savePath, 'affinity');
        load(savePath, 'rawMat');
        return;
    end
    listOfFeatFile = dir(fullfile(conf.featDir, cls, '*.mat')); % at most 10 outliers
    graphCnt = length(listOfFeatFile);
    nodeCnt = conf.numInlier + conf.numOutlier;
    M = getTrans(nodeCnt);
    affinity.K = cell(graphCnt, graphCnt);
    affinity.GT = cell(graphCnt, graphCnt);
    rawMat = cell(graphCnt, graphCnt);
    
    for idxGraph1 = 1:graphCnt
        viewInfo1 = load(fullfile(listOfFeatFile(idxGraph1).folder, listOfFeatFile(idxGraph1).name), 'view');
        rawMat{idxGraph1, idxGraph1} = eye(nodeCnt);
        affinity.GT{idxGraph1, idxGraph1} = eye(nodeCnt);
        for idxGraph2 = idxGraph1+1:graphCnt
            fprintf("matching %s to %s...\n", listOfFeatFile(idxGraph1).name, listOfFeatFile(idxGraph2).name);
            viewInfo2 = load(fullfile(listOfFeatFile(idxGraph2).folder, listOfFeatFile(idxGraph2).name), 'view');
            cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;
            cdata.view(1) = viewInfo1.view;  cdata.view(2) = viewInfo2.view;
            res = wrapper_matchingBasic(cdata);
            % affinity matrix
            affinity.K{idxGraph1, idxGraph2} = res.affinityMatrix;
            affinity.K{idxGraph2, idxGraph1} = M*(res.affinityMatrix*M);
            % ground_truth
            gtFileName = fullfile(conf.gtDir, cls, sprintf("%s_%s_%s.mat", cls, ...
                        listOfFeatFile(idxGraph1).name(end-7:end-4), listOfFeatFile(idxGraph2).name(end-7:end-4)));
            try
                load(gtFileName, 'groundTruth');
            catch
                fprintf("error:idxGraph1 = %d, idxGraph2 = %d\n", idxGraph1, idxGraph2);
                pause;
            end
                
            Xgt = groundTruth.assign;
            
            if conf.numInlier > 0
                Xgt(conf.numInlier+1:nodeCnt, conf.numInlier+1:nodeCnt) = 0;
            end
            
            affinity.GT{idxGraph1, idxGraph2} = Xgt;
            affinity.GT{idxGraph2, idxGraph1} = Xgt';
            X = reshape(res.X, nodeCnt, nodeCnt);
            rawMat{idxGraph1, idxGraph2} = X;
            rawMat{idxGraph2, idxGraph1} = X';
=======
    savePath = fullfile(affinityDir, cls, sprintf('in%02d_out%02d.mat', nInlier, nOutlier));
    if isfile(savePath)
        load(savePath, 'affinity');
    else
        listOfFeatFile = dir(fullfile(featDir, cls, '*.mat'));
        for gc = 1:graphCnt
            feat = load(fullfile(listOfFeatFile(gc).folder, listOfFeatFile(gc).name);
            pointSet = feat.frames(1:2, :);
            [X, Y] = meshgrid(1:nodeCnt, 1:nodeCnt);
            LL = [X(:), Y(:)];
            edges = pointSet(LL(:,1),:)-pointSet(LL(:,2),:); %edge# \times 2
            edgeLen = sqrt(edges(:, 1).^2 + edges(:, 2).^2); %edge# \times 1
            mask = (edgeLen ~= 0);
            edgeAngle = zeros(size(edgeLen));
            edgeAngle(mask) = acos(edges(mask,1)./edgeLen(mask));

            affinity.edgeLen{gc} = edgeLen'/max(edgeLen);
            affinity.edgeAngle{gc} = edgeAngle';
            affinity.pointFeat{gc} = feat.desc;
        end
        for xview = 1:graphCnt
            if affinity.BiDir
                yviewSet = [1:xview-1, xview+1:graphCnt];
            else
                yviewSet = xview+1:graphCnt;
            end
            for yview = yviewSet
                if bUnaryEnable
                    affinity.KP{xview,yview} = exp(-conDst(affinity.pointFeat{xview}, affinity.pointFeat{yview},0) / Sacle_2D);
                else
                    affinity.KP{xview,yview} = zeros(nodeCnt,nodeCnt);
                end
                if bEdgeEnable
                    KQ = exp(-conDst(affinity.edgeLen{xview}, affinity.edgeLen{yview},0) / Sacle_2D);
                else
                    KQ = zeros(nodeCnt^2,nodeCnt^2);
                end
                if bAngleEnable
                    
                end
                KP = affinity.KP{xview, yview}(LL(:,1))
                affinity.K{xview, yview} = affinityRatio* KQ + (1-affinityRatio)*Kangle;
            end
>>>>>>> cd249f82b8ab6061fd7653eae1cdc4b4f6b775cf
        end
    end
    
end