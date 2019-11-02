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
    
    savePath = fullfile(affinityDir, cls, sprintf('in%02d_out%02d.mat', nInlier, nOutlier));
    if isfile(savePath)
        load(savePath, 'affinity');
        load(savePath, 'rawMat');
    else
        listOfFeatFile = dir(fullfile(featDir, cls, '*.mat'));
        M = getTrans(nodeCnt);
        affinity.K = cell(graphCnt, graphCnt);
        affinity.GT = cell(graphCnt, graphCnt);
        rawMat = cell(graphCnt, graphCnt);
        
        for idxGraph1 = 1:graphCnt
            viewInfo1 = load(fullfile(listOfFeatFile(idxGraph1).folder, listOfFeatFile(idxGraph1).name), 'view');
            rawMat{idxGraph1, idxGraph1} = eye(nodeCnt);
            affinity.GT{idxGraph1, idxGraph1} = eye(nodeCnt);
            affinity.pointFeat{idxGraph1} = viewInfo1.view.desc;
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
                gtFileName = fullfile(gtDir, cls, sprintf("%s_%s_%s.mat", cls, ...
                            listOfFeatFile(idxGraph1).name(end-7:end-4), listOfFeatFile(idxGraph2).name(end-7:end-4)));
                try
                    load(gtFileName, 'groundTruth');
                catch
                    fprintf("error:idxGraph1 = %d, idxGraph2 = %d\n", idxGraph1, idxGraph2);
                    pause;
                end
                Xgt = zeros(nodeCnt);
                Xgt(1:nInlier, 1:nInlier) = groundTruth.assign;
                affinity.GT{idxGraph1, idxGraph2} = Xgt;
                affinity.GT{idxGraph2, idxGraph1} = Xgt';
                X = reshape(res.X, nodeCnt, nodeCnt);
                rawMat{idxGraph1, idxGraph2} = X;
                rawMat{idxGraph2, idxGraph1} = X';
            end
        end
        if ~exist(fullfile(affinityDir, cls),'dir')
        mkdir(fullfile(affinityDir, cls));
        end
        save(savePath, 'rawMat', 'affinity');
    end
    % disrupt the order
    randOrder = randperm(graphCnt);
    affinity.K = affinity.K(randOrder, randOrder);
    affinity.GT = cell2mat(affinity.GT(randOrder, randOrder));
    rawMat = cell2mat(rawMat(randOrder, randOrder));
    affinity.pointFeat = affinity.pointFeat(1, randOrder);
    affinity.BiDir = 1;
    affinity.edgeAffinityWeight = 1;
    affinity.angleAffinityWeight = 0;
    affinity.nodeCnt = nodeCnt;
    affinity.graphCnt = graphCnt;
    affinity.EG = cell(graphCnt, 1);
    [r,c] = find(ones(nodeCnt, nodeCnt));
    for i = 1:graphCnt
        affinity.EG{i} = [r, c]';
        affinity.nP{i} = nodeCnt;
        affinity.nE{i} = nodeCnt * nodeCnt;
    end
end