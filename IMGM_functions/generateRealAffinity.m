function [affinity, rawMat] = generateRealAffinity
    % generate affinity matrix
    global target;
    
    savePath = fullfile(conf.affinityDir, cls, sprintf('in%02d_out%02d.mat', conf.numInlier, conf.numOutlier));
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
        affinity.GT{idxGraph1, idxGraph1} = eye(conf.numInlier);
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
            affinity.GT{idxGraph1, idxGraph2} = Xgt;
            affinity.GT{idxGraph2, idxGraph1} = Xgt';
            X = reshape(res.X, nodeCnt, nodeCnt);
            rawMat{idxGraph1, idxGraph2} = X;
            rawMat{idxGraph2, idxGraph1} = X';
        end
    end
    rawMat = cell2mat(rawMat);
    affinity.GT = cell2mat(affinity.GT);
    if ~exist(fullfile(conf.affinityDir, cls),'dir')
       mkdir(fullfile(conf.affinityDir, cls));
    end
    save(savePath, 'rawMat', 'affinity');
end