function load_target_data
    global target

    switch target.config.database
    case "willow"
        target.config.dataDir = 'C:\Users\xzh87\OneDrive\Reference\Lab\code\WILLOW-ObjectClass_dataset';
        target.config.imgDir=[target.config.dataDir '/WILLOW-ObjectClass'];
        target.config.annoDir=[target.config.dataDir '/WILLOW-ObjectClass'];
        target.config.gtDir = [target.config.dataDir '/ground_truth'];
        target.config.resDir = './res';
        target.config.tmpDir = './tmp';
        target.config.class = 'Duck';%'Car','Winebottle','Duck'
        target.config.category = 'outlier';%'deform','outlier'
        switch target.config.category
            case 'deform' % same setting with 5th row in Table 1 in the PAMI paper 
                target.config.nInlier = 10;
                target.config.nOutlier = 0;
                target.config.featDir = [target.config.dataDir '\feature'];
                target.config.complete = 1;
                target.config.inCntType = 'all';% set 'all' for "only a few outlier case", e.g. Fig.1&2&3&4
            case 'outlier'
                target.config.nInlier = 10;
                target.config.nOutlier = 4;
                target.config.featDir = [target.config.dataDir '\feature'];
                target.config.complete = 1;
                target.config.inCntType = 'all';% set 'all' for "only a few outlier case", e.g. Fig.1&2&3&4
        end
        
        switch target.config.class
            case 'Duck'
                target.config.graphMaxCnt = min(50, target.config.graphMaxCnt);
            case 'Motorbike'
                target.config.graphMaxCnt = min(40, target.config.graphMaxCnt);
            case 'Car'
                target.config.graphMaxCnt = min(40, target.config.graphMaxCnt);
            case 'Face'
                target.config.graphMaxCnt = min(109, target.config.graphMaxCnt);
            case 'Winebottle'
                target.config.graphMaxCnt = min(66, target.config.graphMaxCnt);
        end

        target.config.totalCnt = target.config.graphMaxCnt;
        totalCnt = target.config.totalCnt;
        featDir = target.config.featDir;
        cls = target.config.class;
        nInlier = target.config.nInlier;
        nOutlier = target.config.nOutlier;
        nodeCnt = nInlier + nOutlier;
        bUnaryEnable = target.config.bUnaryEnable;
        bEdgeEnable = target.config.bEdgeEnable;
        target.GT = cell(totalCnt, totalCnt);
        target.data = cell(totalCnt, 1);
        listOfFeatFile = dir(fullfile(featDir, cls, '*.mat'));
        permutation = zeros(totalCnt, nInlier);

        for x = 1:totalCnt
            rawFeat = load(fullfile(listOfFeatFile(x).folder, listOfFeatFile(x).name));
            np = 1:nInlier;
            if bEdgeEnable
                rawFeatPoint = rawFeat.frames(1:2, 1:nodeCnt);
                target.data{x}.point = rawFeatPoint; % [2, nodeCnt]
                target.data{x}.point = target.data{x}.point';
            end
            if bUnaryEnable
                rawFeatFeat = rawFeat.descriptors(:, 1:nodeCnt);
                target.data{x}.feat = rawFeatFeat; % [128, nodeCnt]
                target.data{x}.feat = target.data{x}.feat';
            end
        end
    
    case "synthetic"
        target.config.bGraphMatch = 1;
        target.config.category = 'deform';%'deform','outlier','density','complete'
        target.config.inCntType = "all";
        switch target.config.category
        case 'deform' % same setting with 5th row in Table 1 in the PAMI paper 
            target.config.nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.15;
            target.config.density = .9;target.config.complete = 1;
        case 'outlier' % same setting with 6th row in Table 1 in the PAMI paper 
            target.config.nInlier = 6;target.config.nOutlier = 4;target.config.deform = 0;
            target.config.density = 1;target.config.complete = 1;
        case 'density' % same setting with 7th row in Table 1 in the PAMI paper 
            target.config.nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.0;
            target.config.density = 0.5;target.config.complete = 1;
        case 'complete' % same setting with 8th row in Table 1 in the PAMI paper 
            target.config.nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.05;
            target.config.density = 1;target.config.complete = 0.1;     
        end
    case "CMU-sequence"
        error("not finished!\n");
    end
end

