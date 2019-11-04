function load_real_data
    global target affinity
    set_param_GM;
    affinity.BiDir = target.config.affinityBiDir;
    affinity.edgeAffinityWeight = target.config.edgeAffinityWeight;
    affinity.angleAffinityWeight = target.config.angleAffinityWeight;
    nInlier = target.config.nInlier;
    nOutlier = target.config.nOutlier;
    featDir = target.config.featDir;
    gtDir = target.config.gtDir;
    cls = target.config.class;
    nodeCnt = nOutlier + nInlier;
    graphCnt = target.config.graphCnt;
    affinity.graphCnt = graphCnt;
    affinity.nodeCnt = nodeCnt;

    affinity.GT = zeros(nodeCnt*graphCnt);
    target.data = cell(graphCnt, 1);
    listOfFeatFile = dir(fullfile(featDir, cls, '*.mat'));
    % target.data.feat 一行是一个点的feature
    % target.data.point [nodeCnt, 2]
    for x = 1:graphCnt
        rawFeat = load(fullfile(listOfFeatFile(x).folder, listOfFeatFile(x).name));
        target.data{x}.point = rawFeat.frames(1:2,:)';
        target.feat{x}.feat = rawFeat.descriptor';
        nameXgt = listOfFeatFile(x).name;
        xview_gt = (x-1)*nodeCnt+1:(x-1)*nodeCnt+nInlier;
        for y = x+1:graphCnt
            nameYgt = listOfFeatFile(y).name;
            yview_gt = (y-1)*nodeCnt+1:(y-1)*nodeCnt+nInlier;
            gtPath = fullfile(gtDir, cls, sprintf("%s_%s_%s.mat", cls, nameXgt(end-7:end-4), nameYgt(end-7:end-4)));
            gt = load(gtPath);
            affinity.GT(xview_gt, yview_gt) = gt.groundTruth.assign;
        end
    end
    
end

