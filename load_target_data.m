function load_target_data
    global target
    set_param_GM;
    nInlier = target.config.nInlier;
    nOutlier = target.config.nOutlier;
    nodeCnt = nInlier + nOutlier;
    featDir = target.config.featDir;
    gtDir = target.config.gtDir;
    cls = target.config.class;
    totalCnt = target.config.totalCnt;

    target.GT = cell(totalCnt, totalCnt);
    target.data = cell(totalCnt, 1);
    listOfFeatFile = dir(fullfile(featDir, cls, '*.mat'));
    % target.data.feat �?行是�?个点的feature
    % target.data.point [nodeCnt, 2]
    for x = 1:totalCnt
        rawFeat = load(fullfile(listOfFeatFile(x).folder, listOfFeatFile(x).name));
        target.data{x}.point = rawFeat.frames(1:2,1:nodeCnt)';
        target.data{x}.feat = double(rawFeat.descriptors(:, 1:nodeCnt)');
        nameXgt = listOfFeatFile(x).name;
        for y = x+1:totalCnt
            nameYgt = listOfFeatFile(y).name;
            gtPath = fullfile(gtDir, cls, sprintf("%s_%s_%s.mat", cls, nameXgt(end-7:end-4), nameYgt(end-7:end-4)));
            gt = load(gtPath);
            assert(size(gt.groundTruth.assign, 1) == nInlier, "error size(gt.groundTruth.assign, 1) ~= nInlier\n");
            target.GT{x, y} = gt.groundTruth.assign;
            target.GT{y, x} = gt.groundTruth.assign';
        end
    end
    
end

