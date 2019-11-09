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
    % target.data.feat �??行是�??个点的feature
    % target.data.point [nodeCnt, 2]
    permutation = zeros(totalCnt, nodeCnt);
    for gc = 1:totalCnt
        permutation(gc, :) = randperm(nodeCnt);
    end
    for x = 1:totalCnt
        rawFeat = load(fullfile(listOfFeatFile(x).folder, listOfFeatFile(x).name));
        rawFeatPoint = rawFeat.frames_anno(1:2, 1:nodeCnt);
        rawFeatFeat = rawFeat.descriptors(:, 1:nodeCnt);
        target.data{x}.point = zeros(size(rawFeatPoint)); % [2, nodeCnt]
        target.data{x}.feat = zeros(size(rawFeatFeat)); % [128, nodeCnt]
        np = 1:nodeCnt;
        target.data{x}.point(:, permutation(x, np)) = rawFeatPoint(:, np);
        target.data{x}.feat(:, permutation(x, np)) = rawFeatFeat(:, np);
        target.data{x}.point = target.data{x}.point';
        target.data{x}.feat = target.data{x}.feat';
        nameXgt = listOfFeatFile(x).name;
        for y = x+1:totalCnt
            nameYgt = listOfFeatFile(y).name;
            gtPath = fullfile(gtDir, cls, sprintf("%s_%s_%s.mat", cls, nameXgt(end-7:end-4), nameYgt(end-7:end-4)));
            gt = load(gtPath);
            assert(size(gt.groundTruth.assign, 1) == nInlier, "error size(gt.groundTruth.assign, 1) ~= nInlier\n");
            Xgt = eye(nodeCnt);
            Xgt(1:nInlier, 1:nInlier) = gt.groundTruth.assign;
            perXgt = zeros(nodeCnt);
            [row, col]= find(Xgt);
            for m = 1:length(row)
                r = row(m);
                c = col(m);
                perXgt(permutation(x, r), permutation(y, c)) = 1;
            end
            target.GT{x, y} = perXgt;
            target.GT{y, x} = perXgt';
        end
    end
    
end

