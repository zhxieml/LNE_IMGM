function [P, numPairMatch] = ANC_IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, param)
    % batch version of ANC_IMGM
    % put multiple graphs at one time
    % calculate adaptive graph order to promote accuracy
    %for debug
    assert(baseGraphCnt >= 0);
    assert(batchSize >= 1, "error: batch size < 1\n");
    assert(size(rawMat, 1) == nodeCnt*(baseGraphCnt+batchSize), "error in IMGM_batch\n");

    param.n = nodeCnt;
    aptOrder = 1:batchSize;
    graphCnt = size(rawMat, 1) / nodeCnt;
    batchSize = min(batchSize, graphCnt - baseGraphCnt);
    numPairMatch = 0;

    if useAptOrder && batchSize > 1
        % calculate adaptive order 
        dimN = baseGraphCnt*nodeCnt + 1: (baseGraphCnt+batchSize)*nodeCnt;
        aptOrder = cal_adaptive_graph_order(rawMat(dimN, dimN),nodeCnt,batchSize);
        aptOrder = aptOrder + baseGraphCnt;
        aptOrder = [1:baseGraphCnt, aptOrder']; % length(aptOrder) = baseGraphCnt + batchSize
        % reorder rawMat, affinity and target
        rawMatReOrder = crop_rawMat(aptOrder, rawMat, nodeCnt);
        affinityReOrder = crop_affinity(aptOrder, affinity);
        targetReOrder = crop_target(aptOrder, target);
    else
        rawMatReOrder = rawMat;
        affinityReOrder = affinity;
        targetReOrder = target;
    end

    % perform batch increment by repetting single step ANC_IMGM
    % baseGraphCnt + 1 -> baseGraphCnt + batchSize;
    for bs = 1:batchSize
        param.N = baseGraphCnt + bs - 1;
        if bs == 1
            prevMatching = rawMatReOrder(1:baseGraphCnt*nodeCnt, 1:baseGraphCnt*nodeCnt);
        end
        increCnt = param.N + 1;
        matTmp = rawMatReOrder(1:nodeCnt*increCnt, 1:nodeCnt*increCnt);
        if ~isempty(prevMatching)
            matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching;
        end
        affScore = cal_pair_graph_inlier_score_local(affinityReOrder, matTmp, nodeCnt, increCnt, nodeCnt);
        [increMatching, numPairMatchSingle] = ANC_IMGM(affinityReOrder, affScore, rawMatReOrder, targetReOrder, param);
        % TBIMGM(globalVar, affScore, rawMat, target, param)
        numPairMatch = numPairMatch + numPairMatchSingle;
        prevMatching = increMatching;
    end

    % recover the order
    if useAptOrder
        [~, rcvrOrder] = sort(aptOrder, 'ascend');
        P = crop_rawMat(rcvrOrder, increMatching, nodeCnt);
    else
        P = increMatching;
    end
end