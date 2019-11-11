function [P, numPairMatch] = TBIMGM_offline(rawMat,nodeCnt,graphCnt,param)
    global affinity
    
    if param.adaptive_order
        rawMat = cal_adaptive_graph_order(rawMat, nodeCnt, graphCnt);
    end
    
    baseGraphCnt = 2;
    paraCnt = graphCnt - baseGraphCnt;
    graphStep = 1;
    
    numPairMatch = 0;
    
    % # outliers = 0
%     scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,nodeCnt);
    
    prevMatching = rawMat(1:nodeCnt*baseGraphCnt, 1:nodeCnt*baseGraphCnt);
    
    param.n = nodeCnt;
    param.graphStep = graphStep;

    for parak = 1:paraCnt
        param.N = baseGraphCnt + (parak-1)*graphStep;
        
        matTmp = rawMat(1:nodeCnt*(param.N+graphStep), 1:nodeCnt*(param.N+graphStep));
        matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N) = prevMatching;
        
        scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp,affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
        
        simAP = scrDenomMatInCntTmp;
        param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));

        [increMatching, increNumPairMatching] = TBIMGM(affinity, simAP, matTmp, param);
        prevMatching = increMatching;
        numPairMatch = numPairMatch + increNumPairMatching;
    end
    
    P = increMatching;
end