function [P, numPairMatch] = IMGM_offline(rawMat,nodeCnt,graphCnt,param)
    global affinity
        
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

        increMatching = IMGM_old(simAP, matTmp, param);
        prevMatching = increMatching;
        numPairMatch = numPairMatch + param.N + graphStep;
    end
    
    P = increMatching;
end