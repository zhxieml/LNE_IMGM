function P = TBIMGM_offline(rawMat,nodeCnt,graphCnt,param)
    global affinity
    
    orderedRawMat = cal_adaptive_graph_order(rawMat, nodeCnt, graphCnt);
    
    baseGraphCnt = 20;
    paraCnt = graphCnt - baseGraphCnt;
    graphStep = 1;
    sigma = 0;
    
    % # outliers = 0
    scrDenomMatInCnt = cal_pair_graph_inlier_score(orderedRawMat,affinity.GT,nodeCnt,graphCnt,nodeCnt);
    conDenomMatInCnt = cal_pair_graph_consistency(orderedRawMat,nodeCnt,graphCnt,0);
    
    prevMatching = orderedRawMat(1:nodeCnt*baseGraphCnt, 1:nodeCnt*baseGraphCnt);
    
    param.n = nodeCnt;
    param.graphStep = graphStep;

    for parak = 1:paraCnt
        param.N = baseGraphCnt + (parak-1)*graphStep;
        
        matTmp = orderedRawMat(1:nodeCnt*(param.N+graphStep), 1:nodeCnt*(param.N+graphStep));
        matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N) = prevMatching;
        
        scrDenomMatInCntTmp = scrDenomMatInCnt(1:param.N+graphStep, 1:param.N+graphStep);
        conDenomMatInCntTmp = conDenomMatInCnt(1:param.N+graphStep, 1:param.N+graphStep);
        
        simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
        param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp));

        increMatching = TBIMGM(affinity, simAP, matTmp, param);
        prevMatching = increMatching;
    end
    
    P = increMatching;
end