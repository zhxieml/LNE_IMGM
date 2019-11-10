function P = TBIMGM_offline(rawMat,nodeCnt,graphCnt,param)
    global affinity
    
    if param.adaptive_order
        rawMat = cal_adaptive_graph_order(rawMat, nodeCnt, graphCnt);
    end
    
    baseGraphCnt = 20;
    paraCnt = graphCnt - baseGraphCnt;
    graphStep = 1;
    sigma = 0;
    
    % # outliers = 0
    scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,nodeCnt);
    conDenomMatInCnt = cal_pair_graph_consistency(rawMat,nodeCnt,graphCnt,0);
    
    prevMatching = rawMat(1:nodeCnt*baseGraphCnt, 1:nodeCnt*baseGraphCnt);
    
    param.n = nodeCnt;
    param.graphStep = graphStep;

    for parak = 1:paraCnt
        param.N = baseGraphCnt + (parak-1)*graphStep;
        
        matTmp = rawMat(1:nodeCnt*(param.N+graphStep), 1:nodeCnt*(param.N+graphStep));
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