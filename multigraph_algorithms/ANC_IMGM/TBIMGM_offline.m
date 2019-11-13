function [P, numPairMatch] = TBIMGM_offline(rawMat,nodeCnt,graphCnt,batchSize,useAptOrder,param)
    global affinity target
    
%     if param.adaptive_order
%         rawMat = cal_adaptive_graph_order(rawMat, nodeCnt, graphCnt);
%     end
    
    baseGraphCnt = 0;
    graphStep = batchSize;
    paraCnt = (graphCnt - baseGraphCnt) / graphStep;
    
%     if paraCnt == 0
%         matTmp = rawMat(1:nodeCnt*(param.N+graphStep), 1:nodeCnt*(param.N+graphStep));
%         
%         scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp,affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
%         
%         param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
%         
%         [increMatching, increNumPairMatching] = ANC_IMGM_batch(affinity, target, matTmp, nodeCnt, param.N, batchSize, useAptOrder, param);
%         
%         prevMatching = increMatching;
%         numPairMatch = numPairMatch + increNumPairMatching;
%         
%     end

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
        
%         simAP = scrDenomMatInCntTmp;
        param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
        
        [increMatching, increNumPairMatching] = ANC_IMGM_batch(affinity, target, matTmp, nodeCnt, param.N, batchSize, useAptOrder, param);
        
        prevMatching = increMatching;
        numPairMatch = numPairMatch + increNumPairMatching;


%         [increMatching, increNumPairMatching] = TBIMGM(affinity, simAP, matTmp, param);
%         prevMatching = increMatching;
%         numPairMatch = numPairMatch + increNumPairMatching;
    end
    
    P = increMatching;
end