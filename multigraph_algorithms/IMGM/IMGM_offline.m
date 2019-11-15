function [P, numPairMatch] = IMGM_offline(rawMat,nodeCnt,graphCnt,param)
    global affinity target
        
    numPairMatch = 0;
    
    baseGraphCnt = 0;
    graphStep = 1;
    paraCnt = (graphCnt - baseGraphCnt) / graphStep;
    
    prevMatching = rawMat(1:nodeCnt*baseGraphCnt, 1:nodeCnt*baseGraphCnt);
    
    param.n = nodeCnt;
    param.graphStep = graphStep;

    for parak = 1:paraCnt
        param.N = baseGraphCnt + (parak-1)*graphStep;
        
        matTmp = rawMat(1:nodeCnt*(param.N+graphStep), 1:nodeCnt*(param.N+graphStep));
        matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N) = prevMatching;
        
%         scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp,affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
        
%         simAP = scrDenomMatInCntTmp;

                                   %         IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, param)
        [increMatching, increNumPairMatch] = IMGM_batch(affinity, target, matTmp, nodeCnt, param.N, graphStep, 0, param);
        prevMatching = increMatching;
        numPairMatch = numPairMatch + increNumPairMatch;

%         increMatching = IMGM_old(simAP, matTmp, param);
%         prevMatching = increMatching;
%         numPairMatch = numPairMatch + param.N + graphStep;
    end
    
    P = increMatching;
end