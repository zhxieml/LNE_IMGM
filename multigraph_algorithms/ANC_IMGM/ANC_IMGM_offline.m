function [P, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,batchSize,useAptOrder,opt,param)
    global affinity target
    
    baseGraphCnt = 0;
    paraCnt = graphCnt - baseGraphCnt;
    
    aptOrder = 1:graphCnt;

    numPairMatch = 0;
    
    % # outliers = 0
%     scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,nodeCnt);
    
    prevMatching = rawMat(1:nodeCnt*baseGraphCnt, 1:nodeCnt*baseGraphCnt);        
    
    param.n = nodeCnt;
%     param.graphStep = graphStep;

    for parak = 1:paraCnt
        param.N = baseGraphCnt + parak - 1;
        
%         matTmp = rawMat(1:nodeCnt*(param.N+graphStep), 1:nodeCnt*(param.N+graphStep));
%         matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N) = prevMatching;
%         
%         scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp,affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
%         
% %         simAP = scrDenomMatInCntTmp;
%         param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
%         
%         [increMatching, increNumPairMatching] = ANC_IMGM_batch(affinity, target, matTmp, nodeCnt, param.N, batchSize, useAptOrder, opt, param);
%         
%         prevMatching = increMatching;
%         numPairMatch = numPairMatch + increNumPairMatching;

        if useAptOrder
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, opt);
                subOrder = subOrder + param.N;
                aptOrder(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder = crop_rawMat(aptOrder, rawMat, nodeCnt);
                affinityReOrder = crop_affinity(aptOrder, affinity);
                targetReOrder = crop_target(aptOrder, target);
            end

            matTmp = rawMatReOrder(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching;
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder, matTmp, nodeCnt, param.N+1, nodeCnt);
            [increMatching, increNumPairMatching] = ANC_IMGM(affinityReOrder, affScore, matTmp, targetReOrder, param);
            prevMatching = increMatching;
            numPairMatch = numPairMatch + increNumPairMatching;
        else
            matTmp = rawMat(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching;
            affScore = cal_pair_graph_inlier_score_local(affinity, matTmp, nodeCnt, param.N+1, nodeCnt);
            [increMatching, increNumPairMatching] = ANC_IMGM(affinity, affScore, matTmp, target, param);
            prevMatching = increMatching;
            numPairMatch = numPairMatch + increNumPairMatching;
        end
    end
    
    P = increMatching;
end