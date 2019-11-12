function P = IMGM_batch(rawMat, prevMatching, nodeCnt, baseGraphCnt, graphStep, param)
    global affinity;
    param.n = nodeCnt;
    param.graphStep = graphStep;
    assert(size(prevMatching, 1) == nodeCnt*baseGraphCnt, "error in IMGM_batch\n");
    for parak = 1:graphStep
        param.N = baseGraphCnt+parak-1;
        graphCnt = param.N+1;
        matTmp = rawMat(1:nodeCnt*graphCnt,1:nodeCnt*graphCnt);
        matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching;
        simAP = cal_pair_graph_inlier_score(matTmp,affinity.GT(1:nodeCnt*graphCnt,1:nodeCnt*graphCnt),nodeCnt,graphCnt,nodeCnt);
        increMatching = IMGM_old(simAP, matTmp, param);
        prevMatching = increMatching;
    end
    P = increMatching;
end