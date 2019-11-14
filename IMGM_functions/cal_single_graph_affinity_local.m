function afPerGraph = cal_single_graph_affinity_local(affinity,X,nodeCnt,graphCnt)
    afPerGraph = zeros(graphCnt, 1);
    for r = 1:graphCnt
        rscope = (r-1)*nodeCnt+1:r*nodeCnt;
        iview = [1:r-1, r+1:graphCnt];
        for i = iview
            iscope = (i-1)*nodeCnt+1:i*nodeCnt;
            p = mat2vec(X(rscope, iscope));
            afPerGraph(r) = afPerGraph(r) + p'*affinity.K{r, i}*p;
        end
    end
    afPerGraph = afPerGraph ./max(afPerGraph);
end