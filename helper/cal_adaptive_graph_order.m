function order = cal_adaptive_graph_order(affinity, X, nodeCnt, graphCnt, opt)
    assert(graphCnt >= 1);
    order = 1:graphCnt;
    if graphCnt <= 2
        return
    end
    unaryListConsistency = cal_single_graph_consistency(X,nodeCnt,graphCnt,0,0);
    unaryListAffinity = cal_single_graph_affinity_local(affinity,X,nodeCnt,graphCnt);
    sigma = 1;
    sortList = sigma * unaryListAffinity + (1-sigma) * unaryListConsistency;
    switch opt
    case 'd'
        [~, order] = sort(sortList, 'descend');
    case 'a'
        [~, order] = sort(sortList, 'ascend');
    otherwise
        error("undefined adaptive graph order\n");
    end
    
end