function order = cal_adaptive_graph_order(X,nodeCnt,graphCnt, opt)
    assert(graphCnt >= 1);
    order = 1:graphCnt;
    if graphCnt <= 2
        return
    end
    unaryListConsistency = cal_single_graph_consistency(X,nodeCnt,graphCnt,0,0);
    switch opt
    case 'd'
        [~, order] = sort(unaryListConsistency, 'descend');
    case 'a'
        [~, order] = sort(unaryListConsistency, 'ascend');
    otherwise
        error("undefined adaptive graph order\n");
    end
    
end