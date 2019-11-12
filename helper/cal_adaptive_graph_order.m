function order = cal_adaptive_graph_order(X,nodeCnt,graphCnt)
    unaryListConsistency = cal_single_graph_consistency(X,nodeCnt,graphCnt,0,0);
    [~, order] = sort(unaryListConsistency, 'descend');
end