function [newX, order] = cal_adaptive_graph_order(X,nodeCnt,graphCnt)
    unaryListConsistency = cal_single_graph_consistency(X,nodeCnt,graphCnt,0,0);
    
    [~, order] = sort(unaryListConsistency, 'descend');
    
    order = order';
    
    subIndies = getSubIndices(order, nodeCnt);
    newX = X(subIndies, subIndies);
end