function Ku = dfs(K, X, hyperGraph, iu, ir, graphOf, propRate, minPropRate, nodeCnt)
% depth first search on graph with depth < maxDepth
    len = length(graphOf);
    Xu = zeros(nodeCnt, nodeCnt, len);
    stk = zeros(1, len);
    stkPropRate = zeros(1, len);
    Ku = zeros(nodeCnt*nodeCnt);
    visited = zeros(1, len);
    top = 1;
    visited(ir) = true;
    stk(top) = iu;
    stkPropRate(top) = propRate;
    Xu(:, :, iu) = eye(nodeCnt);
    while top >= 1
        iT = stk(top);
        visited(iT) = true;
        notVisited = hyperGraph(graphOf(iT), graphOf) & (~visited);
        if (~nnz(notVisited)) || (stkPropRate(top) < minPropRate)
            top = top - 1;
        else
            % iT is the father point, iF is the adjecent point
            iF = find(notVisited, 1);
            f = (graphOf(iF)-1)*nodeCnt;
            t = (graphOf(iT)-1)*nodeCnt;
            Xu(:, :, iF) = X(f+1:f+nodeCnt, t+1:t+nodeCnt)*Xu(:, :, iT); % Xu = Xft * Xtu
            % compute affinity matrix
            Ffu = kron(Xu(:, :, iF), eye(nodeCnt));
            Ku = Ku + stkPropRate(top)*(Ffu'*(K{graphOf(ir), graphOf(iF)}*Ffu)); 
            % update stack
            stk(top+1) = iF;
            stkPropRate(top+1) = stkPropRate(top)*propRate;
            top = top + 1;
        end
    end
end
