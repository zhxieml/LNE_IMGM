function Ku = dfs(K, X, MST, gama, subSet, iu, ir, nodeCnt)
% depth first search on graph with depth < maxDepth
    len = length(subSet);
    Xu = zeros(nodeCnt, nodeCnt, len);
    stk = zeros(1, len);
    Ku = zeros(nodeCnt*nodeCnt);
    visited = zeros(1, len);
    top = 1;
    visited(ir) = true;
    stk(top) = iu;
    Xu(:, :, iu) = eye(nodeCnt);
    while top >= 1
        iT = stk(top);
        visited(iT) = true;
        notVisited = MST(iT, :) & (~visited);
        if (~nnz(notVisited))
            top = top - 1;
        else
            % iT is the father point, iF is the adjecent point
            iF = find(notVisited, 1);
            f = (subSet(iF)-1)*nodeCnt;
            t = (subSet(iT)-1)*nodeCnt;
            Xu(:, :, iF) = X(f+1:f+nodeCnt, t+1:t+nodeCnt)*Xu(:, :, iT); % Xu = Xft * Xtu
            % compute affinity matrix
            Ffu = kron(Xu(:, :, iF), eye(nodeCnt));
            Ku = Ku + gama(iu, iF)*(Ffu'*(K{subSet(ir), subSet(iF)})*Ffu); 
            % update stack
            stk(top+1) = iF;
            top = top + 1;
        end
    end
end
