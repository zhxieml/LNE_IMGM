function K = depthFirstSearch(tree, pointSet, iu, ir, np)
% depth first search on MST on u side
% compute the affinity matrix
% X(f) means X_fu
    global assignMat;
    global affinityMatrix;

    r = pointSet(ir);
    K = zeros(np*np);
    nCand = size(tree, 1); % # of points in this cluster
    X = zeros(np, np, nCand);
    
    stk = zeros(1, nCand, 'uint32');
    top = 1;
    visited = zeros(1, nCand, 'logical');
    visited(ir) = true;
    stk(top) = iu;
    X(:, :, iu) = eye(np);
    while top >= 1
        visited(stk(top)) = true;
        adjSet = tree(stk(top), :) & (~visited);
        ff = find(adjSet, 1);
        if isempty(ff)
            top = top - 1;
        else
            % t is the father point, f is the adjecent point
            t = pointSet(stk(top));
            f = pointSet(ff);
            Xft = reshape(assignMat(:, t, f), np, np);
            X(:, :, ff) = Xft * X(:, :, stk(top)); % Xfu = Xft * Xtu
            
            % compute affinity matrix
            Ffu = kron(eye(np), X(:, :, ff));
            Krf = affinityMatrix(:, :, r, f);
            K = K + Ffu'*(Krf*Ffu); 
            % update stack
            top = top + 1;
            stk(top) = ff;
        end
    end
end
