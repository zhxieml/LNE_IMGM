function X = make_consistent(MST, X)
    flag = false;
    graphCnt = int32(size(MST, 1));
    nodeCnt = idivide(size(X, 1), graphCnt);
    mask = ~eye(graphCnt);
    idx = nodeCnt*(0:graphCnt-1);
    while ~flag
        flag = true;
        for n1 = 1:graphCnt
            connections = mask(n1, :) & MST(n1, :);
            unconnections = mask(n1, :) & ~MST(n1, :);
            n3 = find(unconnections, 1);
            if isempty(n3)
                continue;
            end
            for n2 = find(connections)
                if MST(n2, n3)
                    X13 = X(idx(n1)+1:idx(n1)+nodeCnt, idx(n2)+1:idx(n2)+nodeCnt) * X(idx(n2)+1:idx(n2)+nodeCnt, idx(n3)+1:idx(n3)+nodeCnt);
                    X(idx(n1)+1:idx(n1)+nodeCnt, idx(n3)+1:idx(n3)+nodeCnt) = X13;
                    X(idx(n3)+1:idx(n3)+nodeCnt, idx(n1)+1:idx(n1)+nodeCnt) = X13';
                    %fprintf('connect %d and %d\n', n1, n3);
                    MST(n1, n3) = 1;
                    MST(n3, n1) = 1;
                    break;
                end
            end
        end
    end
end