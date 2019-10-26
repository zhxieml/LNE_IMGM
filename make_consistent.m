function res = make_consistent(MST, X, varargin)
    flag = false;
    graphCnt = int32(size(MST, 1));
    mask = ~eye(graphCnt);
    if nargin == 3
        nodeCnt = idivide(size(X, 1), graphCnt);
        idx = nodeCnt*(varargin{1}-1);
        res = X;
    else
        res = MST;
    end
    
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
                    if nargin == 3
                        X13 = res(idx(n1)+1:idx(n1)+nodeCnt, idx(n2)+1:idx(n2)+nodeCnt) * res(idx(n2)+1:idx(n2)+nodeCnt, idx(n3)+1:idx(n3)+nodeCnt);
                        res(idx(n1)+1:idx(n1)+nodeCnt, idx(n3)+1:idx(n3)+nodeCnt) = X13;
                        res(idx(n3)+1:idx(n3)+nodeCnt, idx(n1)+1:idx(n1)+nodeCnt) = X13';
                    else
                        res(n1, n3) = res(n1, n2)*res(n2, n3);
                        res(n3, n1) = res(n1, n3);
                    end
                    %fprintf('connect %d and %d\n', n1, n3);
                    MST(n1, n3) = 1;
                    MST(n3, n1) = 1;
                    break;
                end
            end
        end
    end
end