function res = make_consistent(MST, X, varargin)
    flag = false;
    graphCnt = int32(size(MST, 1));
    mask = ~eye(graphCnt);
    if nargin == 4
        nodeCnt = varargin{1,1};
        base = varargin{1,2}-1;
        base = base*nodeCnt;
        res = X;
    else
        res = double(MST)*X;
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
            flag = false;
            for n2 = find(connections)
                if MST(n2, n3)
                    if nargin == 4
                        X13 = res(base(n1)+1:base(n1)+nodeCnt, base(n2)+1:base(n2)+nodeCnt) * res(base(n2)+1:base(n2)+nodeCnt, base(n3)+1:base(n3)+nodeCnt);
                        res(base(n1)+1:base(n1)+nodeCnt, base(n3)+1:base(n3)+nodeCnt) = X13;
                        res(base(n3)+1:base(n3)+nodeCnt, base(n1)+1:base(n1)+nodeCnt) = X13';
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