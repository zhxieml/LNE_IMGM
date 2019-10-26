function [X, hyperGraph] = IMGM_single_step(globalVar, affScore, rawMat, hyperGraph, param)
    %%  single step of Incremental Multi Graph Matching
    % 1. in this algorithm, all graph must have equal # of keypoints
    % 3. �?要修改的地方
    %    1.hyperGraph 是图之间的连接，用于维护�?颗MST
    %    2.globalVar 是全�?变量，即原代码中的global affinity,
    %    3. param 中需要加入以下参数， propRate(0.5), minPropRate(0.21), bVerbose(0), maxNumSearch(20) 
    %       param.N 是已经匹配好的图的个数，本次匹配第N+1�?(见第20�?)

    graphCnt = size(affScore, 1); % # of graphs
    nodeCnt = param.n; % # of keypoints per graph
    % initialize hyper graph
    isCenter = sum(hyperGraph) > 1;
    % calculate transF for further use
    transF = getTrans(nodeCnt);
    % calculate group1 group2 for pairwise matching
    [X, Y] = meshgrid(1:nodeCnt, 1:nodeCnt);
    matchlist = [Y(:), X(:)]';
    [group1, group2] = make_group12(matchlist);

    iNewGraph = param.N + 1;
    %% match members of all center
    centerScore = affScore(iNewGraph, :).*double(isCenter);
    [~, bestCenter] = max(centerScore);
    %% match members of all edge points of best center
    nodeScore = affScore(iNewGraph, :).*double(~isCenter);
    nodeScore(bestCenter) = affScore(iNewGraph, bestCenter);
    [~, iBestMatch] = max(nodeScore);

    % connect new graph with best match
    hyperGraph(iBestMatch, iNewGraph) = 1;
    hyperGraph(iNewGraph, iBestMatch) = 1;
    if iNewGraph <= 2
        return;
    end
    if param.bVerbose
        fprintf('got best match with %d\n', iBestMatch);
        fprintf('before improvment, match score = %.3f\n', affScore(iNewGraph, iBestMatch));
    end
    visited = bfs(hyperGraph, iNewGraph, param.maxNumSearch);
    if param.bVerbose
        nSearch = nnz(visited);
        fprintf("bfs find %d graphs\n", nSearch);
    end

    graphSet = find(visited);
    
    %% apply Prim algorithm to find MST
    hyperGraph(visited, visited) = 0;
    searched = zeros(1, graphCnt, 'logical');
    [~, I] = sort(affScore(:), 'descend');
    
    I = int32(I);
    graphCnt = int32(graphCnt);
    
    a = mod(I - 1, graphCnt) + 1;
    b = idivide(I - 1, graphCnt) + 1;
    searched(a(1)) = true;
    for ii = 1:length(I)
        if searched(a(ii)) && visited(b(ii)) && ~searched(b(ii))
            searched(b(ii)) = 1;
            hyperGraph(a(ii), b(ii)) = 1;
            hyperGraph(b(ii), a(ii)) = 1;
        end
    end

    [row, col] = find(triu(hyperGraph(visited, visited)));
    X = rawMat;
    for ii = 1:length(row)
        %% solve composed QAP sub problem
        iu = row(ii);
        ir = col(ii);
        u = graphSet(iu);
        r = graphSet(ir);
        Ku = dfs(globalVar.K, X, hyperGraph, iu, ir, graphSet, param.propRate, param.minPropRate, nodeCnt);
        Kr = dfs(globalVar.K, X, hyperGraph, iu, ir, graphSet, param.propRate, param.minPropRate, nodeCnt);
        K = globalVar.K{u, r} + Kr + transF*(Ku*transF);
        Xur_raw = RRWM(K, group1, group2);
        Xur_dis = greedyMapping(Xur_raw, group1, group2);
        Xur = reshape(Xur_dis, nodeCnt, nodeCnt);
        U = (u-1)*nodeCnt;
        R = (r-1)*nodeCnt;
        X(U+1:U+nodeCnt, R+1:R+nodeCnt) = Xur;
        % reverse tranform
        X(R+1:R+nodeCnt, U+1:U+nodeCnt)= Xur';
        affScore(u, r) = Xur_dis'*(globalVar.K{u, r}*Xur_dis);
        affScore(r, u) = affScore(u, r);
    end
    % make consistent
    stk = zeros(1, graphCnt);
    visited(:) = 0;
    top = 1;
    stk(top) = iNewGraph;
    while(top >= 1)
        t = stk(top);
        visited(t) = true;
        notVisited = hyperGraph(t, :) & (~visited);
        if (~nnz(notVisited))
            top = top - 1;
        else
            % t is the father point, f is the adjecent point
            f = find(notVisited, 1);
            R = (iNewGraph-1)*nodeCnt; %root
            T = (t-1)*nodeCnt;
            F = (f-1)*nodeCnt;
            Xgf = X(R+1:R+nodeCnt, T+1:T+nodeCnt)*X(T+1:T+nodeCnt, F+1:F+nodeCnt);
            X(R+1:R+nodeCnt, F+1:F+nodeCnt) = Xgf;
            X(F+1:F+nodeCnt, R+1:R+nodeCnt) = Xgf';
            % update stack
            stk(top+1) = f;
            top = top + 1;
        end
    end
    
end