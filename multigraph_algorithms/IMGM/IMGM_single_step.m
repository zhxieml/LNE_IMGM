function [X, numPairMatch] = IMGM_single_step(globalVar, affScore, rawMat, param)
    %%  single step of Incremental Multi Graph Matching
    % 1. in this algorithm, all graph must have equal # of keypoints
    % affScore is a param.N * param.N matrix
    % param.subMethodParam.name = 'DPMC'
    % param.subMethodParam.useCstDecay = 1
    % param.subMethodParam.cstDecay  = 0.7
    % param.subMethodParam.useWeightedDecay  = 0
    % param.subMethodParam.iterMax = 5
    newGraph = param.N + 1;
    MST = zeros(newGraph, 'logical');
    assert(newGraph == size(affScore,1), 'error:param.N != size(affScore,1)\n');
    % find MST
    
    MST(1:param.N, 1:param.N) = Prim(affScore(1:param.N, 1:param.N));
    nodeCnt = param.n; % # of keypoints per graph
    % initialize hyper graph
    isCenter = sum(MST) > 1;
    %% match members of all center
    centerScore = affScore(newGraph, :).*double(isCenter);
    [~, bestCenter] = max(centerScore);
    %% match members of all edge points of best center
    isConsidered = MST(bestCenter, :) & ~isCenter;
    isConsidered(bestCenter) = true;
    nodeScore = affScore(newGraph, :).*double(isConsidered);
    [~, bestMatch] = max(nodeScore);

    %% connect new graph with best match
    MST(bestMatch, newGraph) = 1;
    MST(newGraph, bestMatch) = 1;
    X = rawMat;
    if newGraph <= 2
        return;
    end
    if param.bVerbose
        fprintf('got best match with %d\n', bestMatch);
        fprintf('before improvment, match score = %.3f\n', affScore(newGraph, bestMatch));
    end
    
    %% apply multigraph algorithms to solve matches in included
    nSubSet = param.maxNumSearch;
    % breadth first search to find nSubSet nearest graphs
    isInSubSet = bfs(MST, newGraph, nSubSet);
    if param.bVerbose
        nSearch = nnz(isInSubSet);
        fprintf("bfs find %d graphs\n", nSearch);
    end
    
    % number of matches
    numPairMatch = sum(isCenter | isConsidered | isInSubSet);

    included = find(isInSubSet);
    excluded = find(~isInSubSet);
    % pairwise match included
    affScore(newGraph, ~isInSubSet) = 0;
    affScore(~isInSubSet, newGraph) = 0;
    % apply multigraph algorithms
    method = param.subMethodParam;
    switch method.name
    case 'CAO'
        % TODO:not finished
        % function P = CAO_test(rawMat,nodeCnt,graphCnt,iterMax,scrDenom,optType,useCstInlier)
        subIndies = getSubIndices(included);
        X(subIndies, subIndies) = CAO(rawMat(subIndies, subIndies),nodeCnt,length(included),method.iterMax,method.scrDenom,method.optType,method.useCstInlier);
    case 'DPMC'
        X = DPMC(globalVar.K, rawMat, affScore, nodeCnt, included, method);
    otherwise
        error('Unexpected sub-multigraph-matching method\n');
    end

    %% make consistent
    stk = zeros(1, newGraph);
    consistent = MST;
    for root = included
        stk(:) = 0;
        visited = isInSubSet;
        top = 1;
        stk(top) = root;
        while(top >= 1)
            t = stk(top);
            visited(t) = true;
            notVisited = MST(t, :) & (~visited);
            if (~nnz(notVisited))
                top = top - 1;
            else
                % t is the father point, f is the adjecent point
                f = find(notVisited, 1);
                R = (root-1)*nodeCnt; %root
                T = (t-1)*nodeCnt;
                F = (f-1)*nodeCnt;
                Xrf = X(R+1:R+nodeCnt, T+1:T+nodeCnt)*X(T+1:T+nodeCnt, F+1:F+nodeCnt);
                Srf = Xrf(:)'*(globalVar.K{root, f}*Xrf(:));
                if Srf > affScore(root, f) 
                    X(R+1:R+nodeCnt, F+1:F+nodeCnt) = Xrf;
                    X(F+1:F+nodeCnt, R+1:R+nodeCnt) = Xrf';
                    affScore(root, f) = Srf;
                    affScore(f, root) = Srf;
                end
                % update consistent
                consistent(root, f) = 1;
                consistent(f, root) = 1;
                % update stack
                stk(top+1) = f;
                top = top + 1;
            end
        end
    end
    [row, col] = find(~consistent(excluded, included));
    for ii = 1:length(row)
        a = excluded(row(ii));
        b = included(col(ii));
        c = included(find(consistent(a, included),1));
        assert(~isempty(c), 'error: cannot find inter point c\n');
        base_a = (a-1)*nodeCnt;
        base_b = (b-1)*nodeCnt;
        base_c = (c-1)*nodeCnt;
        Xac = X(base_a+1:base_a+nodeCnt, base_c+1:base_c+nodeCnt);
        Xcb = X(base_c+1:base_c+nodeCnt, base_b+1:base_b+nodeCnt);
        Xab = Xac*Xcb;
        Sab = Xab(:)'*(globalVar.K{a, b}*Xab(:));
        if Sab > affScore(a, b)
            X(base_a+1:base_a+nodeCnt, base_b+1:base_b+nodeCnt) = Xab;
            X(base_b+1:base_b+nodeCnt, base_a+1:base_a+nodeCnt) = Xab';
            affScore(a, b) = Sab;
            affScore(b, a) = Sab;
        end
    end
end