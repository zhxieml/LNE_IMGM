function [X] = IMGM_single_step(globalVar, affScore, rawMat, param)
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
    nodeScore = affScore(newGraph, :).*double(~isCenter);
    nodeScore(bestCenter) = affScore(newGraph, bestCenter);
    [~, bestMatch] = max(nodeScore);

    %% connect new graph with best match
    MST(bestMatch, newGraph) = 1;
    MST(newGraph, bestMatch) = 1;
    if newGraph <= 2
        X = rawMat;
        return;
    end
    if param.bVerbose
        fprintf('got best match with %d\n', bestMatch);
        fprintf('before improvment, match score = %.3f\n', affScore(newGraph, bestMatch));
    end
    
    %% apply multigraph algorithms to solve matches in subSet
    nSubSet = param.maxNumSearch;
    % breadth first search to find nSubSet nearest graphs
    isInSubSet = bfs(MST, newGraph, nSubSet);
    if param.bVerbose
        nSearch = nnz(isInSubSet);
        fprintf("bfs find %d graphs\n", nSearch);
    end
    
    subSet = find(isInSubSet);
    % apply multigraph algorithms
    method = param.subMethodParam;
    switch method.name
    case 'CAO'
        % TODO:not finished
        X = CAO(rawMat,nodeCnt,nSubSet,param.iterMax,param.scrDenom,param.optType,param.useCstInlier);
    case 'DPMC'
        X = DPMC(globalVar.K, rawMat, affScore, nodeCnt, subSet, method);
    otherwise
        error('Unexpected sub-multigraph-matching method\n');
    end

    %% make consistent
    stk = zeros(1, newGraph);
    for root = subSet
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
                Srf = Xrf(:)'*(globalVar.K{r, f}*Xrf(:));
                if affScore(r, f) < Srf
                    X(R+1:R+nodeCnt, F+1:F+nodeCnt) = Xrf;
                    X(F+1:F+nodeCnt, R+1:R+nodeCnt) = Xrf';
                end
                % update stack
                stk(top+1) = f;
                top = top + 1;
            end
        end
    end
end