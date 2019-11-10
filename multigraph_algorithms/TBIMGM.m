function [X, numPairMatch] = TBIMGM(globalVar, affScore, rawMat, param)
    %%  single step of Incremental Multi Graph Matching
    global target
    graphCnt = param.N + param.graphStep;
    MST = zeros(graphCnt, 'logical');
    % find MST
    
    MST(1:param.N, 1:param.N) = Prim(affScore(1:param.N, 1:param.N));
    nodeCnt = param.n; % # of keypoints per graph
    % initialize hyper graph
    isCenter = sum(MST) > 1;
    %% match members of all center
    centerScore = affScore(graphCnt, :).*double(isCenter);
    [~, bestCenter] = max(centerScore);
    %% match members of all edge points of best center
    isConsidered = MST(bestCenter, :) & ~isCenter;
    isConsidered(bestCenter) = true;
    nodeScore = affScore(graphCnt, :).*double(isConsidered);
    [~, bestMatch] = max(nodeScore);

    %% connect new graph with best match
    MST(bestMatch, graphCnt) = 1;
    MST(graphCnt, bestMatch) = 1;
    X = rawMat;
    if graphCnt <= 2
        return;
    end
    if param.bVerbose
        fprintf('got best match with %d\n', bestMatch);
        fprintf('before improvment, match score = %.3f\n', affScore(graphCnt, bestMatch));
    end
    
    %% apply multigraph algorithms to solve matches in included
    nSubSet = param.maxNumSearch;
    % breadth first search to find nSubSet nearest graphs
    isInSubSet = bfs(MST, graphCnt, nSubSet);
    if param.bVerbose
        nSearch = nnz(isInSubSet);
        fprintf("bfs find %d graphs\n", nSearch);
    end
    
    % number of matches
    numPairMatch = sum(isCenter | isConsidered | isInSubSet);

    included = find(isInSubSet);
    excluded = find(~isInSubSet);
    % pairwise match included
    affScore(graphCnt, ~isInSubSet) = 0;
    affScore(~isInSubSet, graphCnt) = 0;
    % apply multigraph algorithms
    method = param.subMethodParam;
    subIndies = getSubIndices(included, nodeCnt);
    switch method.name
    case 'CAO'
        %%%%%%%%%%%%%%%% crop the affinity %%%%%%%%%%%%%%%%%%%%%
        affinityCrop = crop_affinity(included', globalVar);
        %%%%%%%%%%%%%%%% crop the target %%%%%%%%%%%%%%%%%%%%%%%
        targetCrop = crop_target(included', target);
        
        X(subIndies, subIndies) = CAO_local(rawMat(subIndies, subIndies),nodeCnt,length(included),method.iterMax,method.scrDenom,affinityCrop,targetCrop,method.optType,method.useCstInlier);
%         X(subIndies, subIndies) = CAO(rawMat(subIndies, subIndies),nodeCnt,length(included),method.iterMax,method.scrDenom,method.optType,method.useCstInlier);
    case 'quickmatch'
%         pointFeat = globalVar.pointFeat(included);
%         X(subIndies, subIndies) = quickmatch(pointFeat, nodeCnt, method);
        nFeature = ones(1, length(included)) * nodeCnt;
        M_out = QuickMatch(rawMat(subIndies, subIndies), nFeature);
        X(subIndies, subIndies) = double(full(M_out));

    case 'matchALS'
        nFeature = ones(1, length(included)) * nodeCnt;
        M_out = mmatch_CVX_ALS(rawMat(subIndies, subIndies), nFeature, 'verbose', false, 'univsize', length(subIndies));
        X(subIndies, subIndies) = double(full(M_out));
    otherwise
        error('Unexpected sub-multigraph-matching method\n');
    end

    %% make consistent

    
    % r = graphCnt;
    % rview = (r-1)*nodeCnt+1:r*nodeCnt;
    % for x = included(1:end-1)
    %     % if x == r
    %     %     continue;
    %     % end
    %     xview = (x-1)*nodeCnt+1:x*nodeCnt;
    %     for y = excluded
    %         yview = (y-1)*nodeCnt+1:y*nodeCnt;
    %         Xry = X(rview, xview)*X(xview, yview);
    %         Sry = mat2vec(Xry)'*(globalVar.K{r, y}*mat2vec(Xry));
    %         if Sry > affScore(r, y)
    %             X(rview, yview) = Xry;
    %             X(yview, rview) = Xry';
    %             affScore(r, y) = Sry;
    %             affScore(y, r) = Sry;
    %         end
    %     end
    % end

    stk = zeros(1, graphCnt);
    consistent = MST;
    for r = included
        stk(:) = 0;
        visited = isInSubSet;
        top = 1;
        stk(top) = r;
        while(top >= 1)
            t = stk(top);
            visited(t) = true;
            notVisited = MST(t, :) & (~visited);
            if (~nnz(notVisited))
                top = top - 1;
            else
                % t is the father point, f is the adjecent point
                f = find(notVisited, 1);
                view_r = (r-1)*nodeCnt+1:r*nodeCnt;
                view_t = (t-1)*nodeCnt+1:t*nodeCnt;
                view_f = (f-1)*nodeCnt+1:f*nodeCnt;
                Xrf = X(view_r, view_t)*X(view_t, view_f);
                Srf = mat2vec(Xrf)'*(globalVar.K{r, f}*mat2vec(Xrf));
                if Srf > affScore(r, f) 
                    X(view_r, view_f) = Xrf;
                    X(view_f, view_r) = Xrf';
                    affScore(r, f) = Srf;
                    affScore(f, r) = Srf;
                end
                % update consistent
                consistent(r, f) = 1;
                consistent(f, r) = 1;
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
        view_a = (a-1)*nodeCnt+1:a*nodeCnt;
        view_b = (b-1)*nodeCnt+1:b*nodeCnt;
        view_c = (c-1)*nodeCnt+1:c*nodeCnt;
        Xab = X(view_a, view_c)*X(view_c, view_b);
        Sab = mat2vec(Xab)'*(globalVar.K{a, b}*mat2vec(Xab));
        if Sab > affScore(a, b)
            X(view_a, view_b) = Xab;
            X(view_b, view_a) = Xab';
            affScore(a, b) = Sab;
            affScore(b, a) = Sab;
        end
    end

end


