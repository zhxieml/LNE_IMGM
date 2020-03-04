function X = DPMC(K, X, affScore, nodeCnt, subSet, param)
%TPMC AAAI20 paper multi-graph algorithms
% 
    S = affScore(subSet, subSet);
    MST = Prim(S);
    M = getTrans(nodeCnt);
    [xDim, yDim] = meshgrid(1:nodeCnt, 1:nodeCnt);
    matchlist = [xDim(:), yDim(:)]';
    [group1, group2] = make_group12(matchlist);
    if param.useCstDecay
        gama = make_consistent(MST, param.cstDecay);
    elseif param.useWeightedDecay
        minWeight = min(S(MST));
        maxWeight = max(S(MST));
        gama = (S-minWeight) / (maxWeight-minWeight);
    end
        
    [row, col] = find(triu(MST));
    for itr = 1:param.iterMax
        for ii = 1:length(row)
            %% solve composed QAP sub problem
            iu = row(ii);
            ir = col(ii);
            u = subSet(iu);
            r = subSet(ir);
            Ku = DFS(K, X, MST, gama, subSet, iu, ir, nodeCnt);
            Kr = DFS(K, X, MST, gama, subSet, ir, iu, nodeCnt);
            Kur = K{u, r} + Kr + M*(Ku*M);
            Xur_raw = RRWM(Kur, group1, group2);
            Xur_flat = greedyMapping(Xur_raw, group1, group2);
            Xur = reshape(Xur_flat, nodeCnt, nodeCnt);
            U = (u-1)*nodeCnt;
            R = (r-1)*nodeCnt;
            X(U+1:U+nodeCnt, R+1:R+nodeCnt) = Xur;
            % reverse tranform
            X(R+1:R+nodeCnt, U+1:U+nodeCnt)= Xur';
        end
    end
    % make full consistent
    X = make_consistent(MST, X, nodeCnt, subSet);
end