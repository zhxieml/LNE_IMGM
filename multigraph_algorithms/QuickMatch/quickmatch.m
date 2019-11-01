function X = quickmatch(pointFeat, nodeCnt, param)
    graphCnt = length(pointFeat);
    nodeCnt = int32(nodeCnt);
    data = cell2mat(pointFeat);
    clusters = quickmatch_cluster(data', nodeCnt, param);
    X = zeros(graphCnt*nodeCnt);
    
    inGraph = cellfun(@(x) idivide(x-1,nodeCnt) + 1,clusters, 'UniformOutput', false);
    inPoint = cellfun(@(x) mod(x-1,nodeCnt) + 1,clusters, 'UniformOutput', false);
    for cst = 1:length(clusters)
        if length(clusters{cst}) == 1
            continue;
        end
        for ii = 1:graphCnt-1
            iscope = (ii-1)*nodeCnt+1:ii*nodeCnt;
            for jj = ii+1:graphCnt
                Xij = zeros(nodeCnt);
                jscope = (jj-1)*nodeCnt+1:jj*nodeCnt;
                inPointii = inPoint{cst}(inGraph{cst}==ii);
                inPointjj = inPoint{cst}(inGraph{cst}==jj);
                Xij(inPointii, inPointjj) = 1;
                X(iscope, jscope) = X(iscope, jscope) | Xij;
            end
        end
    end
    X = X + X' + eye(nodeCnt*graphCnt,nodeCnt*graphCnt);
end