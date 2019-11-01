function X = quickmatch(pointFeat, nodeCnt, param)
    graphCnt = int32(length(pointFeat));
    data = cell2mat(pointFeat);
    data = data(:);
    clusters = quickmatch_cluster(data', nodeCnt, param);
    X = zeros(graphCnt*nodeCnt);
    
    inGraph = cellfun(@(x) idivide(x-1,graphCnt) + 1,clusters);
    inPoint = cellfun(@(x) mod(x-1,graphCnt) + 1,clusters);
    for cst = clusters
        for ii = 1:graphCnt-1
            iscope = (ii-1)*nodeCnt+1:ii*nodeCnt;
            for jj = ii+1:graphCnt
                Xij = zeros(nodeCnt);
                jscope = (jj-1)*nodeCnt+1:jj*nodeCnt;
                inPointii = inPoint{cst}(inGraph{cst}==ii);
                inPointjj = inPoint{cst}(inGraph{cst}==jj);
                Xij(inPointii, inPointjj) = 1;
                X(iscope, jscope) = Xij;
            end
        end
    end
    X = X + X' + eye(nodeCnt*graphCnt,nodeCnt*graphCnt);
end