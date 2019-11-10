function [affinity]= generateAffinity(testk)
global target
    if target.config.database == "synthetic"
        affinity = generateRandomAffinity(testk);
        graphCnt = target.config.graphCnt;
        nodeCnt = target.config.nodeCnt;
        affinity.GT = repmat(eye(nodeCnt), graphCnt, graphCnt);
    else
        affinity = generateRealAffinity(testk);
    end
end
