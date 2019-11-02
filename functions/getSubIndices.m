function subIndices = getSubIndices(subset, nodeCnt)
    subLength = length(subset);
    subIndices = nodeCnt * ones(nodeCnt, subLength) * diag(subset - ones(1, subLength)) + repmat((1:nodeCnt)', 1, subLength);
    subIndices = subIndices(:);
end