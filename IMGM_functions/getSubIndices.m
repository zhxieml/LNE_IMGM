function subIndices = getSubIndices(subset, nodeCnt)
    % subLength = length(subset);
    % subIndices = nodeCnt * ones(nodeCnt, subLength) * diag(subset - ones(1, subLength)) + repmat((1:nodeCnt)', 1, subLength);
    % subIndices = subIndices(:);
    subIndices = [];
    for ii = subset
        subIndices = [subIndices, (ii-1)*nodeCnt+1:ii*nodeCnt];
    end
end