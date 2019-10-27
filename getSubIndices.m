function subIndices = getSubIndices(subset)
    subLength = length(subset);
    subIndices = 10 * ones(10, subLength) * diag(subset - ones(1, subLength)) + repmat((1:10)', 1, subLength);
    subIndices = subIndices(:);
end