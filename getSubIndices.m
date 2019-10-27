function subIndices = getSubIndices(subset)
    subLength = length(subset);
    subIndices = 10 * diag(subset - ones(1, subLength)) * ones(subLength, 10) + repmat(1:10, subLength, 1);
end