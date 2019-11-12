function [newX] = crop_rawMat(tmpClusterPosition, rawMat, nodeCnt)
    subIndies = getSubIndices(tmpClusterPosition, nodeCnt);
    newX = rawMat(subIndies, subIndies);
end