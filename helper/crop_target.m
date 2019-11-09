% a simple function to crop the a target of a specific cluster
function targetCrop = crop_target(tmpClusterPosition, target)
    targetCrop = target;
    tmp = double([1:length(tmpClusterPosition)]);
    targetCrop.config.selectGraphMask = {tmp};
    targetCrop.config.graphCnt = length(tmpClusterPosition);
    tmpMask = target.pairwiseMask{1};
    nDivide = ones([1 target.config.graphCnt])*targetCrop.config.nodeCnt;
    tmpMaskCell = mat2cell(tmpMask,nDivide,nDivide);
    tmpMaskCellCrop = cell(length(tmpClusterPosition), length(tmpClusterPosition));
    for i = 1:length(tmpClusterPosition)
        for j = 1:length(tmpClusterPosition)
            tmpMaskCellCrop(i,j) = tmpMaskCell(tmpClusterPosition(i),tmpClusterPosition(j));
        end
    end
    maskCrop = cell2mat(tmpMaskCellCrop);
    targetCrop.pairwiseMask = {maskCrop};
end