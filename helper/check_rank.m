function rankMat = check_rank(M, nodeCnt, graphCnt)
    rankMat = zeros(graphCnt, graphCnt);
    
    for xview = 1:graphCnt
        for yview = 1:graphCnt
            xscope = (xview-1)*nodeCnt+1:xview*nodeCnt;
            yscope = (yview-1)*nodeCnt+1:yview*nodeCnt;
            
            rankMat(xview, yview) = rank(M(xscope, yscope));
        end
    end
end