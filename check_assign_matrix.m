function is_valid = check_assign_matrix(assignMat, nodeCnt)
    matSize = size(assignMat);
    graphCnt = matSize(1) / nodeCnt;
    
    if size(unique(assignMat >= 0), 1) ~= 1 || size(unique(assignMat <= 1), 1) ~= 1
        pause;
        error("Invalid input: the assignment matrix has error entries");
    end
    
    if matSize(1) ~= matSize(2) || mod(matSize(1), nodeCnt) ~= 0
        error("Invalid input: the assignment matrix and the number of nodes don't match");
    end
    
    for i = 1:graphCnt
        indies = (i-1)*nodeCnt+1:i*nodeCnt;
        
        oneMat = ones(nodeCnt, 1);
        cropMat = assignMat(indies, indies);
        
        % the doubly stochastic constraints
        constraint_11 = size(unique(cropMat*oneMat <= 1), 1) == 1;
        constraint_12 = size(unique(cropMat*oneMat >= 0), 1) == 1;
        
        constraint_21 = size(unique(cropMat'*oneMat <= 1), 1) == 1;
        constraint_22 = size(unique(cropMat'*oneMat >= 0), 1) == 1;
        
        if constraint_11 && constraint_12 && constraint_21 && constraint_22
            continue;
        end            
        
        fprintf("(%d, %d) pair does not satisfy the doubly stochastic constraints\n");
        pause;
    end
end