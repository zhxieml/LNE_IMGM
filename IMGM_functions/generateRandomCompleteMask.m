function mask = generateRandomCompleteMask(nodeCnt,graphCnt,complete)
mask = zeros(graphCnt*nodeCnt,graphCnt*nodeCnt);
for x = 1:graphCnt-1
    xscope = (x-1)*nodeCnt+1:x*nodeCnt;
    mask(xscope,xscope) = 0.5;
    for y = x+1:graphCnt
        if rand(1,1)<complete % some are set to 1
            yscope = (y-1)*nodeCnt+1:y*nodeCnt;
            mask(xscope,yscope) = 1;% set to 1, if they are complete. Otherwise, set to 0 as they are not provided by a pairwise solver
        end
    end
end
mask = (mask+mask');