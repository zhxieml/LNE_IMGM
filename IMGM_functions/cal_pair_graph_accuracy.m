function accPair = cal_pair_graph_accuracy(X,GT,nOutlier,nodeCnt,graphCnt)
% 输入X只要上半部分完整就行viewx<viewy
accPair = zeros(graphCnt,graphCnt);
 for viewx=1:graphCnt
     xscope = (viewx-1)*nodeCnt+1:viewx*nodeCnt;
    for viewy = viewx+1:graphCnt
        
        yscope = (viewy-1)*nodeCnt+1:viewy*nodeCnt;
        accPair(viewx,viewy) = cal_acc(X(xscope,yscope),nOutlier,GT(xscope,yscope));
    end
 end
 accPair = accPair + accPair' + eye(graphCnt);%这会使得整体精度偏高，而且图的数目少的时候更明显