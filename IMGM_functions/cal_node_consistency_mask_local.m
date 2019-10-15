% Input
% X: the pairwise matching configuration
% nodeCnt: total number of nodes for each graph, assuming each graph is of equal size for simplicity
% graphCnt: number of graphs for matching
% inCnt: the number of inliers, which is specified by a certain means
% 
% Output
% nodeCstPerGraph: each row is a node, each column is a graph, for each
% column, only inCnt elements are 1, others are 0
function nodeCstPerGraph = cal_node_consistency_mask_local(X,nodeCnt,graphCnt,inCnt)
% node-wise consistency of each node in each graph, 
nodeCstPerGraph = zeros(nodeCnt,graphCnt);
% compute all nodes' node-wise consistency in a batch for each graph 
for ref = 1:graphCnt
    rscope = (ref-1)*nodeCnt+1:ref*nodeCnt;
    for i = 1:graphCnt-1
        iscope = (i-1)*nodeCnt+1:i*nodeCnt;
        for j = i+1:graphCnt
            jscope = (j-1)*nodeCnt+1:j*nodeCnt;
            Xrij=X(rscope,iscope)*X(iscope,jscope);
            nodeCstPerGraph(:,ref) = nodeCstPerGraph(:,ref) + (1-sum(abs(Xrij-X(rscope,jscope)),2)/2);
        end
    end
end
% normalize the summation value
nodeCstPerGraph = nodeCstPerGraph/(graphCnt*(graphCnt-1)/2);
% sort the 
[~,IX] = sort(nodeCstPerGraph,1,'descend');
nodeCstPerGraph2 = zeros(nodeCnt,graphCnt);
for ref = 1:graphCnt
    nodeCstPerGraph2(IX(1:inCnt,ref),ref) = 1;
end
nodeCstPerGraph = nodeCstPerGraph2;