function nodeAfnPerGraph= cal_node_affinity_mask(X,nodeCnt,graphCnt,inCnt)
global affinity
nodeAfnPerGraph = zeros(nodeCnt,graphCnt);
mask = zeros(nodeCnt^2,nodeCnt);% mask operation by 'and'
oneI = ones(nodeCnt,1);
for n = 1:nodeCnt
    mask((n-1)*nodeCnt+1:n*nodeCnt,n) = oneI;
end
for r = 1:graphCnt
    rscope = (r-1)*nodeCnt+1:r*nodeCnt;
    for i = 1:r-1
            K = affinity.K{i,r};%bidirectional affinity matrix set i,k<->k,i
            iscope = (i-1)*nodeCnt+1:i*nodeCnt;
            p = mat2vec(X(iscope,rscope));
            Kp = K*p;
            for n=1:nodeCnt
                q = p.*mask(:,n);
                nodeAfnPerGraph(n,r) = nodeAfnPerGraph(n,r) + q'*Kp;
            end
    end
    for i = r+1:graphCnt
            K = affinity.K{r,i};%bidirectional affinity matrix set i,k<->k,i
            iscope = (i-1)*nodeCnt+1:i*nodeCnt;
            pri = mat2vec(X(rscope,iscope));
            p = mat2vec(X(iscope,rscope));
            Kpri = K*pri;
            for n=1:nodeCnt
                q = p.*mask(:,n);
                q = vec2mat(q,nodeCnt,nodeCnt)';
                nodeAfnPerGraph(n,r) = nodeAfnPerGraph(n,r) + q(:)'*Kpri;
            end
    end
end
nodeAfnPerGraph = nodeAfnPerGraph/(graphCnt-1);
[~,IX] = sort(nodeAfnPerGraph,1,'descend');
nodeAfnPerGraph2 = zeros(nodeCnt,graphCnt);
for r = 1:graphCnt
    nodeAfnPerGraph2(IX(1:inCnt,r),r) = 1;%binarization, to pick the the top inCnt node as inliers
end
nodeAfnPerGraph = nodeAfnPerGraph2;
