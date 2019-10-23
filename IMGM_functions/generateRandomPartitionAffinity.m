function affinity = generateRandomPartitionAffinity(halfNodeCnt, param)
% function generating a graph with two partitions for CutMatch
% param.saveGraph: if save the graph
% param.deform: deformation factor 
% param.gap: gap between two partitions
    fLength = 128; % feature length on nodes
    delta = 5;
    isLength = 1; % generating distance map for graph cut 
    basePartitionSet = randn(2, halfNodeCnt); % generate base graph nodes
    refPartitionSet = basePartitionSet + param.deform*randn(2, halfNodeCnt); % generate reference graph nodes
    minL = max(basePartitionSet(1,:))-min(refPartitionSet(1,:));
    refPartitionSet = refPartitionSet + repmat([minL+param.gap;0],[1 halfNodeCnt]);
    GraphSet = [basePartitionSet, refPartitionSet]; % generate graph nodes
    tmpZeroMatch = zeros(halfNodeCnt,halfNodeCnt);
    tmpEyeMatch = eye(halfNodeCnt);
    affinity.GT = [tmpZeroMatch,tmpEyeMatch;tmpEyeMatch,tmpZeroMatch];
    if param.saveGraph
        affinity.graph = GraphSet;
    else
        affinity.graph = [];
    end
    if param.isTri
        tri = delaunayTriangulation(GraphSet');
        triEdge = edges(tri);
    end
    
    GraphSet = GraphSet';
    [X, Y] = meshgrid(1:2*halfNodeCnt, 1:2*halfNodeCnt);
    LL = [X(:),Y(:)];
    G = GraphSet(LL(:,1),:)-GraphSet(LL(:,2),:);
    G = sqrt(G(:,1).^2+G(:,2).^2);
    G = reshape(G, [2*halfNodeCnt 2*halfNodeCnt]);
    if param.isTri
        tri = delaunayTriangulation(GraphSet);
        triEdge = edges(tri);
        mask = zeros([2*halfNodeCnt,2*halfNodeCnt]);
        for k = 1:size(triEdge,1)
            mask(triEdge(k,1),triEdge(k,2))=1;
        end
        mask = mask + mask';
    end
    affinity.edge = mask; % edge matrix of graph
    % a = 1;
    
    [X, Y]=meshgrid(1:(2*halfNodeCnt)^2,1:(2*halfNodeCnt)^2);
    X = X(:);  Y = Y(:);
    G = G(:);
    
    aff = exp(-(G(X) - G(Y)).^2/param.delta);
    affinity.aff = reshape(aff, [((2*halfNodeCnt)^2),((2*halfNodeCnt)^2)]); % initialize affinity matrix
    if param.isTri
        maskE = zeros([((2*halfNodeCnt)^2),((2*halfNodeCnt)^2)]);
        for i = 1:2*halfNodeCnt
            for j = 1:2*halfNodeCnt;
                for k = 1:2*halfNodeCnt
                    for l = 1:2*halfNodeCnt
                        if mask(i,j)==1 && mask(k,l)==1 && (i~=k || j~=l)
                            maskE(sub2ind([2*halfNodeCnt, 2*halfNodeCnt],i,k),sub2ind([2*halfNodeCnt, 2*halfNodeCnt],j,l))=1;
                        end
                    end
                end
            end
        end
        affinity.aff = affinity.aff.*maskE;
        affinity.tri = tri.ConnectivityList;
    end
    
    % generate node features
    if param.isNodeFeature == 1
        affinity.nodeFeature = zeros(fLength, 2*halfNodeCnt);
        affinity.nodeFeature(:,1:halfNodeCnt)=rand(fLength,halfNodeCnt);
        affinity.nodeFeature(:,halfNodeCnt+1:end)=affinity.nodeFeature(:,1:halfNodeCnt)+param.featureDeform*randn(fLength,halfNodeCnt); % add noise to the ref graph node features
    end
    % disturb one half-graph
    if exist('param.variationProb')
        perm = randperm(halfNodeCnt);
        for i = 1:2:halfNodeCnt*param.variationProb
            tmp = affinity.nodeFeature(:,halfNodeCnt+perm(i));
            affinity.nodeFeature(:,halfNodeCnt+perm(i))=affinity.nodeFeature(:,halfNodeCnt+perm(i+1));
            affinity.nodeFeature(:,halfNodeCnt+perm(i+1))=tmp;
        end
    end
    % disturb the whole graph
    if exist('param.variationConProb');
        perm = randperm(2*halfNodeCnt);
        for i = 1:2:2*halfNodeCnt*param.variationConProb
            tmp = affinity.nodeFeature(:,perm(i));
            affinity.nodeFeature(:,perm(i))=affinity.nodeFeature(:,perm(i+1));
            affinity.nodeFeature(:,perm(i+1))=tmp;
        end
    end
    if param.isFeatureWeight
        affinity.featureWeight = zeros(2*halfNodeCnt,2*halfNodeCnt);
        for i = 1:2*halfNodeCnt
            for j=1:2*halfNodeCnt
                if affinity.edge(i,j)==1
                    affinity.featureWeight(i,j)=exp(-(norm(affinity.nodeFeature(:,i)-affinity.nodeFeature(:,j)))^2/delta); % kernel function: similarity
                    if isLength
                        affinity.featureWeight(i,j) = affinity.featureWeight(i,j) + exp(-norm(GraphSet(i,:)-GraphSet(j,:))^2/0.5);
                    end
                    % affinity.featureWeight(i,j)=2*(1/(1+exp(-abs(0.01*norm(affinity.nodeFeature(:,i)-affinity.nodeFeature(:,j),1))))-0.5); % sigmoid function: dissimilarity
                end
            end
        end
    end
end