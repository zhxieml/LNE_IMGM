function [ matchScore, hyperGraph, clusterMember, interClusterGraph, assignMat] = IMGM(method, threshold, viewInfo)
    %%  Incremental Multi Graph Matching
    set_param_GM;
    nGraph = length(viewInfo); % # of graphs
    
    interClusterGraph = zeros(nGraph, nGraph);
    % interClusterGraph(n, p) means connection graph from cluster n

    % preallocated memory for clusters
    clusterMember = zeros(nGraph, nGraph);
    clusterSize = zeros(1, nGraph);
    clusterCenter = zeros(1, nGraph);

    % initialize cluster
    clusterMember(1, 1) = 1;
    clusterCenter(1) = 1;
    clusterSize(1) = 1;
    nCluster = 1;

    hyperGraph = zeros(nGraph, nGraph,'int32');

    np = size(viewInfo(1).desc, 2); 
    % here is a restriction: two graph must have equal number of Keyt points

    assignMat = cell(nGraph, nGraph);
    affinityMatrix = cell(nGraph, nGraph);
    % need a way to save the storage
    
    matchScore = zeros(nGraph, nGraph);

    % calculate transF for further use
    transF = getTrans(np);
    
    % calculate group1 group2 for further use
    [X, Y] = meshgrid(1:np, 1:np);
    matchlist = [X(:), Y(:)]';
    [group1, group2] = make_group12(matchlist);

    % preallocated memory for depth first search
    stk = zeros(1, nGraph, 'uint32');
    Xfu = zeros(np, np, nGraph);
    Ku = zeros(np*np);
    visited = zeros(1, nGraph, 'logical');

    for iNewGraph = 2:nGraph
        % fprintf("graph NO.%d: %s\n",iNewGraph, viewInfo(iNewGraph).fileName);
        maxScore = 0;
        bestCluster = -1;
        %% match members of all cluster center
        for iCluster = 1:nCluster
            % TODO: need to eliminate the second loop in order to improve the speed
            iCenter = clusterCenter(iCluster);
            cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;
            cdata.view(1) = viewInfo(iCenter);  cdata.view(2) = viewInfo(iNewGraph);
            restmp = wrapper_matchingBasic(method, cdata);
            %TODO: need to refine pairgraph matching process, wrapper_matchingBasic is too complex
            assignMat{iCenter, iNewGraph} = restmp.X;
            % suppose the assignment matrix X(iCenter -> iNewGraph) is
            %                               [x11, x12, x13, ..., x1n
            %                                x21, x22, .............
            %    X(iCenter -> iNewGraph) =   .......................
            %                                .......................
            %                                xn1, .............. xnn]
            % note that restmp.X = [x11, x12, ..., x1n, x21, x22, ....x2n, ...., xnn]
            % after reshape restmp.X becomes
                            %               [x11, x21, x31, ..., xn1
                            %                x12, x22, .............
                            %    restmp.X =  .......................
                            %                .......................
                            %                x1n, .............. xnn]
            % which is the transpose of X(iCenter -> iNewGraph)
            % due to the special matlab property, we don't do transpose to restmp, because
            % if we transpose it, 
            % we need to transpose it again in order to get restmp.X
            % whose form is convenient to calculate the affinity and consistency in future works.

            matchScore(iCenter, iNewGraph) = restmp.score_GM;
            matchScore(iNewGraph, iCenter) = restmp.score_GM;
            if maxScore < matchScore(iCenter, iNewGraph)
                bestCluster = iCluster;
                maxScore = matchScore(iCenter, iNewGraph);
            end
        end
        %fprintf("max score: %.2f\n", maxScore);
        if maxScore > threshold
            %% match new graph with all members of best cluster
            bestScore = maxScore; % initialize
            bestGraph = clusterCenter(bestCluster);
            for iGraph = clusterMember(bestCluster, 1:clusterSize(bestCluster))
                cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;
                cdata.view(1) = viewInfo(iGraph);  cdata.view(2) = viewInfo(iNewGraph);
                restmp = wrapper_matchingBasic(method, cdata);
                % all I want is just a score
                % assignMat(:, iGraph, iNewGraph) = restmp.X;
                affinityMatrix{iGraph, iNewGraph} = restmp.affinityMatrix;
                affinityMatrix{iNewGraph, iGraph} = transF*(restmp.affinityMatrix*transF);
                if bestScore < restmp.score_GM
                    bestGraph = iGraph;
                    bestScore = restmp.score_GM;
                end
            end

            u = bestGraph;
            r = iNewGraph;
            %connect new graph with best match

            hyperGraph(u, r) = 1;
            hyperGraph(r, u) = 1;

            clusterMember(bestCluster, clusterSize(bestCluster)+1) = iNewGraph;
            clusterSize(bestCluster) = clusterSize(bestCluster)+1;

            pointSet = clusterMember(bestCluster, 1:clusterSize(bestCluster));
            % pointSet is the Set of graphs of bestCluster
            MST = hyperGraph(pointSet, pointSet);
            % Strictly speaking, this is not a Maximum Spanning Tree.
            ir = clusterSize(bestCluster);
            iu = find(MST(ir, :), 1);

            %% solve sub-problems as a standard QAP problem 
            %% Depth first search for composed affinity matrix
            nCand = clusterSize(bestCluster);
            Ku(:) = 0;
            Xfu(:) = 0;
            stk(:) = 0;
            visited(:) = 0;
            top = 1;
            visited(ir) = true;
            stk(top) = iu;
            Xfu(:, :, iu) = eye(np);
            while top >= 1
                visited(stk(top)) = true;
                adjSet = MST(stk(top), :) & (~visited(1:nCand));
                iff = find(adjSet, 1);
                if isempty(iff)
                    top = top - 1;
                else
                    % t is the father point, f is the adjecent point
                    t = pointSet(stk(top));
                    f = pointSet(iff);
                    Xft = reshape(assignMat{t, f}, np, np);
                    Xfu(:, :, iff) = Xft*Xfu(:, :, stk(top)); % Xfu = Xft * Xtu
                    
                    % compute affinity matrix
                    Ffu = kron(eye(np), Xfu(:, :, iff));
                    Krf = affinityMatrix{r, f};
                    Ku = Ku + Ffu'*(Krf*Ffu); 
                    % update stack
                    top = top + 1;
                    stk(top) = iff;
                end
            end


            K = affinityMatrix{u, r} + transF*(Ku*transF); % compose affinity matrix
            %% solve composed QAP problem
            Xur_raw = RRWM(K, group1, group2);
            Xur = greedyMapping(Xur_raw, group1, group2);
            assignMat{u, r} = Xur;
            matchScore(u, r) = Xur'*(affinityMatrix{u, r}*Xur);
            % reverse tranform
            
            assignMat{r, u}= transF*Xur;
            matchScore(r, u) = matchScore(u, r);
            % TODO: update cluster center, most edge connected
            [~, I]= max(sum(MST, 1));
            clusterCenter(bestCluster) = pointSet(I(1));
            
        else
            %% create new cluster
            nCluster = nCluster + 1;
            clusterMember(nCluster, 1) = iNewGraph;
            clusterCenter(nCluster) = iNewGraph;
            clusterSize(nCluster) = 1;
            bestCluster = nCluster;
        end
        %% post consistency method;
        for iCluster = 1:nCluster
            if iCluster == bestCluster
                continue
            end
            st = interClusterGraph(bestCluster, iCluster);
            ed = interClusterGraph(iCluster, bestCluster);
            if st == 0 || matchScore(iNewGraph, clusterCenter(iCluster)) > matchScore(st, ed)
                interClusterGraph(bestCluster, iCluster) = iNewGraph;
                interClusterGraph(iCluster, bestCluster) = clusterCenter(iCluster);
            end
        end
        % after IMGM, need a Kruscal method to find MST between clusters
    end
end