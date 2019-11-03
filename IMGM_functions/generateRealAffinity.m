function [affinity] = generateRealAffinity(testk)
    % generate affinity matrix
    global target;
    set_param_GM;
    affinity.BiDir = target.config.affinityBiDir;
    affinity.edgeAffinityWeight = target.config.edgeAffinityWeight;
    affinity.angleAffinityWeight = target.config.angleAffinityWeight;
    affinityDir = target.config.affinityDir;
    nInlier = target.config.nInlier;
    nOutlier = target.config.nOutlier;
    featDir = target.config.featDir;
    gtDir = target.config.gtDir;
    cls = target.config.class;
    complete = target.config.complete;
    nodeCnt = nOutlier + nInlier;
    graphCnt = target.config.graphCnt;
    weight = target.config.weight;
    Sacle_2D = target.config.Sacle_2D;
    affinity.graphCnt = graphCnt;
    affinity.nodeCnt = nodeCnt;
    affinity.EG = cell(graphCnt,1);
    M = getTrans(nodeCnt);
    if complete<1
        target.pairwiseMask{testk} = generateRandomCompleteMask(nodeCnt,graphCnt,target.config.complete);
    else % fully complete, hence all are ones
        target.pairwiseMask{testk} = ones(graphCnt*nodeCnt,graphCnt*nodeCnt);
    end
    
    savePath = fullfile(affinityDir, cls, sprintf('in%02d_out%02d.mat', nInlier, nOutlier));
    if isfile(savePath)
        load(savePath, 'affinity');
    else
        listOfFeatFile = dir(fullfile(featDir, cls, '*.mat'));
        for gc = 1:graphCnt
            affinity.nP{gc} = nodeCnt;
            feat = load(fullfile(listOfFeatFile(gc).folder, listOfFeatFile(gc).name));
            pointSet = feat.frames(1:2, :);
            pointSet = pointSet';
            [X, Y] = meshgrid(1:nodeCnt, 1:nodeCnt);
            LL = [X(:), Y(:)];
            G = pointSet(LL(:,1),:)-pointSet(LL(:,2),:); %edge# \times 2
            edgeLen = sqrt(G(:, 1).^2 + G(:, 2).^2); %edge# \times 1
            mask = (edgeLen ~= 0);
            edgeAngle = zeros(size(edgeLen));
            edgeAngle(mask) = acos(G(mask,1)./edgeLen(mask));
            edgeAngle = reshape(edgeAngle, nodeCnt, nodeCnt);
            edgeLen = edgeLen ./ max(edgeLen);
            edgeLen = reshape(edgeLen, [nodeCnt nodeCnt]);
            affinity.edge{gc} = edgeLen';
            affinity.edgeAngle{gc} = edgeAngle'; % nodeCnt \time nodeCnt
            affinity.pointFeat{gc} = feat.descriptors;
            affinity.nE{gc} = sum(sum(~isnan(edgeLen)));
            [r,c]=find(~isnan(edgeLen));
            affinity.EG{gc}=[r,c]';% size is 2 \times Data.nE{1} (i.e. number of edges), the fist row is one ending point, the second row is the other
            % incidence matrix
            affinity.G{gc} = zeros(nodeCnt,affinity.nE{gc});
            for c = 1 : affinity.nE{gc}
                affinity.G{gc}(affinity.EG{gc}(:, c), c) = 1;
            end
            % augumented adjacency
            affinity.H{gc} = [affinity.G{gc}, eye(nodeCnt)];
        end
        affinity.GT = zeros(graphCnt*nodeCnt);
        mask = eye(nodeCnt^2, 'logical');
        for x = 1:graphCnt          
            nameXgt = listOfFeatFile(x).name;
            xview = (x-1)*nodeCnt+1:(x-1)*nodeCnt+nInlier;
            for y = x+1:graphCnt
                nameYgt = listOfFeatFile(y).name;
                yview = (y-1)*nodeCnt+1:(y-1)*nodeCnt+nInlier;
                if target.config.bUnaryEnable
                    KU = exp(-conDst(affinity.pointFeat{x}, affinity.pointFeat{y},0) / Sacle_2D);
                    KU = KU';
                    KU = diag(KU(:));
                else
                    KU = zeros(nodeCnt^2);
                end
                if target.config.bEdgeEnable
                    A = kron(affinity.edge{x}, ones(nodeCnt));
                    B = kron(ones(nodeCnt), affinity.edge{y});
                    KE = exp(-((A-B).^2)/Sacle_2D);
                    KE(mask) = 0;
                else
                    KE = zeros(nodeCnt^2);
                end
                if target.config.bAngleEnable
                    A = kron(affinity.edgeAngle{x}, ones(nodeCnt));
                    B = kron(ones(nodeCnt), affinity.edgeAngle{y});
                    KA = abs(A - B);
                else
                    KA = zeros(nodeCnt^2);
                end
                affinity.K{x, y} = weight(1)*KU + weight(2)* KE + weight(3)*KA;
                affinity.K{y, x} = M*(affinity.K{x, y}*M);
                gtPath = fullfile(gtDir, cls, sprintf("%s_%s_%s.mat", cls, nameXgt(end-7:end-4), nameYgt(end-7:end-4)));
                gt = load(gtPath);
                affinity.GT(xview, yview) = gt.groundTruth.assign;
            end
        end
        affinity.GT = affinity.GT + affinity.GT' + eye(nodeCnt*graphCnt);
        if ~isfolder(fullfile(affinityDir, cls))
            mkdir(fullfile(affinityDir, cls));
        end
        save(savePath, 'affinity');
    end
    
end