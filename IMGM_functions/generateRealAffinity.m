function [affinity, rawMat] = generateRealAffinity
    % generate affinity matrix
    global target;
    set_param_GM;
    affinityDir = target.config.affinityDir;
    nInlier = target.config.nInlier;
    nOutlier = target.config.nOutlier;
    featDir = target.config.featDir;
    gtDir = target.config.gtDir;
    cls = target.config.class;
    nodeCnt = nOutlier + nInlier;
    graphCnt = target.config.graphCnt;
    affinityRatio = target.config.affinityRatio;
    
    savePath = fullfile(affinityDir, cls, sprintf('in%02d_out%02d.mat', nInlier, nOutlier));
    if isfile(savePath)
        load(savePath, 'affinity');
    else
        listOfFeatFile = dir(fullfile(featDir, cls, '*.mat'));
        for gc = 1:graphCnt
            feat = load(fullfile(listOfFeatFile(gc).folder, listOfFeatFile(gc).name);
            pointSet = feat.frames(1:2, :);
            [X, Y] = meshgrid(1:nodeCnt, 1:nodeCnt);
            LL = [X(:), Y(:)];
            edges = pointSet(LL(:,1),:)-pointSet(LL(:,2),:); %edge# \times 2
            edgeLen = sqrt(edges(:, 1).^2 + edges(:, 2).^2); %edge# \times 1
            mask = (edgeLen ~= 0);
            edgeAngle = zeros(size(edgeLen));
            edgeAngle(mask) = acos(edges(mask,1)./edgeLen(mask));

            affinity.edgeLen{gc} = edgeLen'/max(edgeLen);
            affinity.edgeAngle{gc} = edgeAngle';
            affinity.pointFeat{gc} = feat.desc;
        end
        for xview = 1:graphCnt
            if affinity.BiDir
                yviewSet = [1:xview-1, xview+1:graphCnt];
            else
                yviewSet = xview+1:graphCnt;
            end
            for yview = yviewSet
                if bUnaryEnable
                    affinity.KP{xview,yview} = exp(-conDst(affinity.pointFeat{xview}, affinity.pointFeat{yview},0) / Sacle_2D);
                else
                    affinity.KP{xview,yview} = zeros(nodeCnt,nodeCnt);
                end
                if bEdgeEnable
                    KQ = exp(-conDst(affinity.edgeLen{xview}, affinity.edgeLen{yview},0) / Sacle_2D);
                else
                    KQ = zeros(nodeCnt^2,nodeCnt^2);
                end
                if bAngleEnable
                    
                end
                KP = affinity.KP{xview, yview}(LL(:,1))
                affinity.K{xview, yview} = affinityRatio* KQ + (1-affinityRatio)*Kangle;
            end
        end
    end
    
end