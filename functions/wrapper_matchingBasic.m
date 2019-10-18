function [ res ] = wrapper_matchingBasic( method, cdata )

aparam.bOneToOne = 1;

% % command to call a matching function
% strFunc = ['feval(@' func2str(method.fhandle)];
% for j = 1:length(method.variable)
%     strFunc = [ strFunc ', cdata.' method.variable{j} ];
% end
% strFunc = [strFunc ')'];

% res.score_GM = -9999;   res.score_raw = -9999;

cdata.GT = [];
%% Find initial matches 
[ matchInfo ]= make_initialmatches_mcho2( cdata.view, cdata.mparam, 'verbose', false );

cand_matchlist = matchInfo.match;
cand_matchdist = matchInfo.dist;

[ uniq_feat1, tmp, new_feat1 ] = unique(cand_matchlist(1,:));    
[ uniq_feat2, tmp, new_feat2 ] = unique(cand_matchlist(2,:));
matchlist_tmp = [new_feat1 new_feat2]';

edgeAttr1 = computeEdgeAttr( cdata.view(1).frame(:,uniq_feat1) );
edgeAttr2 = computeEdgeAttr( cdata.view(2).frame(:,uniq_feat2) );

featBin1 = makeFeatBin(edgeAttr1);
featBin2 = makeFeatBin(edgeAttr2);

[ eSimVal ] = computeDotProdSimilarity_sym( int32(matchlist_tmp), featBin1, featBin2); % symmetric affinities

cdata.affinityMatrix = double(eSimVal);

%% Make the overlapping groups of initial matches
[ cdata.group1 cdata.group2 ] = make_group12(matchlist_tmp(1:2,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cdata.group1 = full(cdata.group1);
cdata.group2 = full(cdata.group2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ cdata.L12, cdata.E12, cdata.nP1, cdata.nP2 ] = updateL12E12( matchlist_tmp' );

unaryWeight = 10.0*ones(size(cand_matchdist));
vSimVal = max( 0, 0.8 - cand_matchdist );

cdata.affinityMatrix(1:(size(cdata.affinityMatrix,1)+1):end)...
    = unaryWeight(matchlist_tmp(1,:)) .* vSimVal;

cdata.affinityMatrix = cdata.affinityMatrix.*getNotConflictMatrix(cdata.group1, cdata.group2);
% cdata.affinityMatrix = cdata.affinityMatrix.*~getConflictMatrix(cdata.group1, cdata.group2);

%% ------------------ Perform graph matching
[ X_raw ] = RRWM(cdata.affinityMatrix, cdata.group1, cdata.group2);
% [ X_raw ] = feval(@RRWM, cdata.affinityMatrix, cdata.group1, cdata.group2);

if method.postDiscretization 
    X_sol = greedyMapping(X_raw, cdata.group1, cdata.group2);
    score_raw = X_sol'*cdata.affinityMatrix*X_sol;
else
    X_sol = X_raw;
end

score_GM = X_sol'*cdata.affinityMatrix*X_sol;

%% ----------------- Evaluate the solutions
res.X = X_sol;
res.X_raw = X_raw;
res.score_raw = score_raw;
res.score_GM = score_GM;
res.cand_matchlist = cand_matchlist;
res.affinityMatrix = cdata.affinityMatrix;
res.eSimVal = eSimVal;
res.vSimVal = vSimVal;
