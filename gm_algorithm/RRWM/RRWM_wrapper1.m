%% Implementaion of PageRank matching algortihm
function [ X score ] = RRWM_wrapper1( problem, AlgConfig )
% use asymmetrically normalized transition matrix

% note that RRWM_slck and RRWM_v4 differs only in using normalization indexes...
%[ score ] = RRWM_slack( problem.affinityMatrix, problem.group1, problem.group2 );
[ score ] = RRWM_v4( problem.affinityMatrix, problem.group1, problem.group2 );
%[ score ] = RRWM_v3( problem.affinityMatrix, problem.conflictMatrix1, problem.conflictMatrix2 );
% Discretize using the overlap matrix 
X = greedyMapping(score, problem.group1, problem.group2);
%X = greedyMapping(score, problem.conflictMatrix1 | problem.conflictMatrix2);

end