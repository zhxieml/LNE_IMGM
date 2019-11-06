addpath('functions');
addpath('IMGM_functions');
addpath('gm_algorithm\RRWM\'); % Reweighted Random Walks Matching by Cho et al. ECCV2010
addpath(genpath('multigraph_algorithms'));
% addpath('multigraph_algorithms\');
% addpath('multigraph_algorithms\CAO');
% addpath('multigraph_algorithms\DPMC');
% addpath('multigraph_algorithms\QuickMatch');
% addpath('multigraph_algorithms\MatchALS');
% addpath('multigraph_algorithms\Helpers');
run('dpp\startup.m');
run('..\vlfeat-0.9.21\toolbox\vl_setup');