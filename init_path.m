addpath('functions')
addpath('svm_struct_mex');

%% GM algorithms
% You should edit "set_alg.m" to add other graph matching algorithms
addpath('gm_algorithm\RRWM\'); % Reweighted Random Walks Matching by Cho et al. ECCV2010
addpath('gm_algorithm\SM\'); % Spectral Matching by Leordeanu and Hebert. ICCV2005

run('..\vlfeat-0.9.21\toolbox\vl_setup');