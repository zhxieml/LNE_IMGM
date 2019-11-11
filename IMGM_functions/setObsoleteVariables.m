target.config.bUnaryEnable = 0;%bUnaryEnable=1 use point-wise unary similarity, otherwise not
target.config.bEdgeEnable = 1;%bEdgeEnable=1 use edge-wise 2nd-order similarity, otherwise not
target.config.bSaveRandom = 0;% not to save the random graphs. it is used when reproducing the exact same results
target.config.bCreateRandom = 1;% if set to 0, will read the random graphs from files
target.config.affinityBiDir = 1;% used for mOpt
target.config.bPermute = 1;% set the ground truth by identity matrix, other choices are not used in demo
% some old parameters are used for debug and other tests, less relevant to the algorithm
% by default set to 0 for target.config.useCstInlier (affinity-driven), no use in formal test type, only useful in massiveOutlier mode 
% target.config.useCstInlier = 0;% set to 1 if use consistency to mask inliers, otherwise use affinity metrics, see Sec 3.5 in PAMI paper
% random graph test, note not the random point set test as used in mpm
% testType: different modes for tests, e.g. formal (Fig.3&4), case(Fig.1), iter(Fig.2), massOutlier(Fig.5&6)
% target.config.testType = 'formal';% for logic simplicity, this demo code involves only formal case for random graphs Fig.3, another is massOutlier.