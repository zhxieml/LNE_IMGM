target.config.bUnaryEnable = 0;%bUnaryEnable=1 use point-wise unary similarity, otherwise not
target.config.bEdgeEnable = 1;%bEdgeEnable=1 use edge-wise 2nd-order similarity, otherwise not
target.config.bSaveRandom = 0;% not to save the random graphs. it is used when reproducing the exact same results
target.config.bCreateRandom = 1;% if set to 0, will read the random graphs from files
target.config.affinityBiDir = 1;% used for mOpt
target.config.bPermute = 1;% set the ground truth by identity matrix, other choices are not used in demo
