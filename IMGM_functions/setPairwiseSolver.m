function algpar = setPairwiseSolver
algpar.iterMax1 = 10;%meta算法的内部迭代次数 但不包括fgm 10和20的结果一样 都是提前收敛的不用到迭代最大次数
algpar.bDisc = 1;%采取离散化处理再做大矩阵
algpar.iccvIterMax = 1;
% algpar.uptRatio = 1;
algpar.algMethod = 'RRWM_old';
