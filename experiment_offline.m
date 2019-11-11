%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;

setPlotColor;
setObsoleteVariables;% some old parameters are used for debug and other tests, less relevant to the algorithm
% by default set to 0 for target.config.useCstInlier (affinity-driven), no use in formal test type, only useful in massiveOutlier mode 
% target.config.useCstInlier = 0;% set to 1 if use consistency to mask inliers, otherwise use affinity metrics, see Sec 3.5 in PAMI paper
% random graph test, note not the random point set test as used in mpm
% testType: different modes for tests, e.g. formal (Fig.3&4), case(Fig.1), iter(Fig.2), massOutlier(Fig.5&6)
% target.config.testType = 'formal';% for logic simplicity, this demo code involves only formal case for random graphs Fig.3, another is massOutlier.
target.config.testType = 'formal';% massOutlier
algpar = setPairwiseSolver();
mpmAlgPar = setMPMAlgPar;

varyMinGrhCnt=20; varyMaxGrhCnt=50; grhTestCnt = 10;% 
target.config.database = 'synthetic';% only synthetic test is allowed here
target.config.Sacle_2D = 0.05;
iterRange = 6;
graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
algNameSepSpace = '                    ';
algSet.algNameSet = {'cao_pc','cao_pc_raw','imgm_d','imgm_r','tbimgm_cao','tbimgm_cao_pc','tbimgm_cao_uc','tbimgm_qm','tbimgm_matchALS', 'tbimgm_cao_pc_ordered'};
algSet.algEnable =  [ 1,        0,           0,       0,          0,           1,              0,              0,            0,                     1];
algSet.algColor = {cao_pcClr,cao_pc_rawClr,imgm_dClr,imgm_rClr, tbimgm_caoClr, tbimgm_cao_pcClr, tbimgm_cao_ucClr, tbimgm_qmClr, tbimgm_matchALSClr, tbimgm_cao_pcClr};
algSet.algLineStyle = {'--','--','-','--','-','--','-','--','-'};
algSet.algMarker = {'.','.','.','.','.','.','.','.','.'};
target.config.bGraphMatch = 1;
target.config.inCntType = 'all';% set 'all' for "only a few outlier case", e.g. Fig.1&2&3&4
target.config.category = 'deform';%'deform','outlier','density','complete'
switch target.config.category
    case 'deform' % same setting with 5th row in Table 1 in the PAMI paper 
        nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.15;
        target.config.density = .9;target.config.complete = 1;
        graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
    case 'outlier' % same setting with 6th row in Table 1 in the PAMI paper 
        nInlier = 6;target.config.nOutlier = 4;target.config.deform = 0;
        target.config.density = 1;target.config.complete = 1;
        graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
    case 'density' % same setting with 7th row in Table 1 in the PAMI paper 
        nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.0;
        target.config.density = 0.5;target.config.complete = 1;
        graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
    case 'complete' % same setting with 8th row in Table 1 in the PAMI paper 
        nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.05;
        target.config.density = 1;target.config.complete = 0.1;     
end

graphStep = 1;
baseGraphCnt = graphMinCnt;
graphRange = baseGraphCnt:graphStep:graphMaxCnt-graphStep;
target.config.baseGraphCnt = baseGraphCnt;
target.config.graphRange = graphRange;
target.config.initConstWeight = .2; % initial weight for consitency regularizer, suggest 0.2-0.25
target.config.constStep = 1.1;% inflate parameter, suggest 1.1-1.2
target.config.constWeightMax = 1;
target.config.constIterImmune = 2; % in early iterations, not involve consistency, suggest 1-3
target.config.edgeAffinityWeight = 1;% in random graphs, only edge affinity is used, angle is meaningless
target.config.angleAffinityWeight = 1 - target.config.edgeAffinityWeight;
target.config.selectNodeMask = 1:1:nInlier+target.config.nOutlier;
target.config.selectGraphMask{1} = 1:graphMaxCnt;
paraCnt = length(graphRange);


[~,cao_pcIdx] = ismember('cao_pc',algSet.algNameSet);
[~,cao_pc_rawIdx] = ismember('cao_pc_raw',algSet.algNameSet);
[~,imgm_dIdx] = ismember('imgm_d',algSet.algNameSet);
[~,imgm_rIdx] = ismember('imgm_r',algSet.algNameSet);
[~,tbimgm_caoIdx] = ismember('tbimgm_cao', algSet.algNameSet);
[~,tbimgm_cao_pcIdx] = ismember('tbimgm_cao_pc', algSet.algNameSet);
[~,tbimgm_cao_ucIdx] = ismember('tbimgm_cao_uc', algSet.algNameSet);
[~,tbimgm_qmIdx] = ismember('tbimgm_qm', algSet.algNameSet);
[~,tbimgm_matchALSIdx] = ismember('tbimgm_matchALS', algSet.algNameSet);
[~,tbimgm_cao_pc_orderedIdx] = ismember('tbimgm_cao_pc_ordered', algSet.algNameSet);


algCnt = length(algSet.algNameSet);
X=cell(algCnt,1);
target.config.nodeCnt = length(target.config.selectNodeMask);
% target.config.graphCnt = min(max(graphRange),length(target.config.selectGraphMask{1}));
target.config.graphCnt = 60;
nodeCnt = target.config.nodeCnt;
graphCnt = target.config.graphCnt;

% paraCnt: iterate over graph #
% algCnt: iterate over algorithms
% testCnt: iterate over tests
timAve = zeros(algCnt,testCnt);timAveFull = zeros(algCnt);
accAve = zeros(algCnt,testCnt);accAveFull = zeros(algCnt);
scrAve = zeros(algCnt,testCnt);scrAveFull = zeros(algCnt);
conPairAve = zeros(algCnt,testCnt);conPairAveFull = zeros(algCnt);
countPairAve = zeros(algCnt,testCnt);countPairAveFull = zeros(algCnt);
accStd = zeros(algCnt,testCnt);accStdFull = zeros(algCnt);
scrStd = zeros(algCnt,testCnt);scrStdFull = zeros(algCnt);
conPairStd = zeros(algCnt,testCnt);conPairStdFull = zeros(paraCnt,algCnt);

% conMatPairGraph = cell(paraCnt,algCnt,testCnt);
% accMatPairGraph = cell(paraCnt,algCnt,testCnt);
% scrMatPairGraph = cell(paraCnt,algCnt,testCnt);
acc = cell(1, algCnt);
scr = cell(1, algCnt);
con = cell(1, algCnt);
% prevMatching = cell(1, algCnt);
Matching = cell(1, algCnt);
matTmp = cell(1, algCnt);

fidPerf = fopen('results.csv','w');
fprintf(fidPerf, 'testType,bGraphMatch,unaryFeat,edgeFeat,inCntType,database,category,testCnt,iter#,total node#,outlier#,complete,density,deform,graph#,alg#,scale,edgeWeight,initConstWeight,consStep,constWeightMax,iterImmune\n');
fprintf(fidPerf, '%s,%d,%d,%d,%s,%s,%s,%d,%d,%d,%d,%.2f,%.2f,%.2f,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n',...
    target.config.testType,target.config.bGraphMatch,target.config.bUnaryEnable,target.config.bEdgeEnable,target.config.inCntType,target.config.database,target.config.category,...
    testCnt,max(iterRange),nodeCnt,target.config.nOutlier,target.config.complete,target.config.density,target.config.deform,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,target.config.edgeAffinityWeight,...
    target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
fprintf('testType=%s, bGraphMatch=%d, unaryEnable=%d, edgeEnable=%d, inCntType=%s, database=%s, category=%s, iter#=%d, test#=%d, node#=%d, outlier#=%d,complete=%.2f, density=%.2f, deform=%.2f, graph#=%d, alg#=%d, edgeWeight=%.2f, scale=%.2f, initW=%.2f, stepW=%.2f, maxW=%.2f,iterImmune=%d\n',...
    target.config.testType,target.config.bGraphMatch,target.config.bUnaryEnable,target.config.bEdgeEnable,target.config.inCntType,target.config.database,target.config.category,max(iterRange),testCnt,...
    nodeCnt,target.config.nOutlier,target.config.complete,target.config.density,target.config.deform,graphCnt,sum(algSet.algEnable),...
    target.config.edgeAffinityWeight,target.config.Sacle_2D,target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
fprintf('\n');fprintf(fidPerf,'\n');
inlierAcc = zeros(paraCnt,testCnt);
estErr = zeros(testCnt,1);


%% Run the test
%%%%%%%%%%%%%%%%%%  experiment 1 different incremental algorithms %%%%%%%%%%%%%%%%%%%%%

for testk = 1:testCnt
    fprintf('Run test in round %d\n', testk);

    affinity = generateRandomAffinity(nInlier,testk); 
    affinity.GT = repmat(eye(nodeCnt,nodeCnt),graphCnt,graphCnt);

    % rrwm pairwise match, once for all graph pairs
    rawMat = generatePairAssignment(algpar,nodeCnt,graphCnt,testk);

    switch target.config.inCntType
        case 'exact' % already known, used in Fig.5 and top two rows in Fig.6
             target.config.inCnt = nodeCnt - target.config.nOutlier;
        case 'all' % in case of few outliers, used in Fig.1,2,3,4
             target.config.inCnt = nodeCnt;
        case 'spec' % specified by user, used in the bottom row of Fig.6
            target.config.inCnt = specNodeCnt;
    end
    
    scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
    scrDenomMatInCntGT = cal_pair_graph_inlier_score(affinity.GT,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
    scrDenomCurrent = max(max(scrDenomMatInCnt));

    %%%%%%%%%%%% calculate the multigraph matching with cao_pc %%%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(cao_pcIdx)
        tStart = tic;
        Matching{cao_pcIdx} = CAO(rawMat, nodeCnt, graphCnt, iterRange, scrDenomCurrent, 'pair', 1);
        tEnd = toc(tStart);
        
        acc{cao_pcIdx} = cal_pair_graph_accuracy(Matching{cao_pcIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{cao_pcIdx} = cal_pair_graph_score(Matching{cao_pcIdx},affinity.GT,nodeCnt,graphCnt);
        con{cao_pcIdx} = cal_pair_graph_consistency(Matching{cao_pcIdx},nodeCnt,graphCnt,0);
        accAve(cao_pcIdx, testk) = mean(acc{cao_pcIdx}(:));
        scrAve(cao_pcIdx, testk) = mean(scr{cao_pcIdx}(:));
        conPairAve(cao_pcIdx, testk) = mean(con{cao_pcIdx}(:));
        timAve(cao_pcIdx, testk) = tEnd;
        countPairAve(cao_pcIdx, testk) = graphCnt;
    end
    
    %%%%%%%%%%%% calculate the multigraph matching with imgm_d %%%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(imgm_dIdx)     
        param.iterMax = iterRange;
        param.visualization = 0;
        param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
        
        tStart = tic;
        Matching{imgm_dIdx} = IMGM_offline(rawMat,nodeCnt,graphCnt,param);
        tEnd = toc(tStart);
        
        acc{imgm_dIdx} = cal_pair_graph_accuracy(Matching{imgm_dIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{imgm_dIdx} = cal_pair_graph_score(Matching{imgm_dIdx},affinity.GT,nodeCnt,graphCnt);
        con{imgm_dIdx} = cal_pair_graph_consistency(Matching{imgm_dIdx},nodeCnt,graphCnt,0);
        accAve(imgm_dIdx, testk) = mean(acc{imgm_dIdx}(:));
        scrAve(imgm_dIdx, testk) = mean(scr{imgm_dIdx}(:));
        conPairAve(imgm_dIdx, testk) = mean(con{imgm_dIdx}(:));
        timAve(imgm_dIdx, testk) = tEnd;
        countPairAve(imgm_dIdx, testk) = graphCnt;
    end
    
    %%%%%%%%%%%% calculate the multigraph matching with imgm_r %%%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(imgm_rIdx)     
        param.iterMax = iterRange;
        param.visualization = 0;
        param.method = 3; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
        
        tStart = tic;
        Matching{imgm_rIdx} = IMGM_offline(rawMat,nodeCnt,graphCnt,param);
        tEnd = toc(tStart);
        
        acc{imgm_rIdx} = cal_pair_graph_accuracy(Matching{imgm_rIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{imgm_rIdx} = cal_pair_graph_score(Matching{imgm_rIdx},affinity.GT,nodeCnt,graphCnt);
        con{imgm_rIdx} = cal_pair_graph_consistency(Matching{imgm_rIdx},nodeCnt,graphCnt,0);
        accAve(imgm_rIdx, testk) = mean(acc{imgm_rIdx}(:));
        scrAve(imgm_rIdx, testk) = mean(scr{imgm_rIdx}(:));
        conPairAve(imgm_rIdx, testk) = mean(con{imgm_rIdx}(:));
        timAve(imgm_rIdx, testk) = tEnd;
        countPairAve(imgm_rIdx, testk) = graphCnt;
    end
    
    %%%%%%%%%%%% calculate the multigraph matching with tbimgm_cao_pc %%%%%%%%%%%%%%%%%%%%%
    % param of tbimgm
    param.bVerbose = 0;
    param.maxNumSearch = 20;
    if algSet.algEnable(tbimgm_cao_pcIdx)
        % param for tbimgm_cao_pc
        param.subMethodParam.name = 'CAO';
        param.subMethodParam.useCstDecay = 1;
        param.subMethodParam.cstDecay  = 0.7;
        param.subMethodParam.useWeightedDecay  = 0;
        param.subMethodParam.iterMax = 5;
        param.subMethodParam.scrDenom = scrDenomCurrent;
        param.subMethodParam.optType = 'pair';
        param.subMethodParam.useCstInlier = 1;
        
        % param for adaptive order
        param.adaptive_order = false;
        
        tStart = tic;
        Matching{tbimgm_cao_pcIdx} = TBIMGM_offline(rawMat,nodeCnt,graphCnt,param);
        tEnd = toc(tStart);
        
        acc{tbimgm_cao_pcIdx} = cal_pair_graph_accuracy(Matching{tbimgm_cao_pcIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{tbimgm_cao_pcIdx} = cal_pair_graph_score(Matching{tbimgm_cao_pcIdx},affinity.GT,nodeCnt,graphCnt);
        con{tbimgm_cao_pcIdx} = cal_pair_graph_consistency(Matching{tbimgm_cao_pcIdx},nodeCnt,graphCnt,0);
        accAve(tbimgm_cao_pcIdx, testk) = mean(acc{tbimgm_cao_pcIdx}(:));
        scrAve(tbimgm_cao_pcIdx, testk) = mean(scr{tbimgm_cao_pcIdx}(:));
        conPairAve(tbimgm_cao_pcIdx, testk) = mean(con{tbimgm_cao_pcIdx}(:));
        timAve(tbimgm_cao_pcIdx, testk) = tEnd;
        countPairAve(tbimgm_cao_pcIdx, testk) = graphCnt;
    end
    
    %%%%%%%%%%%% calculate the multigraph matching with tbimgm_cao_ordered %%%%%%%%%%%%%%%%%%%%%
    % param of tbimgm
    param.bVerbose = 0;
    param.maxNumSearch = 20;
    if algSet.algEnable(tbimgm_cao_pc_orderedIdx)
        % param for tbimgm_cao_pc
        param.subMethodParam.name = 'CAO';
        param.subMethodParam.useCstDecay = 1;
        param.subMethodParam.cstDecay  = 0.7;
        param.subMethodParam.useWeightedDecay  = 0;
        param.subMethodParam.iterMax = 5;
        param.subMethodParam.scrDenom = scrDenomCurrent;
        param.subMethodParam.optType = 'pair';
        param.subMethodParam.useCstInlier = 1;
        
        % param for adaptive order
        param.adaptive_order = true;
        
        tStart = tic;
        Matching{tbimgm_cao_pc_orderedIdx} = TBIMGM_offline(rawMat,nodeCnt,graphCnt,param);
        tEnd = toc(tStart);
        
        acc{tbimgm_cao_pc_orderedIdx} = cal_pair_graph_accuracy(Matching{tbimgm_cao_pc_orderedIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{tbimgm_cao_pc_orderedIdx} = cal_pair_graph_score(Matching{tbimgm_cao_pc_orderedIdx},affinity.GT,nodeCnt,graphCnt);
        con{tbimgm_cao_pc_orderedIdx} = cal_pair_graph_consistency(Matching{tbimgm_cao_pc_orderedIdx},nodeCnt,graphCnt,0);
        accAve(tbimgm_cao_pc_orderedIdx, testk) = mean(acc{tbimgm_cao_pc_orderedIdx}(:));
        scrAve(tbimgm_cao_pc_orderedIdx, testk) = mean(scr{tbimgm_cao_pc_orderedIdx}(:));
        conPairAve(tbimgm_cao_pc_orderedIdx, testk) = mean(con{tbimgm_cao_pc_orderedIdx}(:));
        timAve(tbimgm_cao_pc_orderedIdx, testk) = tEnd;
        countPairAve(tbimgm_cao_pc_orderedIdx, testk) = graphCnt;
    end
    
    %%%%%%%%%%%% calculate the multigraph matching with tbimgm_cao_uc %%%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(tbimgm_cao_ucIdx)
        % param for tbimgm_cao_pc
        param.subMethodParam.name = 'CAO';
        param.subMethodParam.useCstDecay = 1;
        param.subMethodParam.cstDecay  = 0.7;
        param.subMethodParam.useWeightedDecay  = 0;
        param.subMethodParam.iterMax = 5;
        param.subMethodParam.scrDenom = scrDenomCurrent;
        param.subMethodParam.optType = 'unary';
        param.subMethodParam.useCstInlier = 1;
        
        tStart = tic;
        Matching{tbimgm_cao_ucIdx} = TBIMGM_offline(rawMat,nodeCnt,graphCnt,param);
        tEnd = toc(tStart);
        
        acc{tbimgm_cao_ucIdx} = cal_pair_graph_accuracy(Matching{tbimgm_cao_ucIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{tbimgm_cao_ucIdx} = cal_pair_graph_score(Matching{tbimgm_cao_ucIdx},affinity.GT,nodeCnt,graphCnt);
        con{tbimgm_cao_ucIdx} = cal_pair_graph_consistency(Matching{tbimgm_cao_ucIdx},nodeCnt,graphCnt,0);
        accAve(tbimgm_cao_ucIdx, testk) = mean(acc{tbimgm_cao_ucIdx}(:));
        scrAve(tbimgm_cao_ucIdx, testk) = mean(scr{tbimgm_cao_ucIdx}(:));
        conPairAve(tbimgm_cao_ucIdx, testk) = mean(con{tbimgm_cao_ucIdx}(:));
        timAve(tbimgm_cao_ucIdx, testk) = tEnd;
        countPairAve(tbimgm_cao_ucIdx, testk) = graphCnt;
    end
    
    fprintf('test in round %d\n',testk);
    fprintf('%-18s%-18s%-18s%-18s%-18s%-18s\n','field\alg', 'accuracy', 'score', 'consistency', 'time', 'numPairMatch');
    for alg = find(algSet.algEnable)
        fprintf('%-18s%-18f%-18f%-18f%-18f%-18d\n\n',algSet.algNameSet{alg}, accAve(alg, testk), scrAve(alg, testk), conPairAve(alg, testk), timAve(alg, testk), countPairAve(alg, testk));
    end

end % for testCnt

for algk = 1:algCnt
    timAveFull(algk) = mean(timAve(algk,:));
    accAveFull(algk) = mean(accAve(algk,:));
    scrAveFull(algk) = mean(scrAve(algk,:));
    conPairAveFull(algk) = mean(conPairAve(algk,:));
    countPairAveFull(algk) = mean(countPairAve(algk,:));
end
    
i = 0;
savePath = sprintf('./res/off/offline_exp_%d.mat', i);

while isfile(savePath)
    i = i + 1;
    savePath = sprintf('./res/off/offline_exp_%d.mat', i);
end
save(savePath, 'target', 'algSet', 'accAveFull', 'scrAveFull', 'conPairAveFull', 'timAveFull', 'countPairAveFull');

