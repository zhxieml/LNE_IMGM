%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;

setPlotColor;
algpar = setPairwiseSolver();

target.config.bUnaryEnable = 0;%bUnaryEnable=1 use point-wise unary similarity, otherwise not
target.config.bEdgeEnable = 1;%bEdgeEnable=1 use edge-wise 2nd-order similarity, otherwise not
target.config.bAngleEnable = 1;
target.config.weight = [0, 1, 0]; % weight for Unary, Edge and Angle
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


varyMinGrhCnt=30; varyMaxGrhCnt=50 ; grhTestCnt = 1;% 
target.config.Sacle_2D = 0.05;
target.config.database = 'WILLOW-Object-Class';% only synthetic test is allowed here

target.config.dataDir = 'D:\zchen\data\WILLOW-ObjectClass-dataset';
target.config.imgDir=[target.config.dataDir '/WILLOW-ObjectClass'];
target.config.annoDir=[target.config.dataDir '/WILLOW-ObjectClass'];

target.config.gtDir = [target.config.dataDir '/ground_truth'];
target.config.resDir = './res';
target.config.tmpDir = './tmp';
target.config.class = 'Duck';
target.config.category = 'outlier';%'deform','outlier'
target.config.distRatioTrue = 0.15;
target.config.testType = 'formal';% massOutlier
target.config.constIterImmune = 2;% immune from consistency regularization, only use affinity in earl
target.config.constStep = 1.1;% the inflate parameter, e.g. 1.05-1.1
target.config.initConstWeight = 0.2;% initial consistency regularizer weight, e.g 0.2-0.25
target.config.constWeightMax = 1;% the upperbound, always set to 1
switch target.config.category
    case 'deform' % same setting with 5th row in Table 1 in the PAMI paper 
        target.config.nInlier = 10;
        target.config.nOutlier = 0;
        target.config.featDir = [target.config.dataDir '\feature_0'];
        target.config.affinityDir=[target.config.dataDir '\affinity_0'];
        target.config.complete = 1;
    case 'outlier'
        target.config.nInlier = 10;
        target.config.nOutlier = 4;
        target.config.featDir = [target.config.dataDir '\feature_4'];
        target.config.affinityDir=[target.config.dataDir '\affinity_4'];
        target.config.complete = 1;
end
iterRange = 6;
graphMinCnt = varyMinGrhCnt;
testCnt = grhTestCnt;

switch target.config.class
    case 'Duck'
        graphMaxCnt = min(50, varyMaxGrhCnt);
    case 'Motorbike'
        graphMaxCnt = min(40, varyMaxGrhCnt);
    case 'Car'
        graphMaxCnt = min(40, varyMaxGrhCnt);
    case 'Face'
        graphMaxCnt = min(109, varyMaxGrhCnt);
    case 'Winebottle'
        graphMaxCnt = min(66, varyMaxGrhCnt);
end


algNameSepSpace = '                    ';
algSet.algNameSet = {'cao_pc','cao_pc_raw','imgm_d','imgm_r','tbimgm_cao','tbimgm_cao_pc','tbimgm_cao_uc','tbimgm_qm','tbimgm_matchALS'};
algSet.algEnable =  [ 1,        1,           1,       1,       1,           1,              0,              0,          0];
algSet.algColor = {cao_pcClr,cao_pc_rawClr,imgm_dClr,imgm_rClr, tbimgm_caoClr, tbimgm_cao_pcClr, tbimgm_cao_ucClr, tbimgm_qmClr, tbimgm_matchALSClr};
algSet.algLineStyle = {'--','--','-','--','-','--','-','--','-'};
algSet.algMarker = {'.','.','.','.','.','.','.','.','.'};

[~,cao_pcIdx] = ismember('cao_pc',algSet.algNameSet);
[~,cao_pc_rawIdx] = ismember('cao_pc_raw',algSet.algNameSet);
[~,imgm_dIdx] = ismember('imgm_d',algSet.algNameSet);
[~,imgm_rIdx] = ismember('imgm_r',algSet.algNameSet);
[~,tbimgm_caoIdx] = ismember('tbimgm_cao', algSet.algNameSet);
[~,tbimgm_cao_pcIdx] = ismember('tbimgm_cao_pc', algSet.algNameSet);
[~,tbimgm_cao_ucIdx] = ismember('tbimgm_cao_uc', algSet.algNameSet);
[~,tbimgm_qmIdx] = ismember('tbimgm_qm', algSet.algNameSet);
[~,tbimgm_matchALSIdx] = ismember('tbimgm_matchALS', algSet.algNameSet);

algCnt = length(algSet.algNameSet);

graphStep = 1;
baseGraphCnt = graphMinCnt;
graphRange = baseGraphCnt:graphStep:graphMaxCnt-graphStep;
paraCnt = length(graphRange);
nInlier = target.config.nInlier;
nOutlier = target.config.nOutlier;
target.config.totalCnt = graphMaxCnt;
target.config.nodeCnt = nInlier + nOutlier;
target.config.graphCnt = graphRange(end) + graphStep;
nodeCnt = target.config.nodeCnt;
graphCnt = target.config.graphCnt;

target.config.initConstWeight = .2; % initial weight for consitency regularizer, suggest 0.2-0.25
target.config.constStep = 1.1;% inflate parameter, suggest 1.1-1.2
target.config.constWeightMax = 1;
target.config.constIterImmune = 2; % in early iterations, not involve consistency, suggest 1-3
target.config.edgeAffinityWeight = 0.9;% in random graphs, only edge affinity is used, angle is meaningless
target.config.angleAffinityWeight = 1 - target.config.edgeAffinityWeight;
target.config.selectNodeMask = 1:1:nInlier+target.config.nOutlier;
target.config.selectGraphMask{1} = 1:graphMaxCnt;
target.config.connect = 'fc';

% paraCnt: iterate over graph #
% algCnt: iterate over algorithms
% testCnt: iterate over tests
timAve = zeros(paraCnt,algCnt,testCnt);timAveFull = zeros(paraCnt,algCnt);
accAve = zeros(paraCnt,algCnt,testCnt);accAveFull = zeros(paraCnt,algCnt);
scrAve = zeros(paraCnt,algCnt,testCnt);scrAveFull = zeros(paraCnt,algCnt);
conPairAve = zeros(paraCnt,algCnt,testCnt);conPairAveFull = zeros(paraCnt,algCnt);
countPairAve = zeros(paraCnt,algCnt,testCnt);countPairAveFull = zeros(paraCnt,algCnt);

acc = cell(1, algCnt);
scr = cell(1, algCnt);
con = cell(1, algCnt);
prevMatching = cell(1, algCnt);
increMatching = cell(1, algCnt);
matTmp = cell(1, algCnt);


fidPerf = fopen('results.csv','w');
fprintf(fidPerf, 'testType,database,category,testCnt,total node#,outlier#,graph#,alg#,scale,initConstWeight,consStep,constWeightMax,iterImmune\n');
fprintf(fidPerf, '%s,%d,%d,%d,%s,%s,%s,%d,%d,%d,%d,%.2f,%.2f,%.2f,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n',...
    target.config.testType, target.config.database, target.config.category,...
    testCnt,nodeCnt,target.config.nOutlier,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,...
    target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
fprintf('testType=%s, database=%s, category=%s, test#=%d, node#=%d, outlier#=%d,graph#=%d, alg#=%d, scale=%.2f, initW=%.2f, stepW=%.2f, maxW=%.2f,iterImmune=%d\n',...
    target.config.testType,target.config.database,target.config.category,testCnt,...
    nodeCnt,target.config.nOutlier,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,...
    target.config.initConstWeight, target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
fprintf('\n');fprintf(fidPerf,'\n');
%% Run the test
%%%%%%%%%%%%%%%%%%  experiment 1 different incremental algorithms %%%%%%%%%%%%%%%%%%%%%

for testk = 1:testCnt
    fprintf('Run test in round %d\n', testk);

    affinity = generateRealAffinity(testk);
    rawMat = generatePairAssignment(algpar,nodeCnt,graphCnt,testk);
    sigma = 0;
    % rrwm pairwise match, once for all graph pairs
    target.config.inCnt = nodeCnt - target.config.nOutlier;
    target.pairwiseMask = cell(1);
	target.pairwiseMask{1} = ones(graphCnt*nodeCnt,graphCnt*nodeCnt);
    scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
    scrDenomMatInCntGT = cal_pair_graph_inlier_score(affinity.GT,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
    
    for parak = 1:paraCnt 

        param.n = nodeCnt; 
        param.N = baseGraphCnt + (parak-1)*graphStep; % 20
        param.graphStep = graphStep;
        scrDenomCurrent = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
        baseMat = CAO(rawMat(1:nodeCnt*param.N,1:nodeCnt*param.N), nodeCnt, param.N, iterRange,scrDenomCurrent, 'pair',1);
        
        %%%%%%%%%%%% calculate the incremental matching with cao_pc_raw %%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(cao_pc_rawIdx)
            tStart = tic;
            increMatching{cao_pc_rawIdx} = CAO(rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)), nodeCnt, param.N+graphStep, iterRange, scrDenomCurrent, 'pair', 1);
            tEnd = toc(tStart);
            
            acc{cao_pc_rawIdx} = cal_pair_graph_accuracy(increMatching{cao_pc_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{cao_pc_rawIdx} = cal_pair_graph_score(increMatching{cao_pc_rawIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{cao_pc_rawIdx} = cal_pair_graph_consistency(increMatching{cao_pc_rawIdx},nodeCnt,param.N+graphStep,0);
            accAve(parak, cao_pc_rawIdx, testk) = mean(acc{cao_pc_rawIdx}(:));
            scrAve(parak, cao_pc_rawIdx, testk) = mean(scr{cao_pc_rawIdx}(:));
            conPairAve(parak, cao_pc_rawIdx, testk) = mean(con{cao_pc_rawIdx}(:));
            timAve(parak, cao_pc_rawIdx, testk) = tEnd;
            countPairAve(parak, cao_pc_rawIdx, testk) = (param.N+graphStep);
        end

        %%%%%%%%%%% calculate the incremental matching with cao_pc %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % previous matching: caoPrevMatching
        if algSet.algEnable(cao_pcIdx)
            if parak == 1
                prevMatching{cao_pcIdx} = baseMat;
            end
            matTmp{cao_pcIdx} = rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep));
            matTmp{cao_pcIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{cao_pcIdx};
            tStart = tic;
            increMatching{cao_pcIdx} = CAO(matTmp{cao_pcIdx}, nodeCnt, param.N+graphStep , iterRange, scrDenomCurrent, 'pair',1);
            tEnd = toc(tStart);
            prevMatching{cao_pcIdx} = increMatching{cao_pcIdx};

            acc{cao_pcIdx} = cal_pair_graph_accuracy(increMatching{cao_pcIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{cao_pcIdx} = cal_pair_graph_score(increMatching{cao_pcIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{cao_pcIdx} = cal_pair_graph_consistency(increMatching{cao_pcIdx},nodeCnt,param.N+graphStep,0);
            
            accAve(parak, cao_pcIdx, testk) = mean(acc{cao_pcIdx}(:));
            scrAve(parak, cao_pcIdx, testk) = mean(scr{cao_pcIdx}(:));
            conPairAve(parak, cao_pcIdx, testk) = mean(con{cao_pcIdx}(:));
            timAve(parak, cao_pcIdx, testk) = tEnd;
            countPairAve(parak, cao_pcIdx, testk) = (param.N+graphStep);
        end

        %%%%%%%%%%% calculate the incremental matching with imgm_d %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_dIdx)
            if parak == 1
                prevMatching{imgm_dIdx} = baseMat;
            end
            matTmp{imgm_dIdx} = rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep));
            matTmp{imgm_dIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{imgm_dIdx};
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{imgm_dIdx},affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{imgm_dIdx},nodeCnt,param.N+graphStep,0);
            

            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            % param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
            param.iterMax = iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_dIdx} = IMGM_old(simAP, matTmp{imgm_dIdx}, param);
            tEnd = toc(tStart);
            prevMatching{imgm_dIdx} = increMatching{imgm_dIdx};
            
            acc{imgm_dIdx} = cal_pair_graph_accuracy(increMatching{imgm_dIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{imgm_dIdx} = cal_pair_graph_score(increMatching{imgm_dIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{imgm_dIdx} = cal_pair_graph_consistency(increMatching{imgm_dIdx},nodeCnt,param.N+graphStep,0);
            
            accAve(parak, imgm_dIdx, testk) = mean(acc{imgm_dIdx}(:));
            scrAve(parak, imgm_dIdx, testk) = mean(scr{imgm_dIdx}(:));
            conPairAve(parak, imgm_dIdx, testk) = mean(con{imgm_dIdx}(:));
            timAve(parak, imgm_dIdx, testk) = tEnd;
            countPairAve(parak, imgm_dIdx, testk) = (param.N+graphStep);
        end

        %%%%%%%%%%% calculate the incremental matching with imgm_r %%%%%%%%%%%%%%%%%%%%%%%        
        if algSet.algEnable(imgm_rIdx)
            if parak == 1
                prevMatching{imgm_rIdx} = baseMat;
            end
            matTmp{imgm_rIdx} = rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep));
            matTmp{imgm_rIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{imgm_rIdx};
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{imgm_rIdx},affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{imgm_rIdx},nodeCnt,param.N+graphStep,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            % param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
            param.iterMax = iterRange;
            param.visualization = 0;
            param.method = 3; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_rIdx} = IMGM_old(simAP, matTmp{imgm_rIdx}, param);
            tEnd = toc(tStart);
            prevMatching{imgm_rIdx} = increMatching{imgm_rIdx};
            
            acc{imgm_rIdx} = cal_pair_graph_accuracy(increMatching{imgm_rIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{imgm_rIdx} = cal_pair_graph_score(increMatching{imgm_rIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{imgm_rIdx} = cal_pair_graph_consistency(increMatching{imgm_rIdx},nodeCnt,param.N+graphStep,0);
            
            accAve(parak, imgm_rIdx, testk) = mean(acc{imgm_rIdx}(:));
            scrAve(parak, imgm_rIdx, testk) = mean(scr{imgm_rIdx}(:));
            conPairAve(parak, imgm_rIdx, testk) = mean(con{imgm_rIdx}(:));
            timAve(parak, imgm_rIdx, testk) = tEnd;
            countPairAve(parak, imgm_rIdx, testk) = (param.N+graphStep);
        end

        %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_pc %%%%%%%%%%%%%%%%%%%%%%%%%
        % param of tbimgm
        param.bVerbose = 0;
        param.maxNumSearch = int32(0.7*param.N);
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
            % previous matching: prevMatching{tbimgm_cao_pcIdx}
            if parak == 1
                prevMatching{tbimgm_cao_pcIdx} = baseMat;
            end
            matTmp{tbimgm_cao_pcIdx} = rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep));
            matTmp{tbimgm_cao_pcIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_cao_pcIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_cao_pcIdx},affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_cao_pcIdx},nodeCnt,param.N+graphStep,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));

            tStart = tic;
            [increMatching{tbimgm_cao_pcIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_cao_pcIdx}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao_pcIdx} = increMatching{tbimgm_cao_pcIdx};
            
            acc{tbimgm_cao_pcIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_cao_pcIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{tbimgm_cao_pcIdx} = cal_pair_graph_score(increMatching{tbimgm_cao_pcIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{tbimgm_cao_pcIdx} = cal_pair_graph_consistency(increMatching{tbimgm_cao_pcIdx},nodeCnt,param.N+graphStep,0);
            
            accAve(parak, tbimgm_cao_pcIdx, testk) = mean(acc{tbimgm_cao_pcIdx}(:));
            scrAve(parak, tbimgm_cao_pcIdx, testk) = mean(scr{tbimgm_cao_pcIdx}(:));
            conPairAve(parak, tbimgm_cao_pcIdx, testk) = mean(con{tbimgm_cao_pcIdx}(:));
            timAve(parak, tbimgm_cao_pcIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_cao_pcIdx, testk) = numPairMatch;
        end

       %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_caoIdx)
            % param for tbimgm_caoIdx
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.useCstDecay = 1;
            param.subMethodParam.cstDecay  = 0.7;
            param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = 5;
            param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'afnty';
            param.subMethodParam.useCstInlier = 1;
            % previous matching: prevMatching{tbimgm_caoIdx}
            if parak == 1
                prevMatching{tbimgm_caoIdx} = baseMat;
            end
            matTmp{tbimgm_caoIdx} = rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep));
            matTmp{tbimgm_caoIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_caoIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_caoIdx},affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_caoIdx},nodeCnt,param.N+graphStep,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));

            tStart = tic;
            [increMatching{tbimgm_caoIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_caoIdx}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_caoIdx} = increMatching{tbimgm_caoIdx};
            
            acc{tbimgm_caoIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_caoIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{tbimgm_caoIdx} = cal_pair_graph_score(increMatching{tbimgm_caoIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{tbimgm_caoIdx} = cal_pair_graph_consistency(increMatching{tbimgm_caoIdx},nodeCnt,param.N+graphStep,0);
            
            accAve(parak, tbimgm_caoIdx, testk) = mean(acc{tbimgm_caoIdx}(:));
            scrAve(parak, tbimgm_caoIdx, testk) = mean(scr{tbimgm_caoIdx}(:));
            conPairAve(parak, tbimgm_caoIdx, testk) = mean(con{tbimgm_caoIdx}(:));
            timAve(parak, tbimgm_caoIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_caoIdx, testk) = numPairMatch;
        end

        %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_uc %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_cao_ucIdx)
            % param for tbimgm_cao_ucIdx
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.useCstDecay = 1;
            param.subMethodParam.cstDecay  = 0.7;
            param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = 5;
            param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'unary';
            param.subMethodParam.useCstInlier = 1;
            % previous matching: prevMatching{tbimgm_cao_ucIdx}
            if parak == 1
                prevMatching{tbimgm_cao_ucIdx} = baseMat;
            end
            matTmp{tbimgm_cao_ucIdx} = rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep));
            matTmp{tbimgm_cao_ucIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_cao_ucIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_cao_ucIdx},affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_cao_ucIdx},nodeCnt,param.N+graphStep,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));

            tStart = tic;
            [increMatching{tbimgm_cao_ucIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_cao_ucIdx}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao_ucIdx} = increMatching{tbimgm_cao_ucIdx};
            
            acc{tbimgm_cao_ucIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_cao_ucIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{tbimgm_cao_ucIdx} = cal_pair_graph_score(increMatching{tbimgm_cao_ucIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{tbimgm_cao_ucIdx} = cal_pair_graph_consistency(increMatching{tbimgm_cao_ucIdx},nodeCnt,param.N+graphStep,0);
            
            accAve(parak, tbimgm_cao_ucIdx, testk) = mean(acc{tbimgm_cao_ucIdx}(:));
            scrAve(parak, tbimgm_cao_ucIdx, testk) = mean(scr{tbimgm_cao_ucIdx}(:));
            conPairAve(parak, tbimgm_cao_ucIdx, testk) = mean(con{tbimgm_cao_ucIdx}(:));
            timAve(parak, tbimgm_cao_ucIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_cao_ucIdx, testk) = numPairMatch;
        end
        %%%%%%%%%%%% calculate the incremental matching with tbimgm_qm %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_qmIdx)
            % param for tbimgm_qmIdx
            param.subMethodParam.name = 'quickmatch';
            param.subMethodParam.rho = 0.25;
            param.subMethodParam.rho_edge = 0.5;
            % param.subMethodParam.kernel = @exp;
            param.subMethodParam.flag_low_memory = false;
            param.subMethodParam.kernel_support_radius = false;
            % previous matching: prevMatching{tbimgm_qmIdx}
            if parak == 1
                prevMatching{tbimgm_qmIdx} = baseMat;
            end
            matTmp{tbimgm_qmIdx} = rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep));
            matTmp{tbimgm_qmIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_qmIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_qmIdx},affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_qmIdx},nodeCnt,param.N+graphStep,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;

            tStart = tic;
            [increMatching{tbimgm_qmIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_qmIdx}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_qmIdx} = increMatching{tbimgm_qmIdx};
            
            acc{tbimgm_qmIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_qmIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{tbimgm_qmIdx} = cal_pair_graph_score(increMatching{tbimgm_qmIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{tbimgm_qmIdx} = cal_pair_graph_consistency(increMatching{tbimgm_qmIdx},nodeCnt,param.N+graphStep,0);
            
            accAve(parak, tbimgm_qmIdx, testk) = mean(acc{tbimgm_qmIdx}(:));
            scrAve(parak, tbimgm_qmIdx, testk) = mean(scr{tbimgm_qmIdx}(:));
            conPairAve(parak, tbimgm_qmIdx, testk) = mean(con{tbimgm_qmIdx}(:));
            timAve(parak, tbimgm_qmIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_qmIdx, testk) = numPairMatch;
        end     
        %%%%%%%%%%%% calculate the incremental matching with tbimgm_matchALS %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_matchALSIdx)
            % param for tbimgm_matchALSIdx
            param.subMethodParam.name = 'matchALS';

            % previous matching: prevMatching{tbimgm_matchALSIdx}
            if parak == 1
                prevMatching{tbimgm_matchALSIdx} = baseMat;
            end
            matTmp{tbimgm_matchALSIdx} = rawMat(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep));
            matTmp{tbimgm_matchALSIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_matchALSIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_matchALSIdx},affinity.GT(1:nodeCnt*(param.N+graphStep),1:nodeCnt*(param.N+graphStep)),nodeCnt,param.N+graphStep,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_matchALSIdx},nodeCnt,param.N+graphStep,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;

            tStart = tic;
            [increMatching{tbimgm_matchALSIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_matchALSIdx}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_matchALSIdx} = increMatching{tbimgm_matchALSIdx};
            
            acc{tbimgm_matchALSIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_matchALSIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+graphStep);
            scr{tbimgm_matchALSIdx} = cal_pair_graph_score(increMatching{tbimgm_matchALSIdx},affinity.GT,nodeCnt,param.N+graphStep);
            con{tbimgm_matchALSIdx} = cal_pair_graph_consistency(increMatching{tbimgm_matchALSIdx},nodeCnt,param.N+graphStep,0);
            
            accAve(parak, tbimgm_matchALSIdx, testk) = mean(acc{tbimgm_matchALSIdx}(:));
            scrAve(parak, tbimgm_matchALSIdx, testk) = mean(scr{tbimgm_matchALSIdx}(:));
            conPairAve(parak, tbimgm_matchALSIdx, testk) = mean(con{tbimgm_matchALSIdx}(:));
            timAve(parak, tbimgm_matchALSIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_matchALSIdx, testk) = numPairMatch;
        end    
        
        fprintf('test in round %d, Start from %d graphs, %d graphs incremented\n',testk, baseGraphCnt, parak*graphStep);
        fprintf('%-18s%-18s%-18s%-18s%-18s%-18s\n','field\alg', 'accuracy', 'score', 'consistency', 'time', 'numPairMatch');
        for alg = find(algSet.algEnable)
            fprintf('%-18s%-18f%-18f%-18f%-18f%-18d\n\n',algSet.algNameSet{alg}, accAve(parak, alg, testk), scrAve(parak, alg, testk), conPairAve(parak, alg, testk), timAve(parak, alg, testk), countPairAve(parak, alg, testk));
        end
    end

end % for paraCnt
    
% rename the algorithm names to be more friendly
for i=1:length(algSet.algNameSet)
    if strcmp(target.config.testType,'massOutlier')
        algSet.algNameSetDisplay{cao_Idx} = 'cao^{cst}';
        algSet.algNameSetDisplay{i} = strrep(algSet.algNameSet{i},'_s','^{sim}');
        algSet.algNameSetDisplay{i} = strrep(algSet.algNameSetDisplay{i},'o_','o-');
        algSet.algNameSetDisplay{i} = strrep(algSet.algNameSetDisplay{i},'c_','c^{cst}');
    else
        algSet.algNameSetDisplay{i} = strrep(algSet.algNameSet{i},'_','-');
        if algSet.algNameSetDisplay{i}(end)=='-'
            algSet.algNameSetDisplay{i}(end) = '*';
        end
    end
end

fprintf('--------------------------------------------------------------overall performance-------------------------------------------------------------------\n');
algNamePreSpace = '                          ';
fprintf(fidPerf,'overall mean\n');
fprintf(algNamePreSpace);
fprintf(fidPerf,',,');
for algk=1:algCnt
    if algSet.algEnable(algk)==0,continue;end
    fprintf([algSet.algNameSetDisplay{algk},algNameSepSpace]);
    fprintf(fidPerf,[algSet.algNameSetDisplay{algk},',,,,']);
end
fprintf('\n');fprintf(fidPerf,'\n');
fprintf('grh# itr#  ');fprintf(fidPerf,'grh#, itr#');
for algk=1:algCnt
    if algSet.algEnable(algk)==0,continue;end
    fprintf(' acc   scr   con   tim   pair');
    fprintf(fidPerf,', acc,  score, consis, time, pairmatch');
end
fprintf('\n');fprintf(fidPerf,'\n');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
for parak=1:paraCnt
    viewCnt=graphMinCnt + parak - 1;
    for algk = 1:algCnt
        timAveFull(parak,algk) = mean(timAve(parak,algk,:));
        accAveFull(parak,algk) = mean(accAve(parak,algk,:));
        scrAveFull(parak,algk) = mean(scrAve(parak,algk,:));
        conPairAveFull(parak,algk) = mean(conPairAve(parak,algk,:));
        countPairAveFull(parak, algk) = mean(countPairAve(parak, algk, :));
    end
    fprintf(' %02d,  %02d ',viewCnt,grhTestCnt);fprintf(fidPerf,' %02d,  %02d',viewCnt,grhTestCnt);
    for algk=1:algCnt
        if algSet.algEnable(algk)==0,continue;end
        fprintf('| %.3f %.3f %.3f %.3f %4d',accAveFull(parak,algk),scrAveFull(parak,algk),conPairAveFull(parak,algk),timAveFull(parak,algk), countPairAveFull(parak, algk));
        fprintf(fidPerf,', %.3f, %.3f, %.3f, %.3f, %4d',accAveFull(parak,algk),scrAveFull(parak,algk),conPairAveFull(parak,algk),timAveFull(parak,algk), countPairAveFull(parak, algk));% fprintf(cc,  score, consis
    end
    fprintf('\n');fprintf(fidPerf,'\n');
end

legendOff = 0;

current_time = datestr(datetime('now'));
savePath = sprintf('random_exp_all_inliers.mat');
save(savePath, 'target', 'algSet', 'accAveFull', 'scrAveFull', 'conPairAveFull', 'timAveFull', 'countPairAveFull');
ave.accuracy = accAveFull;
ave.score = scrAveFull;
ave.consistency = conPairAveFull;
ave.time = timAveFull;
ave.matchingNumber = countPairAveFull;

fields = fieldnames(ave);
for ifield = 1:length(fields)
    xtag='Arriving graph';ytag=[fields{ifield}, '(',target.config.category,')'];
    plotResult_new(legendOff,graphRange-baseGraphCnt+1, getfield(ave,fields{ifield}), algSet, xtag, ytag);
end
