%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;

target.config.graphMinCnt=20; 
target.config.graphMaxCnt=40; 
target.config.testCnt = 20; % v
target.config.maxNumSearch = 21;
target.config.batchSize = 1;

target.config.database = "willow"; % "willow", "synthetic"
load_target_data;

% set algorithms
algNameSepSpace = '                    ';
algSet.algNameSet = {'cao_pc_inc', 'cao_pc_raw', 'cao_c_inc','cao_c_raw','imgm_d','imgm_r','tbimgm_cao_c','tbimgm_cao_pc','tbimgm_cao_cst', 'tbimgm_qm', 'tbimgm_matchALS', 'tbimgm_cao_c_ada'};
algSet.algEnable =  [ 0,            0,             1,           1,          1,       1,       1,              0,              0,               0,             0,             0];
algSet.algColor = { cao_pcClr, cao_pc_rawClr, cao_cClr,cao_c_rawClr,imgm_dClr,imgm_rClr, anc_imgmClr, tbimgm_cao_pcClr, tbimgm_cao_cstClr, tbimgm_qmClr, tbimgm_matchALSClr, tbimgm_cao_c_adaClr};
algSet.algLineStyle = {'--','--','-','--','-','--','-','--','-', '--', '-', '--', '-'};
algSet.algMarker = {'.','.','.','.','.','.', '.','.','.', '.', '.', '.', '.'};

[~,cao_pcIdx] = ismember('cao_pc_inc',algSet.algNameSet);
[~,cao_pc_rawIdx] = ismember('cao_pc_raw',algSet.algNameSet);
[~,cao_cIdx] = ismember('cao_c_inc',algSet.algNameSet);
[~,cao_c_rawIdx] = ismember('cao_c_raw',algSet.algNameSet);
[~,imgm_dIdx] = ismember('imgm_d',algSet.algNameSet);
[~,imgm_rIdx] = ismember('imgm_r',algSet.algNameSet);
[~,tbimgm_cao_cIdx] = ismember('tbimgm_cao_c', algSet.algNameSet);
[~,tbimgm_cao_pcIdx] = ismember('tbimgm_cao_pc', algSet.algNameSet);
[~,tbimgm_cao_cstIdx] = ismember('tbimgm_cao_cst', algSet.algNameSet);
[~,tbimgm_qmIdx] = ismember('tbimgm_qm', algSet.algNameSet);
[~,tbimgm_matchALSIdx] = ismember('tbimgm_matchALS', algSet.algNameSet);
[~,tbimgm_cao_cadaIdx] = ismember('tbimgm_cao_c_ada', algSet.algNameSet);

baseGraphCnt = target.config.graphMinCnt;
batchSize = target.config.batchSize;
target.config.graphRange = baseGraphCnt:batchSize:target.config.graphMaxCnt-batchSize;
target.config.baseGraphCnt = baseGraphCnt;
nInlier = target.config.nInlier;
nOutlier = target.config.nOutlier;
target.config.nodeCnt = nInlier + nOutlier;
target.config.graphCnt = target.config.graphRange(end) + batchSize;
nodeCnt = target.config.nodeCnt;
graphCnt = target.config.graphCnt;


% algorithms for affinity and CAO
target.config.Sacle_2D = 0.05;
target.config.iterRange = 6;
target.config.distRatioTrue = 0.15;
target.config.testType = 'all';% massOutlier
target.config.constStep = 1.05;% the inflate parameter, e.g. 1.05-1.1
target.config.constWeightMax = 1;% the upperbound, always set to 1
target.config.initConstWeight = 0.2; % initial weight for consitency regularizer, suggest 0.2-0.25
target.config.constIterImmune = 2; % in early iterations, not involve consistency, suggest 1-3
target.config.edgeAffinityWeight = 0.9;% in random graphs, only edge affinity is used, angle is meaningless
target.config.angleAffinityWeight = 1 - target.config.edgeAffinityWeight;
target.config.selectNodeMask = 1:1:nInlier+target.config.nOutlier;
target.config.selectGraphMask{1} = 1:target.config.graphMaxCnt;
target.config.connect = 'nfc';

% data for experiment
paraCnt = length(target.config.graphRange);% paraCnt: iterate over graph #
algCnt = length(algSet.algNameSet);% algCnt: iterate over algorithms
testCnt = target.config.testCnt;% testCnt: iterate over tests
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
fprintf(fidPerf, '%s,      %s,      %s,      %d,     %d,         %d,      %d,    %d,  %f,   %f,             %f,      %.2f,          %d\n',...
    target.config.testType, target.config.database, target.config.category,...
    testCnt,nodeCnt,target.config.nOutlier,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,...
    target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
fprintf('testType=%s, database=%s, category=%s, test#=%d, node#=%d, outlier#=%d,graph#=%d, alg#=%d, scale=%.2f, initW=%.2f, stepW=%.2f, maxW=%.2f,iterImmune=%d\n',...
    target.config.testType,target.config.database,target.config.category,testCnt,...
    nodeCnt,target.config.nOutlier,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,...
    target.config.initConstWeight, target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
fprintf('\n');fprintf(fidPerf,'\n');


for testk = 1:testCnt
    fprintf('Run test in round %d/%d\n', testk, testCnt);
    
    affinity = generateAffinity(testk);
    rawMat = generatePairAssignment(algpar,nodeCnt,graphCnt,testk);
    %debug
    accuracy_raw = cal_pair_graph_accuracy(rawMat, affinity.GT, target.config.nOutlier, nodeCnt, graphCnt);
    accuracy_ave = mean(accuracy_raw(:));
    fprintf("raw accuracy = %f\n", accuracy_ave);
    %debug
    sigma = 0;
    % rrwm pairwise match, once for all graph pairs
    switch target.config.inCntType
        case 'exact' % already known, used in Fig.5 and top two rows in Fig.6
            target.config.inCnt = nodeCnt - target.config.nOutlier;
        case 'all' % in case of few outliers, used in Fig.1,2,3,4
            target.config.inCnt = nodeCnt;
        case 'spec' % specified by user, used in the bottom row of Fig.6
            target.config.inCnt = specNodeCnt;
    end

    target.pairwiseMask = cell(1);
	target.pairwiseMask{1} = ones(graphCnt*nodeCnt,graphCnt*nodeCnt);
    scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
    scrDenomMatInCntGT = cal_pair_graph_inlier_score(affinity.GT,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
    param.n = nodeCnt; 
    param.batchSize = batchSize;
    scrDenomCurrent = max(max(scrDenomMatInCnt(1:baseGraphCnt,1:baseGraphCnt)));
    baseMat = CAO(rawMat(1:nodeCnt*baseGraphCnt,1:nodeCnt*baseGraphCnt), nodeCnt, baseGraphCnt, target.config.iterRange,scrDenomCurrent, 'pair',1);
    %baseMat = rawMat(1:nodeCnt*baseGraphCnt,1:nodeCnt*baseGraphCnt);
    for parak = 1:paraCnt 
        param.N = baseGraphCnt + (parak-1)*batchSize; % 20
        %%%%%%%%%%%% calculate the incremental matching with cao_c_raw %%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(cao_c_rawIdx)
            tStart = tic;
            increMatching{cao_c_rawIdx} = CAO(rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)), nodeCnt, param.N+batchSize, target.config.iterRange, scrDenomCurrent, 'exact', 1);
            tEnd = toc(tStart);
            
            acc{cao_c_rawIdx} = cal_pair_graph_accuracy(increMatching{cao_c_rawIdx}, affinity.GT, target.config.nOutlier, nodeCnt, param.N+batchSize);
            scr{cao_c_rawIdx} = cal_pair_graph_score(increMatching{cao_c_rawIdx},affinity.GT,nodeCnt,param.N+batchSize);
            %scr{cao_c_rawIdx} = cal_pair_graph_inlier_score(increMatching{cao_c_rawIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{cao_c_rawIdx} = cal_pair_graph_consistency(increMatching{cao_c_rawIdx},nodeCnt,param.N+batchSize,0);
            accAve(parak, cao_c_rawIdx, testk) = mean(acc{cao_c_rawIdx}(:));
            scrAve(parak, cao_c_rawIdx, testk) = mean(scr{cao_c_rawIdx}(:));
            conPairAve(parak, cao_c_rawIdx, testk) = mean(con{cao_c_rawIdx}(:));
            timAve(parak, cao_c_rawIdx, testk) = tEnd;
            countPairAve(parak, cao_c_rawIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end

        %%%%%%%%%%%% calculate the incremental matching with cao_pc_raw %%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(cao_pc_rawIdx)
            tStart = tic;
            increMatching{cao_pc_rawIdx} = CAO(rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)), nodeCnt, param.N+batchSize, target.config.iterRange, scrDenomCurrent, 'pair', 1);
            tEnd = toc(tStart);
            
            acc{cao_pc_rawIdx} = cal_pair_graph_accuracy(increMatching{cao_pc_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{cao_pc_rawIdx} = cal_pair_graph_score(increMatching{cao_pc_rawIdx},affinity.GT,nodeCnt,param.N+batchSize);
            %scr{cao_pc_rawIdx} = cal_pair_graph_inlier_score(increMatching{cao_pc_rawIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{cao_pc_rawIdx} = cal_pair_graph_consistency(increMatching{cao_pc_rawIdx},nodeCnt,param.N+batchSize,0);
            accAve(parak, cao_pc_rawIdx, testk) = mean(acc{cao_pc_rawIdx}(:));
            scrAve(parak, cao_pc_rawIdx, testk) = mean(scr{cao_pc_rawIdx}(:));
            conPairAve(parak, cao_pc_rawIdx, testk) = mean(con{cao_pc_rawIdx}(:));
            timAve(parak, cao_pc_rawIdx, testk) = tEnd;
            countPairAve(parak, cao_pc_rawIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end
        
        %%%%%%%%%%% calculate the incremental matching with cao_c_inc %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % previous matching: caoPrevMatching
        if algSet.algEnable(cao_cIdx)
            if parak == 1
                prevMatching{cao_cIdx} = baseMat;
            end
            
            matTmp{cao_cIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{cao_cIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{cao_cIdx};
            tStart = tic;
            increMatching{cao_cIdx} = CAO(matTmp{cao_cIdx}, nodeCnt, param.N+batchSize , target.config.iterRange, scrDenomCurrent, 'exact',1);
            tEnd = toc(tStart);
            prevMatching{cao_cIdx} = increMatching{cao_cIdx};

            acc{cao_cIdx} = cal_pair_graph_accuracy(increMatching{cao_cIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{cao_cIdx} = cal_pair_graph_score(increMatching{cao_cIdx},affinity.GT,nodeCnt,param.N+batchSize);
            %scr{cao_cIdx} = cal_pair_graph_inlier_score(increMatching{cao_cIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{cao_cIdx} = cal_pair_graph_consistency(increMatching{cao_cIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, cao_cIdx, testk) = mean(acc{cao_cIdx}(:));
            scrAve(parak, cao_cIdx, testk) = mean(scr{cao_cIdx}(:));
            conPairAve(parak, cao_cIdx, testk) = mean(con{cao_cIdx}(:));
            timAve(parak, cao_cIdx, testk) = tEnd;
            countPairAve(parak, cao_cIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end

        %%%%%%%%%%% calculate the incremental matching with cao_pc_inc %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % previous matching: caoPrevMatching
        if algSet.algEnable(cao_pcIdx)
            if parak == 1
                prevMatching{cao_pcIdx} = baseMat;
            end
            
            matTmp{cao_pcIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{cao_pcIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{cao_pcIdx};
            tStart = tic;
            increMatching{cao_pcIdx} = CAO(matTmp{cao_pcIdx}, nodeCnt, param.N+batchSize , target.config.iterRange, scrDenomCurrent, 'pair',1);
            tEnd = toc(tStart);
            prevMatching{cao_pcIdx} = increMatching{cao_pcIdx};

            acc{cao_pcIdx} = cal_pair_graph_accuracy(increMatching{cao_pcIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{cao_pcIdx} = cal_pair_graph_score(increMatching{cao_pcIdx},affinity.GT,nodeCnt,param.N+batchSize);
            %scr{cao_pcIdx} = cal_pair_graph_inlier_score(increMatching{cao_pcIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{cao_pcIdx} = cal_pair_graph_consistency(increMatching{cao_pcIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, cao_pcIdx, testk) = mean(acc{cao_pcIdx}(:));
            scrAve(parak, cao_pcIdx, testk) = mean(scr{cao_pcIdx}(:));
            conPairAve(parak, cao_pcIdx, testk) = mean(con{cao_pcIdx}(:));
            timAve(parak, cao_pcIdx, testk) = tEnd;
            countPairAve(parak, cao_pcIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end

        %%%%%%%%%%% calculate the incremental matching with imgm_d %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_dIdx)
            if parak == 1
                prevMatching{imgm_dIdx} = baseMat;
            end
            matTmp{imgm_dIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{imgm_dIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{imgm_dIdx};
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{imgm_dIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{imgm_dIdx},nodeCnt,param.N+batchSize,0);
            

            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            % param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_dIdx} = IMGM_local(affinity, simAP, matTmp{imgm_dIdx}, target, param);
                                                     % IMGM_local(affinity, affScore, rawMat, target, param)
            tEnd = toc(tStart);
            prevMatching{imgm_dIdx} = increMatching{imgm_dIdx};
            
            acc{imgm_dIdx} = cal_pair_graph_accuracy(increMatching{imgm_dIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{imgm_dIdx} = cal_pair_graph_score(increMatching{imgm_dIdx},affinity.GT,nodeCnt,param.N+batchSize);
            %scr{imgm_dIdx} = cal_pair_graph_inlier_score(increMatching{imgm_dIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{imgm_dIdx} = cal_pair_graph_consistency(increMatching{imgm_dIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, imgm_dIdx, testk) = mean(acc{imgm_dIdx}(:));
            scrAve(parak, imgm_dIdx, testk) = mean(scr{imgm_dIdx}(:));
            conPairAve(parak, imgm_dIdx, testk) = mean(con{imgm_dIdx}(:));
            timAve(parak, imgm_dIdx, testk) = tEnd;
            countPairAve(parak, imgm_dIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end

        %%%%%%%%%%% calculate the incremental matching with imgm_r %%%%%%%%%%%%%%%%%%%%%%%        
        if algSet.algEnable(imgm_rIdx)
            if parak == 1
                prevMatching{imgm_rIdx} = baseMat;
            end
            matTmp{imgm_rIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{imgm_rIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{imgm_rIdx};
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{imgm_rIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{imgm_rIdx},nodeCnt,param.N+batchSize,0);
            

            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            % param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 3; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_rIdx} = IMGM_local(affinity, simAP, matTmp{imgm_rIdx}, target, param);
                                                     % IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, param)
            tEnd = toc(tStart);
            prevMatching{imgm_rIdx} = increMatching{imgm_rIdx};
            
            acc{imgm_rIdx} = cal_pair_graph_accuracy(increMatching{imgm_rIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{imgm_rIdx} = cal_pair_graph_score(increMatching{imgm_rIdx},affinity.GT,nodeCnt,param.N+batchSize);
            %scr{imgm_rIdx} = cal_pair_graph_inlier_score(increMatching{imgm_rIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{imgm_rIdx} = cal_pair_graph_consistency(increMatching{imgm_rIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, imgm_rIdx, testk) = mean(acc{imgm_rIdx}(:));
            scrAve(parak, imgm_rIdx, testk) = mean(scr{imgm_rIdx}(:));
            conPairAve(parak, imgm_rIdx, testk) = mean(con{imgm_rIdx}(:));
            timAve(parak, imgm_rIdx, testk) = tEnd;
            countPairAve(parak, imgm_rIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end
       %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_c %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_cao_cIdx)
            % param for tbimgm_cao_cIdx
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.iterMax = target.config.iterRange;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N+batchSize,1:param.N+batchSize)));
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            % previous matching: prevMatching{tbimgm_cao_cIdx}
            if parak == 1
                prevMatching{tbimgm_cao_cIdx} = baseMat;
            end
            matTmp{tbimgm_cao_cIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{tbimgm_cao_cIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_cao_cIdx};

            tStart = tic;
            [increMatching{tbimgm_cao_cIdx}, numPairMatch] = ANC_IMGM_batch(affinity, target, matTmp{tbimgm_cao_cIdx}, nodeCnt, param.N, batchSize, 0, 'a', param);
                                                            % ANC_IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, opt, param)
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao_cIdx} = increMatching{tbimgm_cao_cIdx};
            
            acc{tbimgm_cao_cIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_cao_cIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{tbimgm_cao_cIdx} = cal_pair_graph_score(increMatching{tbimgm_cao_cIdx},affinity.GT,nodeCnt,param.N+batchSize);
            % scr{tbimgm_cao_cIdx} = cal_pair_graph_inlier_score(increMatching{tbimgm_cao_cIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{tbimgm_cao_cIdx} = cal_pair_graph_consistency(increMatching{tbimgm_cao_cIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, tbimgm_cao_cIdx, testk) = mean(acc{tbimgm_cao_cIdx}(:));
            scrAve(parak, tbimgm_cao_cIdx, testk) = mean(scr{tbimgm_cao_cIdx}(:));
            conPairAve(parak, tbimgm_cao_cIdx, testk) = mean(con{tbimgm_cao_cIdx}(:));
            timAve(parak, tbimgm_cao_cIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_cao_cIdx, testk) = numPairMatch;
        end
       %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_cada %%%%%%%%%%%%%%%%%%%%
       if algSet.algEnable(tbimgm_cao_cadaIdx)
            % param for tbimgm_cao_cadaIdx
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.iterMax = target.config.iterRange;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N+batchSize,1:param.N+batchSize)));
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            % previous matching: prevMatching{tbimgm_cao_cadaIdx}
            if parak == 1
                prevMatching{tbimgm_cao_cadaIdx} = baseMat;
            end
            matTmp{tbimgm_cao_cadaIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{tbimgm_cao_cadaIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_cao_cadaIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_cao_cadaIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_cao_cadaIdx},nodeCnt,param.N+batchSize,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));

            tStart = tic;
            [increMatching{tbimgm_cao_cadaIdx}, numPairMatch] = ANC_IMGM_batch(affinity, target, matTmp{tbimgm_cao_cadaIdx}, nodeCnt, param.N, batchSize, 1, 'a', param);
                                                            %ANC_IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, param)
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao_cadaIdx} = increMatching{tbimgm_cao_cadaIdx};
            
            acc{tbimgm_cao_cadaIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_cao_cadaIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{tbimgm_cao_cadaIdx} = cal_pair_graph_score(increMatching{tbimgm_cao_cadaIdx},affinity.GT,nodeCnt,param.N+batchSize);
            % scr{tbimgm_cao_cadaIdx} = cal_pair_graph_inlier_score(increMatching{tbimgm_cao_cadaIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{tbimgm_cao_cadaIdx} = cal_pair_graph_consistency(increMatching{tbimgm_cao_cadaIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, tbimgm_cao_cadaIdx, testk) = mean(acc{tbimgm_cao_cadaIdx}(:));
            scrAve(parak, tbimgm_cao_cadaIdx, testk) = mean(scr{tbimgm_cao_cadaIdx}(:));
            conPairAve(parak, tbimgm_cao_cadaIdx, testk) = mean(con{tbimgm_cao_cadaIdx}(:));
            timAve(parak, tbimgm_cao_cadaIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_cao_cadaIdx, testk) = numPairMatch;
        end

        %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_pc %%%%%%%%%%%%%%%%%%%%%%%%%
        % param of tbimgm
        if algSet.algEnable(tbimgm_cao_pcIdx)
            % param for tbimgm_cao_pc
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.useCstDecay = 1;
            param.subMethodParam.cstDecay  = 0.7;
            param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = target.config.iterRange;
            param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'pair';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            % previous matching: prevMatching{tbimgm_cao_pcIdx}
            if parak == 1
                prevMatching{tbimgm_cao_pcIdx} = baseMat;
            end
            matTmp{tbimgm_cao_pcIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{tbimgm_cao_pcIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_cao_pcIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_cao_pcIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_cao_pcIdx},nodeCnt,param.N+batchSize,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));

            tStart = tic;
            [increMatching{tbimgm_cao_pcIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_cao_pcIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao_pcIdx} = increMatching{tbimgm_cao_pcIdx};
            
            acc{tbimgm_cao_pcIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_cao_pcIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{tbimgm_cao_pcIdx} = cal_pair_graph_score(increMatching{tbimgm_cao_pcIdx},affinity.GT,nodeCnt,param.N+batchSize);
            % scr{tbimgm_cao_pcIdx} = cal_pair_graph_inlier_score(increMatching{tbimgm_cao_pcIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{tbimgm_cao_pcIdx} = cal_pair_graph_consistency(increMatching{tbimgm_cao_pcIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, tbimgm_cao_pcIdx, testk) = mean(acc{tbimgm_cao_pcIdx}(:));
            scrAve(parak, tbimgm_cao_pcIdx, testk) = mean(scr{tbimgm_cao_pcIdx}(:));
            conPairAve(parak, tbimgm_cao_pcIdx, testk) = mean(con{tbimgm_cao_pcIdx}(:));
            timAve(parak, tbimgm_cao_pcIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_cao_pcIdx, testk) = numPairMatch;
        end

        %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_cst %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_cao_cstIdx)
            % param for tbimgm_cao_cstIdx
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.useCstDecay = 1;
            param.subMethodParam.cstDecay  = 0.7;
            param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = target.config.iterRange;
            param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'cstcy';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            % previous matching: prevMatching{tbimgm_cao_cstIdx}
            if parak == 1
                prevMatching{tbimgm_cao_cstIdx} = baseMat;
            end
            matTmp{tbimgm_cao_cstIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{tbimgm_cao_cstIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_cao_cstIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_cao_cstIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_cao_cstIdx},nodeCnt,param.N+batchSize,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));

            tStart = tic;
            [increMatching{tbimgm_cao_cstIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_cao_cstIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao_cstIdx} = increMatching{tbimgm_cao_cstIdx};
            
            acc{tbimgm_cao_cstIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_cao_cstIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{tbimgm_cao_cstIdx} = cal_pair_graph_score(increMatching{tbimgm_cao_cstIdx},affinity.GT,nodeCnt,param.N+batchSize);
            % scr{tbimgm_cao_cstIdx} = cal_pair_graph_inlier_score(increMatching{tbimgm_cao_cstIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{tbimgm_cao_cstIdx} = cal_pair_graph_consistency(increMatching{tbimgm_cao_cstIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, tbimgm_cao_cstIdx, testk) = mean(acc{tbimgm_cao_cstIdx}(:));
            scrAve(parak, tbimgm_cao_cstIdx, testk) = mean(scr{tbimgm_cao_cstIdx}(:));
            conPairAve(parak, tbimgm_cao_cstIdx, testk) = mean(con{tbimgm_cao_cstIdx}(:));
            timAve(parak, tbimgm_cao_cstIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_cao_cstIdx, testk) = numPairMatch;
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
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            % previous matching: prevMatching{tbimgm_qmIdx}
            if parak == 1
                prevMatching{tbimgm_qmIdx} = baseMat;
            end
            matTmp{tbimgm_qmIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{tbimgm_qmIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_qmIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_qmIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_qmIdx},nodeCnt,param.N+batchSize,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;

            tStart = tic;
            [increMatching{tbimgm_qmIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_qmIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_qmIdx} = increMatching{tbimgm_qmIdx};
            
            acc{tbimgm_qmIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_qmIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{tbimgm_qmIdx} = cal_pair_graph_score(increMatching{tbimgm_qmIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{tbimgm_qmIdx} = cal_pair_graph_consistency(increMatching{tbimgm_qmIdx},nodeCnt,param.N+batchSize,0);
            
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
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            % previous matching: prevMatching{tbimgm_matchALSIdx}
            if parak == 1
                prevMatching{tbimgm_matchALSIdx} = baseMat;
            end
            matTmp{tbimgm_matchALSIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{tbimgm_matchALSIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{tbimgm_matchALSIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_matchALSIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_matchALSIdx},nodeCnt,param.N+batchSize,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;

            tStart = tic;
            [increMatching{tbimgm_matchALSIdx}, numPairMatch] = TBIMGM(affinity, simAP, matTmp{tbimgm_matchALSIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_matchALSIdx} = increMatching{tbimgm_matchALSIdx};
            
            acc{tbimgm_matchALSIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_matchALSIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{tbimgm_matchALSIdx} = cal_pair_graph_score(increMatching{tbimgm_matchALSIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{tbimgm_matchALSIdx} = cal_pair_graph_consistency(increMatching{tbimgm_matchALSIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, tbimgm_matchALSIdx, testk) = mean(acc{tbimgm_matchALSIdx}(:));
            scrAve(parak, tbimgm_matchALSIdx, testk) = mean(scr{tbimgm_matchALSIdx}(:));
            conPairAve(parak, tbimgm_matchALSIdx, testk) = mean(con{tbimgm_matchALSIdx}(:));
            timAve(parak, tbimgm_matchALSIdx, testk) = tEnd;
            countPairAve(parak, tbimgm_matchALSIdx, testk) = numPairMatch;
        end    
        
        fprintf('test in round %d/%d, Start from %d graphs, %d graphs incremented\n',testk, testCnt, baseGraphCnt, parak*batchSize);
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
    viewCnt=target.config.graphMinCnt + parak - 1;
    for algk = 1:algCnt
        timAveFull(parak,algk) = mean(timAve(parak,algk,:));
        accAveFull(parak,algk) = mean(accAve(parak,algk,:));
        scrAveFull(parak,algk) = mean(scrAve(parak,algk,:));
        conPairAveFull(parak,algk) = mean(conPairAve(parak,algk,:));
        countPairAveFull(parak, algk) = mean(countPairAve(parak, algk, :));
    end
    fprintf(' %02d,  %02d ',viewCnt,testCnt);fprintf(fidPerf,' %02d,  %02d',viewCnt,testCnt);
    for algk=1:algCnt
        if algSet.algEnable(algk)==0,continue;end
        fprintf('| %.3f %.3f %.3f %.3f %4d',accAveFull(parak,algk),scrAveFull(parak,algk),conPairAveFull(parak,algk),timAveFull(parak,algk), countPairAveFull(parak, algk));
        fprintf(fidPerf,', %.3f, %.3f, %.3f, %.3f, %4d',accAveFull(parak,algk),scrAveFull(parak,algk),conPairAveFull(parak,algk),timAveFull(parak,algk), countPairAveFull(parak, algk));% fprintf(cc,  score, consis
    end
    fprintf('\n');fprintf(fidPerf,'\n');
end

legendOff = 0;
savePath = sprintf('exp_online_%s_%s_batch_%d.mat', target.config.database, target.config.category, target.config.batchSize);
save(savePath, 'target', 'algSet', 'accAveFull', 'scrAveFull', 'conPairAveFull', 'timAveFull', 'countPairAveFull');
ave.accuracy = accAveFull;
ave.score = scrAveFull;
ave.consistency = conPairAveFull;
ave.time = timAveFull;
ave.matchingNumber = countPairAveFull;
fields = fieldnames(ave);
for ifield = 1:length(fields)
    xtag='Arriving graph';ytag=[fields{ifield}, '(',target.config.category,')'];
    plotResult_new(legendOff,target.config.graphRange-baseGraphCnt+1, getfield(ave,fields{ifield}), algSet, xtag, ytag);
end
