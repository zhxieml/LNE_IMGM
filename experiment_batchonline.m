%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;

target.config.graphMinCnt=0; 
target.config.graphMaxCnt=52; 
target.config.testCnt = 50; % v
target.config.maxNumSearch = 20;
target.config.batchSize = 1;
target.config.database = "synthetic"; % "willow", "synthetic"
load_target_data;

% set algorithms
algNameSepSpace = '                    ';
algSet.algNameSet = {'anc_imgm', 'anc_imgm_a4', 'anc_imgm_a8', 'anc_imgm_a16', 'anc_imgm_a32', 'anc_imgm_d4', 'anc_imgm_d8', 'anc_imgm_d16', 'anc_imgm_d32'};
algSet.algEnable =  [ 1,            1,             1,            1,             1,              1,              1,             1,              1      ];
algSet.algColor = { anc_imgmClr, anc_imgmaClr4, anc_imgmaClr8, anc_imgmaClr16, anc_imgmaClr32, anc_imgmdClr4, anc_imgmdClr8, anc_imgmdClr16, anc_imgmdClr32};
algSet.algLineStyle = {'-',       '--',           '--',            '--',           '--',         '-.',             '-.',         '-.',           '-.'};
algSet.algMarker =    {'.',       '.',            '.',             '.',            '.',          '.',              '.',          '.',            '.'};

[~,anc_imgmIdx] = ismember('anc_imgm',algSet.algNameSet);
[~,anc_imgm_a4Idx] = ismember('anc_imgm_a4',algSet.algNameSet);
[~,anc_imgm_a8Idx] = ismember('anc_imgm_a8',algSet.algNameSet);
[~,anc_imgm_a16Idx] = ismember('anc_imgm_a16',algSet.algNameSet);
[~,anc_imgm_a32Idx] = ismember('anc_imgm_a32',algSet.algNameSet);
[~,anc_imgm_d4Idx] = ismember('anc_imgm_d4',algSet.algNameSet);
[~,anc_imgm_d8Idx] = ismember('anc_imgm_d8', algSet.algNameSet);
[~,anc_imgm_d16Idx] = ismember('anc_imgm_d16', algSet.algNameSet);
[~,anc_imgm_d32Idx] = ismember('anc_imgm_d32', algSet.algNameSet);

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
    scrDenomCurrent = max(max(scrDenomMatInCnt(1:baseGraphCnt,1:baseGraphCnt)));
    baseMat = CAO(rawMat(1:nodeCnt*baseGraphCnt,1:nodeCnt*baseGraphCnt), nodeCnt, baseGraphCnt, target.config.iterRange, scrDenomCurrent, 'pair',1);
    
    for alg = find(algSet.algEnable)
        aptOrder{alg} = 1:graphCnt;
    end

    for parak = 1:paraCnt 
        param.n = nodeCnt; 
        param.N = baseGraphCnt + parak - 1; % 20
        param.subMethodParam.name = 'CAO';
        param.subMethodParam.iterMax = target.config.iterRange;
        param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N+1,1:param.N+1)));
        param.subMethodParam.optType = 'exact';
        param.subMethodParam.useCstInlier = 1;
        param.bVerbose = 0;
        param.maxNumSearch = target.config.maxNumSearch;


        %%%%%%%%%%%% calculate the incremental matching with anc_imgm %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgmIdx)
            if parak == 1
                prevMatching{anc_imgmIdx} = baseMat;
            end
            matTmp{anc_imgmIdx} = rawMat(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgmIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgmIdx};
            affScore = cal_pair_graph_inlier_score_local(affinity, matTmp{anc_imgmIdx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgmIdx}, numPairMatch] = ANC_IMGM(affinity, affScore, matTmp{anc_imgmIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgmIdx} = increMatching{anc_imgmIdx};
            
            acc{anc_imgmIdx} = cal_pair_graph_accuracy(increMatching{anc_imgmIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+1);
            scr{anc_imgmIdx} = cal_pair_graph_score(increMatching{anc_imgmIdx},affinity.GT,nodeCnt,param.N+1);
            % scr{anc_imgmIdx} = cal_pair_graph_inlier_score(increMatching{anc_imgmIdx},affinity.GT,nodeCnt,param.N+batchSize, target.config.inCnt);
            con{anc_imgmIdx} = cal_pair_graph_consistency(increMatching{anc_imgmIdx},nodeCnt,param.N+1,0);
            
            accAve(parak, anc_imgmIdx, testk) = mean(acc{anc_imgmIdx}(:));
            scrAve(parak, anc_imgmIdx, testk) = mean(scr{anc_imgmIdx}(:));
            conPairAve(parak, anc_imgmIdx, testk) = mean(con{anc_imgmIdx}(:));
            timAve(parak, anc_imgmIdx, testk) = tEnd;
            countPairAve(parak, anc_imgmIdx, testk) = numPairMatch;
        end

        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_a4 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_a4Idx)
            if parak == 1
                prevMatching{anc_imgm_a4Idx} = baseMat;
            end
            batchSize = 4;
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, 'a');
                subOrder = subOrder + param.N;
                aptOrder{anc_imgm_a4Idx}(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder{anc_imgm_a4Idx} = crop_rawMat(aptOrder{anc_imgm_a4Idx}, rawMat, nodeCnt);
                affinityReOrder{anc_imgm_a4Idx} = crop_affinity(aptOrder{anc_imgm_a4Idx}, affinity);
                targetReOrder{anc_imgm_a4Idx} = crop_target(aptOrder{anc_imgm_a4Idx}, target);
            end

            matTmp{anc_imgm_a4Idx} = rawMatReOrder{anc_imgm_a4Idx}(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgm_a4Idx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_a4Idx};
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder{anc_imgm_a4Idx}, matTmp{anc_imgm_a4Idx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgm_a4Idx}, numPairMatch] = ANC_IMGM(affinityReOrder{anc_imgm_a4Idx}, affScore, matTmp{anc_imgm_a4Idx}, targetReOrder{anc_imgm_a4Idx}, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgm_a4Idx} = increMatching{anc_imgm_a4Idx};
            
            acc{anc_imgm_a4Idx} = cal_pair_graph_accuracy(increMatching{anc_imgm_a4Idx},affinityReOrder{anc_imgm_a4Idx}.GT,nOutlier,nodeCnt,param.N+1);
            scr{anc_imgm_a4Idx} = cal_pair_graph_score_local(affinityReOrder{anc_imgm_a4Idx},increMatching{anc_imgm_a4Idx},nodeCnt,param.N+1);
            con{anc_imgm_a4Idx} = cal_pair_graph_consistency(increMatching{anc_imgm_a4Idx},nodeCnt,param.N+1,0);

            accAve(parak, anc_imgm_a4Idx, testk) = mean(acc{anc_imgm_a4Idx}(:));
            scrAve(parak, anc_imgm_a4Idx, testk) = mean(scr{anc_imgm_a4Idx}(:));
            conPairAve(parak, anc_imgm_a4Idx, testk) = mean(con{anc_imgm_a4Idx}(:));
            timAve(parak, anc_imgm_a4Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_a4Idx, testk) = numPairMatch;
        end
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_a8 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_a8Idx)
            if parak == 1
                prevMatching{anc_imgm_a8Idx} = baseMat;
            end
            batchSize = 8;
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, 'a');
                subOrder = subOrder + param.N;
                aptOrder{anc_imgm_a8Idx}(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder{anc_imgm_a8Idx} = crop_rawMat(aptOrder{anc_imgm_a8Idx}, rawMat, nodeCnt);
                affinityReOrder{anc_imgm_a8Idx} = crop_affinity(aptOrder{anc_imgm_a8Idx}, affinity);
                targetReOrder{anc_imgm_a8Idx} = crop_target(aptOrder{anc_imgm_a8Idx}, target);
            end

            matTmp{anc_imgm_a8Idx} = rawMatReOrder{anc_imgm_a8Idx}(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgm_a8Idx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_a8Idx};
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder{anc_imgm_a8Idx}, matTmp{anc_imgm_a8Idx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgm_a8Idx}, numPairMatch] = ANC_IMGM(affinityReOrder{anc_imgm_a8Idx}, affScore, matTmp{anc_imgm_a8Idx}, targetReOrder{anc_imgm_a8Idx}, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgm_a8Idx} = increMatching{anc_imgm_a8Idx};
            
            acc{anc_imgm_a8Idx} = cal_pair_graph_accuracy(increMatching{anc_imgm_a8Idx},affinityReOrder{anc_imgm_a8Idx}.GT,nOutlier,nodeCnt,param.N+1);
            scr{anc_imgm_a8Idx} = cal_pair_graph_score_local(affinityReOrder{anc_imgm_a8Idx},increMatching{anc_imgm_a8Idx},nodeCnt,param.N+1);
            con{anc_imgm_a8Idx} = cal_pair_graph_consistency(increMatching{anc_imgm_a8Idx},nodeCnt,param.N+1,0);

            accAve(parak, anc_imgm_a8Idx, testk) = mean(acc{anc_imgm_a8Idx}(:));
            scrAve(parak, anc_imgm_a8Idx, testk) = mean(scr{anc_imgm_a8Idx}(:));
            conPairAve(parak, anc_imgm_a8Idx, testk) = mean(con{anc_imgm_a8Idx}(:));
            timAve(parak, anc_imgm_a8Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_a8Idx, testk) = numPairMatch;
        end
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_a16 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_a16Idx)
            if parak == 1
                prevMatching{anc_imgm_a16Idx} = baseMat;
            end
            batchSize = 16;
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, 'a');
                subOrder = subOrder + param.N;
                aptOrder{anc_imgm_a16Idx}(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder{anc_imgm_a16Idx} = crop_rawMat(aptOrder{anc_imgm_a16Idx}, rawMat, nodeCnt);
                affinityReOrder{anc_imgm_a16Idx} = crop_affinity(aptOrder{anc_imgm_a16Idx}, affinity);
                targetReOrder{anc_imgm_a16Idx} = crop_target(aptOrder{anc_imgm_a16Idx}, target);
            end

            matTmp{anc_imgm_a16Idx} = rawMatReOrder{anc_imgm_a16Idx}(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgm_a16Idx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_a16Idx};
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder{anc_imgm_a16Idx}, matTmp{anc_imgm_a16Idx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgm_a16Idx}, numPairMatch] = ANC_IMGM(affinityReOrder{anc_imgm_a16Idx}, affScore, matTmp{anc_imgm_a16Idx}, targetReOrder{anc_imgm_a16Idx}, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgm_a16Idx} = increMatching{anc_imgm_a16Idx};
            
            acc{anc_imgm_a16Idx} = cal_pair_graph_accuracy(increMatching{anc_imgm_a16Idx},affinityReOrder{anc_imgm_a16Idx}.GT,nOutlier,nodeCnt,param.N+1);
            scr{anc_imgm_a16Idx} = cal_pair_graph_score_local(affinityReOrder{anc_imgm_a16Idx},increMatching{anc_imgm_a16Idx},nodeCnt,param.N+1);
            con{anc_imgm_a16Idx} = cal_pair_graph_consistency(increMatching{anc_imgm_a16Idx},nodeCnt,param.N+1,0);

            accAve(parak, anc_imgm_a16Idx, testk) = mean(acc{anc_imgm_a16Idx}(:));
            scrAve(parak, anc_imgm_a16Idx, testk) = mean(scr{anc_imgm_a16Idx}(:));
            conPairAve(parak, anc_imgm_a16Idx, testk) = mean(con{anc_imgm_a16Idx}(:));
            timAve(parak, anc_imgm_a16Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_a16Idx, testk) = numPairMatch;
        end
        
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_a32 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_a32Idx)
            if parak == 1
                prevMatching{anc_imgm_a32Idx} = baseMat;
            end
            batchSize = 32;
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, 'a');
                subOrder = subOrder + param.N;
                aptOrder{anc_imgm_a32Idx}(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder{anc_imgm_a32Idx} = crop_rawMat(aptOrder{anc_imgm_a32Idx}, rawMat, nodeCnt);
                affinityReOrder{anc_imgm_a32Idx} = crop_affinity(aptOrder{anc_imgm_a32Idx}, affinity);
                targetReOrder{anc_imgm_a32Idx} = crop_target(aptOrder{anc_imgm_a32Idx}, target);
            end

            matTmp{anc_imgm_a32Idx} = rawMatReOrder{anc_imgm_a32Idx}(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgm_a32Idx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_a32Idx};
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder{anc_imgm_a32Idx}, matTmp{anc_imgm_a32Idx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgm_a32Idx}, numPairMatch] = ANC_IMGM(affinityReOrder{anc_imgm_a32Idx}, affScore, matTmp{anc_imgm_a32Idx}, targetReOrder{anc_imgm_a32Idx}, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgm_a32Idx} = increMatching{anc_imgm_a32Idx};
            
            acc{anc_imgm_a32Idx} = cal_pair_graph_accuracy(increMatching{anc_imgm_a32Idx},affinityReOrder{anc_imgm_a32Idx}.GT,nOutlier,nodeCnt,param.N+1);
            scr{anc_imgm_a32Idx} = cal_pair_graph_score_local(affinityReOrder{anc_imgm_a32Idx},increMatching{anc_imgm_a32Idx},nodeCnt,param.N+1);
            con{anc_imgm_a32Idx} = cal_pair_graph_consistency(increMatching{anc_imgm_a32Idx},nodeCnt,param.N+1,0);

            accAve(parak, anc_imgm_a32Idx, testk) = mean(acc{anc_imgm_a32Idx}(:));
            scrAve(parak, anc_imgm_a32Idx, testk) = mean(scr{anc_imgm_a32Idx}(:));
            conPairAve(parak, anc_imgm_a32Idx, testk) = mean(con{anc_imgm_a32Idx}(:));
            timAve(parak, anc_imgm_a32Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_a32Idx, testk) = numPairMatch;
        end
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_d4 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_d4Idx)
            if parak == 1
                prevMatching{anc_imgm_d4Idx} = baseMat;
            end
            batchSize = 4;
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, 'd');
                subOrder = subOrder + param.N;
                aptOrder{anc_imgm_d4Idx}(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder{anc_imgm_d4Idx} = crop_rawMat(aptOrder{anc_imgm_d4Idx}, rawMat, nodeCnt);
                affinityReOrder{anc_imgm_d4Idx} = crop_affinity(aptOrder{anc_imgm_d4Idx}, affinity);
                targetReOrder{anc_imgm_d4Idx} = crop_target(aptOrder{anc_imgm_d4Idx}, target);
            end

            matTmp{anc_imgm_d4Idx} = rawMatReOrder{anc_imgm_d4Idx}(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgm_d4Idx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_d4Idx};
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder{anc_imgm_d4Idx}, matTmp{anc_imgm_d4Idx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgm_d4Idx}, numPairMatch] = ANC_IMGM(affinityReOrder{anc_imgm_d4Idx}, affScore, matTmp{anc_imgm_d4Idx}, targetReOrder{anc_imgm_d4Idx}, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgm_d4Idx} = increMatching{anc_imgm_d4Idx};
            
            acc{anc_imgm_d4Idx} = cal_pair_graph_accuracy(increMatching{anc_imgm_d4Idx},affinityReOrder{anc_imgm_d4Idx}.GT,nOutlier,nodeCnt,param.N+1);
            scr{anc_imgm_d4Idx} = cal_pair_graph_score_local(affinityReOrder{anc_imgm_d4Idx},increMatching{anc_imgm_d4Idx},nodeCnt,param.N+1);
            con{anc_imgm_d4Idx} = cal_pair_graph_consistency(increMatching{anc_imgm_d4Idx},nodeCnt,param.N+1,0);

            accAve(parak, anc_imgm_d4Idx, testk) = mean(acc{anc_imgm_d4Idx}(:));
            scrAve(parak, anc_imgm_d4Idx, testk) = mean(scr{anc_imgm_d4Idx}(:));
            conPairAve(parak, anc_imgm_d4Idx, testk) = mean(con{anc_imgm_d4Idx}(:));
            timAve(parak, anc_imgm_d4Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_d4Idx, testk) = numPairMatch;
        end
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_d8 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_d8Idx)
            if parak == 1
                prevMatching{anc_imgm_d8Idx} = baseMat;
            end
            batchSize = 8;
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, 'd');
                subOrder = subOrder + param.N;
                aptOrder{anc_imgm_d8Idx}(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder{anc_imgm_d8Idx} = crop_rawMat(aptOrder{anc_imgm_d8Idx}, rawMat, nodeCnt);
                affinityReOrder{anc_imgm_d8Idx} = crop_affinity(aptOrder{anc_imgm_d8Idx}, affinity);
                targetReOrder{anc_imgm_d8Idx} = crop_target(aptOrder{anc_imgm_d8Idx}, target);
            end

            matTmp{anc_imgm_d8Idx} = rawMatReOrder{anc_imgm_d8Idx}(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgm_d8Idx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_d8Idx};
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder{anc_imgm_d8Idx}, matTmp{anc_imgm_d8Idx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgm_d8Idx}, numPairMatch] = ANC_IMGM(affinityReOrder{anc_imgm_d8Idx}, affScore, matTmp{anc_imgm_d8Idx}, targetReOrder{anc_imgm_d8Idx}, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgm_d8Idx} = increMatching{anc_imgm_d8Idx};
            
            acc{anc_imgm_d8Idx} = cal_pair_graph_accuracy(increMatching{anc_imgm_d8Idx},affinityReOrder{anc_imgm_d8Idx}.GT,nOutlier,nodeCnt,param.N+1);
            scr{anc_imgm_d8Idx} = cal_pair_graph_score_local(affinityReOrder{anc_imgm_d8Idx},increMatching{anc_imgm_d8Idx},nodeCnt,param.N+1);
            con{anc_imgm_d8Idx} = cal_pair_graph_consistency(increMatching{anc_imgm_d8Idx},nodeCnt,param.N+1,0);

            accAve(parak, anc_imgm_d8Idx, testk) = mean(acc{anc_imgm_d8Idx}(:));
            scrAve(parak, anc_imgm_d8Idx, testk) = mean(scr{anc_imgm_d8Idx}(:));
            conPairAve(parak, anc_imgm_d8Idx, testk) = mean(con{anc_imgm_d8Idx}(:));
            timAve(parak, anc_imgm_d8Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_d8Idx, testk) = numPairMatch;
        end
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_d16 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_d16Idx)
            if parak == 1
                prevMatching{anc_imgm_d16Idx} = baseMat;
            end
            batchSize = 16;
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, 'd');
                subOrder = subOrder + param.N;
                aptOrder{anc_imgm_d16Idx}(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder{anc_imgm_d16Idx} = crop_rawMat(aptOrder{anc_imgm_d16Idx}, rawMat, nodeCnt);
                affinityReOrder{anc_imgm_d16Idx} = crop_affinity(aptOrder{anc_imgm_d16Idx}, affinity);
                targetReOrder{anc_imgm_d16Idx} = crop_target(aptOrder{anc_imgm_d16Idx}, target);
            end

            matTmp{anc_imgm_d16Idx} = rawMatReOrder{anc_imgm_d16Idx}(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgm_d16Idx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_d16Idx};
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder{anc_imgm_d16Idx}, matTmp{anc_imgm_d16Idx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgm_d16Idx}, numPairMatch] = ANC_IMGM(affinityReOrder{anc_imgm_d16Idx}, affScore, matTmp{anc_imgm_d16Idx}, targetReOrder{anc_imgm_d16Idx}, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgm_d16Idx} = increMatching{anc_imgm_d16Idx};
            
            acc{anc_imgm_d16Idx} = cal_pair_graph_accuracy(increMatching{anc_imgm_d16Idx},affinityReOrder{anc_imgm_d16Idx}.GT,nOutlier,nodeCnt,param.N+1);
            scr{anc_imgm_d16Idx} = cal_pair_graph_score_local(affinityReOrder{anc_imgm_d16Idx},increMatching{anc_imgm_d16Idx},nodeCnt,param.N+1);
            con{anc_imgm_d16Idx} = cal_pair_graph_consistency(increMatching{anc_imgm_d16Idx},nodeCnt,param.N+1,0);

            accAve(parak, anc_imgm_d16Idx, testk) = mean(acc{anc_imgm_d16Idx}(:));
            scrAve(parak, anc_imgm_d16Idx, testk) = mean(scr{anc_imgm_d16Idx}(:));
            conPairAve(parak, anc_imgm_d16Idx, testk) = mean(con{anc_imgm_d16Idx}(:));
            timAve(parak, anc_imgm_d16Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_d16Idx, testk) = numPairMatch;
        end
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_d32 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_d32Idx)
            if parak == 1
                prevMatching{anc_imgm_d32Idx} = baseMat;
            end
            batchSize = 32;
            if mod(param.N - baseGraphCnt, batchSize) == 0
                if graphCnt < param.N + batchSize
                    upCnt = graphCnt - param.N;
                else
                    upCnt = batchSize;
                end
                dimN = param.N*nodeCnt+1 : (param.N+upCnt)*nodeCnt;
                rawAffinity = crop_affinity(param.N+1:param.N+upCnt, affinity);
                subOrder = cal_adaptive_graph_order(rawAffinity, rawMat(dimN, dimN), nodeCnt, upCnt, 'd');
                subOrder = subOrder + param.N;
                aptOrder{anc_imgm_d32Idx}(param.N+1:param.N+upCnt) = subOrder;
                rawMatReOrder{anc_imgm_d32Idx} = crop_rawMat(aptOrder{anc_imgm_d32Idx}, rawMat, nodeCnt);
                affinityReOrder{anc_imgm_d32Idx} = crop_affinity(aptOrder{anc_imgm_d32Idx}, affinity);
                targetReOrder{anc_imgm_d32Idx} = crop_target(aptOrder{anc_imgm_d32Idx}, target);
            end

            matTmp{anc_imgm_d32Idx} = rawMatReOrder{anc_imgm_d32Idx}(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{anc_imgm_d32Idx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_d32Idx};
            affScore = cal_pair_graph_inlier_score_local(affinityReOrder{anc_imgm_d32Idx}, matTmp{anc_imgm_d32Idx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{anc_imgm_d32Idx}, numPairMatch] = ANC_IMGM(affinityReOrder{anc_imgm_d32Idx}, affScore, matTmp{anc_imgm_d32Idx}, targetReOrder{anc_imgm_d32Idx}, param);
            tEnd = toc(tStart);
            prevMatching{anc_imgm_d32Idx} = increMatching{anc_imgm_d32Idx};
            
            acc{anc_imgm_d32Idx} = cal_pair_graph_accuracy(increMatching{anc_imgm_d32Idx},affinityReOrder{anc_imgm_d32Idx}.GT,nOutlier,nodeCnt,param.N+1);
            scr{anc_imgm_d32Idx} = cal_pair_graph_score_local(affinityReOrder{anc_imgm_d32Idx},increMatching{anc_imgm_d32Idx},nodeCnt,param.N+1);
            con{anc_imgm_d32Idx} = cal_pair_graph_consistency(increMatching{anc_imgm_d32Idx},nodeCnt,param.N+1,0);

            accAve(parak, anc_imgm_d32Idx, testk) = mean(acc{anc_imgm_d32Idx}(:));
            scrAve(parak, anc_imgm_d32Idx, testk) = mean(scr{anc_imgm_d32Idx}(:));
            conPairAve(parak, anc_imgm_d32Idx, testk) = mean(con{anc_imgm_d32Idx}(:));
            timAve(parak, anc_imgm_d32Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_d32Idx, testk) = numPairMatch;
        end
        fprintf('test in round %d/%d, Start from %d graphs, %d graphs incremented\n',testk, testCnt, baseGraphCnt, parak);
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
savePath = sprintf('exp_batchonline_%s_%s.mat', target.config.database, target.config.category);
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
