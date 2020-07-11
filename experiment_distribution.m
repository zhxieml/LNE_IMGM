%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;

target.config.graphMinCnt=30; 
target.config.graphMaxCnt=50; 
target.config.testCnt = 1; % v
target.config.maxNumSearch = 20;
target.config.batchSize = 1;
target.config.database = "synthetic"; % "willow", "synthetic"
load_target_data;

if strcmp(target.config.database, 'synthetic')
    savePath = sprintf('./exp/exp_distribution_%s_%s.mat', target.config.database, target.config.category);
elseif strcmp(target.config.database, 'willow')
    savePath = sprintf('./exp/exp_distribution_%s_%s_%s.mat', target.config.database, target.config.class, target.config.category);
end

distribution = zeros(target.config.testCnt, target.config.graphMaxCnt - target.config.graphMinCnt, 10);

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
testCnt = target.config.testCnt;% testCnt: iterate over tests

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


        %%%%%%%%%%%% calculate the incremental matching with lne_imgm %%%%%%%%%%%%%%%%%%%%
        if parak == 1
            prevMatching = baseMat;
        end
        matTmp = rawMat(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
        matTmp(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching;
        affScore = cal_pair_graph_inlier_score_local(affinity, matTmp, nodeCnt, param.N+1, nodeCnt);
        tStart = tic;
        [increMatching, numPairMatch] = LNE_IMGM(affinity, affScore, matTmp, target, param);
        tEnd = toc(tStart);
        prevMatching = increMatching;

        acc = cal_pair_graph_accuracy(increMatching,affinity.GT,target.config.nOutlier,nodeCnt,param.N+1);

        h = histogram(acc, 10);
        distribution(testk, parak, :) = h.Values / sum(h.Values);
        disp(distribution(testk, parak, 10));

        scr = cal_pair_graph_score(increMatching,affinity.GT,nodeCnt,param.N+1);
        con = cal_pair_graph_consistency(increMatching,nodeCnt,param.N+1,0);
           
        fprintf('test in round %d/%d, Start from %d graphs, %d graphs incremented\n',testk, testCnt, baseGraphCnt, parak);
        fprintf('%-18s%-18s%-18s%-18s%-18s%-18s\n','field\alg', 'accuracy', 'score', 'consistency', 'time', 'numPairMatch');
        fprintf('%-18s%-18f%-18f%-18f%-18f%-18d\n\n','lne_imgm',mean(acc(:)),mean(scr(:)),mean(con(:)),tEnd,numPairMatch);
    end

end % for paraCnt

save(savePath, 'distribution');
