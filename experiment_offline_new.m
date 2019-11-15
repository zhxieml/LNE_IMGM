%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;

% target.config.graphMinCnt=19; 
% target.config.graphMaxCnt=50; 
target.config.testCnt = 20;% v
target.config.maxNumSearch = 18;
target.config.batchSize = 4;

target.config.graphMinCnt=32; 
target.config.graphMaxCnt=32; 
graphMaxCnt = target.config.graphMaxCnt;
graphMinCnt = target.config.graphMinCnt;
batchSize = target.config.batchSize;

graphCntRange = graphMinCnt:batchSize:graphMaxCnt;
% graphStep = 1;
target.config.database = "synthetic"; % "willow", "synthetic", "CMU-sequence"
load_target_data;

i = 0;
% savePath = sprintf('./exp/exp_offline_%s_%s_%s_%donly_full.mat', target.config.database, target.config.class, target.config.category, target.config.graphMinCnt);
savePath = sprintf('./exp/exp_offline_%s_%s_M=%d.mat', target.config.database, target.config.category, target.config.maxNumSearch);
while isfile(savePath)
    i = i + 1;
    savePath = [savePath(1:end-4) sprintf('_%d', i) savePath(end-3:end)];
end

% set algorithms
algNameSepSpace = '                    ';
algSet.algNameSet = {'cao_c_raw','imgm_d','imgm_r','anc_imgm','anc_imgm_a16','anc_imgm_d16', 'anc_imgm_a32','anc_imgm_d32'};
algSet.algEnable =  [     1,        0,       0,         1,          0,              0,             1,              1];
algSet.algColor = { cao_c_rawClr,imgm_dClr,imgm_rClr,anc_imgmClr,anc_imgmaClr16,anc_imgmdClr16, anc_imgmaClr16,anc_imgmdClr16};
algSet.algLineStyle = {'-',        '--',      '--',      '-.',        '-.',            '-.',          '-.'            '-.'};
algSet.algMarker = {   '.',        '.',       '.',       '.',        '.',            '.',              '.',            '.'};


[~,cao_c_rawIdx] = ismember('cao_c_raw',algSet.algNameSet);
[~,imgm_dIdx] = ismember('imgm_d',algSet.algNameSet);
[~,imgm_rIdx] = ismember('imgm_r',algSet.algNameSet);
[~,anc_imgmIdx] = ismember('anc_imgm',algSet.algNameSet);
[~,anc_imgm_a16Idx] = ismember('anc_imgm_a16',algSet.algNameSet);
[~,anc_imgm_d16Idx] = ismember('anc_imgm_d16',algSet.algNameSet);
[~,anc_imgm_a32Idx] = ismember('anc_imgm_a32',algSet.algNameSet);
[~,anc_imgm_d32Idx] = ismember('anc_imgm_d32',algSet.algNameSet);

% baseGraphCnt = target.config.graphMinCnt;
% target.config.graphRange = baseGraphCnt:graphStep:target.config.graphMaxCnt-graphStep;
% target.config.baseGraphCnt = baseGraphCnt;
nInlier = target.config.nInlier;
nOutlier = target.config.nOutlier;
target.config.nodeCnt = nInlier + nOutlier;
nodeCnt = target.config.nodeCnt;


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
% target.config.selectGraphMask{1} = 1:target.config.graphMaxCnt;
target.config.connect = 'nfc';

% data for experiment
paraCnt = length(graphCntRange);% paraCnt: iterate over graph #
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
matTmp = cell(1, algCnt);
Matching = cell(1, algCnt);

fidPerf = fopen('results.csv','w');
% fprintf(fidPerf, 'testType,database,category,testCnt,total node#,outlier#,graph#,alg#,scale,initConstWeight,consStep,constWeightMax,iterImmune\n');
% fprintf(fidPerf, '%s,      %s,      %s,      %d,     %d,         %d,      %d,    %d,  %f,   %f,             %f,      %.2f,          %d\n',...
%     target.config.testType, target.config.database, target.config.category,...
%     testCnt,nodeCnt,target.config.nOutlier,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,...
%     target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
% fprintf('testType=%s, database=%s, category=%s, test#=%d, node#=%d, outlier#=%d,graph#=%d, alg#=%d, scale=%.2f, initW=%.2f, stepW=%.2f, maxW=%.2f,iterImmune=%d\n',...
%     target.config.testType,target.config.database,target.config.category,testCnt,...
%     nodeCnt,target.config.nOutlier,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,...
%     target.config.initConstWeight, target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
% fprintf('\n');fprintf(fidPerf,'\n');

for parak = 1:paraCnt
    graphCnt = graphCntRange(parak);
    target.config.graphCnt = graphCnt;
    
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
        scrDenomCurrent = max(max(scrDenomMatInCnt));
        
        %%%%%%%%%%%% calculate the incremental matching with cao_c_raw %%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(cao_c_rawIdx)
            tStart = tic;
            Matching{cao_c_rawIdx} = CAO(rawMat, nodeCnt, graphCnt, target.config.iterRange, scrDenomCurrent, 'exact', 1);
            tEnd = toc(tStart);
            
            acc{cao_c_rawIdx} = cal_pair_graph_accuracy(Matching{cao_c_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
            scr{cao_c_rawIdx} = cal_pair_graph_score(Matching{cao_c_rawIdx},affinity.GT,nodeCnt,graphCnt);
            con{cao_c_rawIdx} = cal_pair_graph_consistency(Matching{cao_c_rawIdx},nodeCnt,graphCnt,0);
            accAve(parak, cao_c_rawIdx, testk) = mean(acc{cao_c_rawIdx}(:));
            scrAve(parak, cao_c_rawIdx, testk) = mean(scr{cao_c_rawIdx}(:));
            conPairAve(parak, cao_c_rawIdx, testk) = mean(con{cao_c_rawIdx}(:));
            timAve(parak, cao_c_rawIdx, testk) = tEnd;
            countPairAve(parak, cao_c_rawIdx, testk) = graphCnt;
        end
        
        %%%%%%%%%%% calculate the incremental matching with imgm_d %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_dIdx)
            %             param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            
            tStart = tic;
            [Matching{imgm_dIdx}, numPairMatch] = IMGM_offline(rawMat,nodeCnt,graphCnt,param);
            tEnd = toc(tStart);
            
            acc{imgm_dIdx} = cal_pair_graph_accuracy(Matching{imgm_dIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
            scr{imgm_dIdx} = cal_pair_graph_score(Matching{imgm_dIdx},affinity.GT,nodeCnt,graphCnt);
            con{imgm_dIdx} = cal_pair_graph_consistency(Matching{imgm_dIdx},nodeCnt,graphCnt,0);
            accAve(parak, imgm_dIdx, testk) = mean(acc{imgm_dIdx}(:));
            scrAve(parak, imgm_dIdx, testk) = mean(scr{imgm_dIdx}(:));
            conPairAve(parak, imgm_dIdx, testk) = mean(con{imgm_dIdx}(:));
            timAve(parak, imgm_dIdx, testk) = tEnd;
            countPairAve(parak, imgm_dIdx, testk) = numPairMatch;
        end
        
        %%%%%%%%%%% calculate the incremental matching with imgm_r %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_rIdx)
            %             simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            %             param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 3; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            
            tStart = tic;
            [Matching{imgm_rIdx}, numPairMatch] = IMGM_offline(rawMat,nodeCnt,graphCnt,param);
            tEnd = toc(tStart);
            
            acc{imgm_rIdx} = cal_pair_graph_accuracy(Matching{imgm_rIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
            scr{imgm_rIdx} = cal_pair_graph_score(Matching{imgm_rIdx},affinity.GT,nodeCnt,graphCnt);
            con{imgm_rIdx} = cal_pair_graph_consistency(Matching{imgm_rIdx},nodeCnt,graphCnt,0);
            accAve(parak, imgm_rIdx, testk) = mean(acc{imgm_rIdx}(:));
            scrAve(parak, imgm_rIdx, testk) = mean(scr{imgm_rIdx}(:));
            conPairAve(parak, imgm_rIdx, testk) = mean(con{imgm_rIdx}(:));
            timAve(parak, imgm_rIdx, testk) = tEnd;
            countPairAve(parak, imgm_rIdx, testk) = numPairMatch;
        end
        
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgmIdx)
            % param for tbimgm_cao_cIdx
            param.subMethodParam.name = 'CAO';
%             param.subMethodParam.useCstDecay = 1;
%             param.subMethodParam.cstDecay  = 0.7;
%             param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = target.config.iterRange;
%             param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            
            tStart = tic;
%             [P, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,batchSize,useAptOrder,opt,param)
            [Matching{anc_imgmIdx}, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,1,0,'None',param);
            tEnd = toc(tStart);
            
            acc{anc_imgmIdx} = cal_pair_graph_accuracy(Matching{anc_imgmIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
            scr{anc_imgmIdx} = cal_pair_graph_score(Matching{anc_imgmIdx},affinity.GT,nodeCnt,graphCnt);
            con{anc_imgmIdx} = cal_pair_graph_consistency(Matching{anc_imgmIdx},nodeCnt,graphCnt,0);
            accAve(parak, anc_imgmIdx, testk) = mean(acc{anc_imgmIdx}(:));
            scrAve(parak, anc_imgmIdx, testk) = mean(scr{anc_imgmIdx}(:));
            conPairAve(parak, anc_imgmIdx, testk) = mean(con{anc_imgmIdx}(:));
            timAve(parak, anc_imgmIdx, testk) = tEnd;
            countPairAve(parak, anc_imgmIdx, testk) = numPairMatch;
        end
        
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_a16 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_a16Idx)
            % param for tbimgm_cao_cIdx
            param.subMethodParam.name = 'CAO';
%             param.subMethodParam.useCstDecay = 1;
%             param.subMethodParam.cstDecay  = 0.7;
%             param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = target.config.iterRange;
%             param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            
            tStart = tic;
            [Matching{anc_imgm_a16Idx}, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,16,1,'a',param);
            tEnd = toc(tStart);
            
            acc{anc_imgm_a16Idx} = cal_pair_graph_accuracy(Matching{anc_imgm_a16Idx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
            scr{anc_imgm_a16Idx} = cal_pair_graph_score(Matching{anc_imgm_a16Idx},affinity.GT,nodeCnt,graphCnt);
            con{anc_imgm_a16Idx} = cal_pair_graph_consistency(Matching{anc_imgm_a16Idx},nodeCnt,graphCnt,0);
            accAve(parak, anc_imgm_a16Idx, testk) = mean(acc{anc_imgm_a16Idx}(:));
            scrAve(parak, anc_imgm_a16Idx, testk) = mean(scr{anc_imgm_a16Idx}(:));
            conPairAve(parak, anc_imgm_a16Idx, testk) = mean(con{anc_imgm_a16Idx}(:));
            timAve(parak, anc_imgm_a16Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_a16Idx, testk) = numPairMatch;
        end
        
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_d16 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_d16Idx)
            % param for tbimgm_cao_cIdx
            param.subMethodParam.name = 'CAO';
%             param.subMethodParam.useCstDecay = 1;
%             param.subMethodParam.cstDecay  = 0.7;
%             param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = target.config.iterRange;
%             param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            
            tStart = tic;
%             [P, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,batchSize,useAptOrder, param)
            [Matching{anc_imgm_d16Idx}, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,16,1,'d',param);
            tEnd = toc(tStart);
            
            acc{anc_imgm_d16Idx} = cal_pair_graph_accuracy(Matching{anc_imgm_d16Idx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
            scr{anc_imgm_d16Idx} = cal_pair_graph_score(Matching{anc_imgm_d16Idx},affinity.GT,nodeCnt,graphCnt);
            con{anc_imgm_d16Idx} = cal_pair_graph_consistency(Matching{anc_imgm_d16Idx},nodeCnt,graphCnt,0);
            accAve(parak, anc_imgm_d16Idx, testk) = mean(acc{anc_imgm_d16Idx}(:));
            scrAve(parak, anc_imgm_d16Idx, testk) = mean(scr{anc_imgm_d16Idx}(:));
            conPairAve(parak, anc_imgm_d16Idx, testk) = mean(con{anc_imgm_d16Idx}(:));
            timAve(parak, anc_imgm_d16Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_d16Idx, testk) = numPairMatch;
        end
        
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_a32 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_a32Idx)
            % param for tbimgm_cao_cIdx
            param.subMethodParam.name = 'CAO';
%             param.subMethodParam.useCstDecay = 1;
%             param.subMethodParam.cstDecay  = 0.7;
%             param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = target.config.iterRange;
%             param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            
            tStart = tic;
            [Matching{anc_imgm_a32Idx}, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,32,1,'a',param);
            tEnd = toc(tStart);
            
            acc{anc_imgm_a32Idx} = cal_pair_graph_accuracy(Matching{anc_imgm_a32Idx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
            scr{anc_imgm_a32Idx} = cal_pair_graph_score(Matching{anc_imgm_a32Idx},affinity.GT,nodeCnt,graphCnt);
            con{anc_imgm_a32Idx} = cal_pair_graph_consistency(Matching{anc_imgm_a32Idx},nodeCnt,graphCnt,0);
            accAve(parak, anc_imgm_a32Idx, testk) = mean(acc{anc_imgm_a32Idx}(:));
            scrAve(parak, anc_imgm_a32Idx, testk) = mean(scr{anc_imgm_a32Idx}(:));
            conPairAve(parak, anc_imgm_a32Idx, testk) = mean(con{anc_imgm_a32Idx}(:));
            timAve(parak, anc_imgm_a32Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_a32Idx, testk) = numPairMatch;
        end
        
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_d16 %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_d32Idx)
            % param for tbimgm_cao_cIdx
            param.subMethodParam.name = 'CAO';
%             param.subMethodParam.useCstDecay = 1;
%             param.subMethodParam.cstDecay  = 0.7;
%             param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = target.config.iterRange;
%             param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            
            tStart = tic;
%             [P, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,batchSize,useAptOrder, param)
            [Matching{anc_imgm_d32Idx}, numPairMatch] = ANC_IMGM_offline(rawMat,nodeCnt,graphCnt,32,1,'d',param);
            tEnd = toc(tStart);
            
            acc{anc_imgm_d32Idx} = cal_pair_graph_accuracy(Matching{anc_imgm_d32Idx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
            scr{anc_imgm_d32Idx} = cal_pair_graph_score(Matching{anc_imgm_d32Idx},affinity.GT,nodeCnt,graphCnt);
            con{anc_imgm_d32Idx} = cal_pair_graph_consistency(Matching{anc_imgm_d32Idx},nodeCnt,graphCnt,0);
            accAve(parak, anc_imgm_d32Idx, testk) = mean(acc{anc_imgm_d32Idx}(:));
            scrAve(parak, anc_imgm_d32Idx, testk) = mean(scr{anc_imgm_d32Idx}(:));
            conPairAve(parak, anc_imgm_d32Idx, testk) = mean(con{anc_imgm_d32Idx}(:));
            timAve(parak, anc_imgm_d32Idx, testk) = tEnd;
            countPairAve(parak, anc_imgm_d32Idx, testk) = numPairMatch;
        end
        
        fprintf('%d graphs, test in round %d/%d\n', graphCnt, testk, testCnt);
        fprintf('%-18s%-18s%-18s%-18s%-18s%-18s\n','field\alg', 'accuracy', 'score', 'consistency', 'time', 'numPairMatch');
        for alg = find(algSet.algEnable)
            fprintf('%-18s%-18f%-18f%-18f%-18f%-18d\n\n',algSet.algNameSet{alg}, accAve(parak, alg, testk), scrAve(parak, alg, testk), conPairAve(parak, alg, testk), timAve(parak, alg, testk), countPairAve(parak, alg, testk));
        end
    end
    
    for algk = 1:algCnt
        timAveFull(parak,algk) = mean(timAve(parak,algk,:));
        accAveFull(parak,algk) = mean(accAve(parak,algk,:));
        scrAveFull(parak,algk) = mean(scrAve(parak,algk,:));
        conPairAveFull(parak,algk) = mean(conPairAve(parak,algk,:));
        countPairAveFull(parak, algk) = mean(countPairAve(parak, algk, :));
    end
    
    fprintf('%d graphs overall\n', graphCnt);
    fprintf('%-18s%-18s%-18s%-18s%-18s%-18s\n','field\alg', 'accuracy', 'score', 'consistency', 'time', 'numPairMatch');
    for alg = find(algSet.algEnable)
        fprintf('%-18s%-18f%-18f%-18f%-18f%-18d\n\n',algSet.algNameSet{alg}, accAveFull(parak, alg), scrAveFull(parak, alg), conPairAveFull(parak, alg), timAveFull(parak, alg), countPairAveFull(parak, alg));
    end
end
    



    
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

fprintf(' %02d,  %02d ',graphCnt,testCnt);fprintf(fidPerf,' %02d,  %02d',graphCnt,testCnt);
for algk=1:algCnt
    if algSet.algEnable(algk)==0,continue;end
    fprintf('| %.3f %.3f %.3f %.3f %4d',accAveFull(algk),scrAveFull(algk),conPairAveFull(algk),timAveFull(algk), countPairAveFull(algk));
    fprintf(fidPerf,', %.3f, %.3f, %.3f, %.3f, %4d',accAveFull(algk),scrAveFull(algk),conPairAveFull(algk),timAveFull(algk), countPairAveFull(algk));% fprintf(cc,  score, consis
end
fprintf('\n');fprintf(fidPerf,'\n');

legendOff = 0;

save(savePath, 'target', 'algSet', 'accAveFull', 'scrAveFull', 'conPairAveFull', 'timAveFull', 'countPairAveFull');
ave.accuracy = accAveFull;
ave.score = scrAveFull;
ave.consistency = conPairAveFull;
ave.time = timAveFull;
ave.matchingNumber = countPairAveFull;
fields = fieldnames(ave);
for ifield = 1:length(fields)
    xtag='number of graph';ytag=[fields{ifield}, '(',target.config.category,')'];
    plotResult_new(legendOff,graphCntRange, getfield(ave,fields{ifield}), algSet, xtag, ytag);
end
