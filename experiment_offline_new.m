%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;

target.config.graphMinCnt=19; 
target.config.graphMaxCnt=50; 
target.config.testCnt = 10;% v
target.config.maxNumSearch = 20;
graphStep = 1;
target.config.database = "synthetic"; % "willow", "synthetic", "CMU-sequence"
load_target_data;

% set algorithms
algNameSepSpace = '                    ';
algSet.algNameSet = {'cao_pc_inc', 'cao_pc_raw', 'cao_c_inc','cao_c_raw','imgm_d','imgm_r','tbimgm_cao_c','tbimgm_cao_pc','tbimgm_cao_cst', 'tbimgm_qm', 'tbimgm_matchALS', 'tbimgm_cao_c_ada'};
algSet.algEnable =  [ 0,            0,             0,           1,          1,       1,       1,              0,              0,               0,             0,                     1];
algSet.algColor = { cao_pcClr,cao_pc_rawClr,cao_cClr,cao_c_rawClr,imgm_dClr,imgm_rClr, tbimgm_cao_cClr, tbimgm_cao_pcClr, tbimgm_cao_cstClr, tbimgm_qmClr, tbimgm_matchALSClr, tbimgm_cao_c_adaClr};
algSet.algLineStyle = {'--','--','-','--','-','--','-','--','-', '--', '-', '-'};
algSet.algMarker = {'.','.','.','.','.','.','.','.','.', '.', '.', '.'};

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
[~,tbimgm_cao_c_adaIdx] = ismember('tbimgm_cao_c_ada', algSet.algNameSet);



% baseGraphCnt = target.config.graphMinCnt;
% target.config.graphRange = baseGraphCnt:graphStep:target.config.graphMaxCnt-graphStep;
% target.config.baseGraphCnt = baseGraphCnt;
nInlier = target.config.nInlier;
nOutlier = target.config.nOutlier;
target.config.nodeCnt = nInlier + nOutlier;
target.config.graphCnt = 30;
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
% paraCnt = length(target.config.graphRange);% paraCnt: iterate over graph #
algCnt = length(algSet.algNameSet);% algCnt: iterate over algorithms
testCnt = target.config.testCnt;% testCnt: iterate over tests
timAve = zeros(algCnt,testCnt);timAveFull = zeros(algCnt);
accAve = zeros(algCnt,testCnt);accAveFull = zeros(algCnt);
scrAve = zeros(algCnt,testCnt);scrAveFull = zeros(algCnt);
conPairAve = zeros(algCnt,testCnt);conPairAveFull = zeros(algCnt);
countPairAve = zeros(algCnt,testCnt);countPairAveFull = zeros(algCnt);

acc = cell(1, algCnt);
scr = cell(1, algCnt);
con = cell(1, algCnt);
matTmp = cell(1, algCnt);
Matching = cell(1, algCnt);

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
    scrDenomCurrent = max(max(scrDenomMatInCnt));
        
    %%%%%%%%%%%% calculate the incremental matching with cao_c_raw %%%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(cao_c_rawIdx)
        tStart = tic;
        Matching{cao_c_rawIdx} = CAO(rawMat, nodeCnt, graphCnt, target.config.iterRange, scrDenomCurrent, 'exact', 1);
        tEnd = toc(tStart);
        
        acc{cao_c_rawIdx} = cal_pair_graph_accuracy(Matching{cao_c_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{cao_c_rawIdx} = cal_pair_graph_score(Matching{cao_c_rawIdx},affinity.GT,nodeCnt,graphCnt);
        con{cao_c_rawIdx} = cal_pair_graph_consistency(Matching{cao_c_rawIdx},nodeCnt,graphCnt,0);
        accAve(cao_c_rawIdx, testk) = mean(acc{cao_c_rawIdx}(:));
        scrAve(cao_c_rawIdx, testk) = mean(scr{cao_c_rawIdx}(:));
        conPairAve(cao_c_rawIdx, testk) = mean(con{cao_c_rawIdx}(:));
        timAve(cao_c_rawIdx, testk) = tEnd;
        countPairAve(cao_c_rawIdx, testk) = graphCnt;
    end
    
    %%%%%%%%%%% calculate the incremental matching with imgm_d %%%%%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(imgm_dIdx)
        %             param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
        param.iterMax = target.config.iterRange;
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
    
    %%%%%%%%%%% calculate the incremental matching with imgm_r %%%%%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(imgm_rIdx)
        %             simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
        %             param.subMethodParam.scrDenom = max(max(scrDenomMatInCntTmp(1:param.N,1:param.N)));
        param.iterMax = target.config.iterRange;
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
    
    %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_c %%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(tbimgm_cao_cIdx)
        % param for tbimgm_cao_cIdx
        param.subMethodParam.name = 'CAO';
        param.subMethodParam.useCstDecay = 1;
        param.subMethodParam.cstDecay  = 0.7;
        param.subMethodParam.useWeightedDecay  = 0;
        param.subMethodParam.iterMax = target.config.iterRange;
        param.subMethodParam.scrDenom = scrDenomCurrent;
        param.subMethodParam.optType = 'exact';
        param.subMethodParam.useCstInlier = 1;
        param.bVerbose = 0;
        param.maxNumSearch = target.config.maxNumSearch;
        
        % param for adaptive order
        param.adaptive_order = false;
        
        tStart = tic;
        Matching{tbimgm_cao_cIdx} = TBIMGM_offline(rawMat,nodeCnt,graphCnt,param);
        tEnd = toc(tStart);
        
        acc{tbimgm_cao_cIdx} = cal_pair_graph_accuracy(Matching{tbimgm_cao_cIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{tbimgm_cao_cIdx} = cal_pair_graph_score(Matching{tbimgm_cao_cIdx},affinity.GT,nodeCnt,graphCnt);
        con{tbimgm_cao_cIdx} = cal_pair_graph_consistency(Matching{tbimgm_cao_cIdx},nodeCnt,graphCnt,0);
        accAve(tbimgm_cao_cIdx, testk) = mean(acc{tbimgm_cao_cIdx}(:));
        scrAve(tbimgm_cao_cIdx, testk) = mean(scr{tbimgm_cao_cIdx}(:));
        conPairAve(tbimgm_cao_cIdx, testk) = mean(con{tbimgm_cao_cIdx}(:));
        timAve(tbimgm_cao_cIdx, testk) = tEnd;
        countPairAve(tbimgm_cao_cIdx, testk) = graphCnt;
    end
    
    %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_c_ada %%%%%%%%%%%%%%%%%%%%
    if algSet.algEnable(tbimgm_cao_c_adaIdx)
        % param for tbimgm_cao_cIdx
        param.subMethodParam.name = 'CAO';
        param.subMethodParam.useCstDecay = 1;
        param.subMethodParam.cstDecay  = 0.7;
        param.subMethodParam.useWeightedDecay  = 0;
        param.subMethodParam.iterMax = target.config.iterRange;
        param.subMethodParam.scrDenom = scrDenomCurrent;
        param.subMethodParam.optType = 'exact';
        param.subMethodParam.useCstInlier = 1;
        param.bVerbose = 0;
        param.maxNumSearch = target.config.maxNumSearch;
        
        % param for adaptive order
        param.adaptive_order = true;
        
        tStart = tic;
        Matching{tbimgm_cao_c_adaIdx} = TBIMGM_offline(rawMat,nodeCnt,graphCnt,param);
        tEnd = toc(tStart);
        
        acc{tbimgm_cao_c_adaIdx} = cal_pair_graph_accuracy(Matching{tbimgm_cao_c_adaIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
        scr{tbimgm_cao_c_adaIdx} = cal_pair_graph_score(Matching{tbimgm_cao_c_adaIdx},affinity.GT,nodeCnt,graphCnt);
        con{tbimgm_cao_c_adaIdx} = cal_pair_graph_consistency(Matching{tbimgm_cao_c_adaIdx},nodeCnt,graphCnt,0);
        accAve(tbimgm_cao_c_adaIdx, testk) = mean(acc{tbimgm_cao_c_adaIdx}(:));
        scrAve(tbimgm_cao_c_adaIdx, testk) = mean(scr{tbimgm_cao_c_adaIdx}(:));
        conPairAve(tbimgm_cao_c_adaIdx, testk) = mean(con{tbimgm_cao_c_adaIdx}(:));
        timAve(tbimgm_cao_c_adaIdx, testk) = tEnd;
        countPairAve(tbimgm_cao_c_adaIdx, testk) = graphCnt;
    end
    
    fprintf('test in round %d\n',testk);
    fprintf('%-18s%-18s%-18s%-18s%-18s%-18s\n','field\alg', 'accuracy', 'score', 'consistency', 'time', 'numPairMatch');
    for alg = find(algSet.algEnable)
        fprintf('%-18s%-18f%-18f%-18f%-18f%-18d\n\n',algSet.algNameSet{alg}, accAve(alg, testk), scrAve(alg, testk), conPairAve(alg, testk), timAve(alg, testk), countPairAve(alg, testk));
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
savePath = sprintf('exp_offline_%s_%s.mat', target.config.database, target.config.category);
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
