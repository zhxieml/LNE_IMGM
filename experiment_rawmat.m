%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;

target.config.graphMinCnt=20; 
target.config.graphMaxCnt=50; 
target.config.testCnt = 1; % v
target.config.maxNumSearch = 20;
target.config.batchSize = 1;

target.config.database = "synthetic"; % "willow", "synthetic"
load_target_data;

if strcmp(target.config.database, 'synthetic')
    savePath = sprintf('./exp/exp_rawmat_%s_%s.mat', target.config.database, target.config.category);
elseif strcmp(target.config.database, 'willow')
    savePath = sprintf('./exp/exp_rawmat_%s_%s_%s.mat', target.config.database, target.config.class, target.config.category);
end

% set algorithms
algNameSepSpace = '                    ';
algSet.algNameSet = {'imgm_raw','imgm_cao','lne_imgm_raw','lne_imgm_cao'};
algSet.algEnable =  [ 1,            1,             1,           1       ];
algSet.algColor =   {imgm_rawClr, imgm_caoClr, lne_imgm_rawClr, lne_imgm_caoClr};
algSet.algLineStyle = {'--','-','--','-','-','--','-','--','-', '--', '-', '--', '-'};
algSet.algMarker = {'.','.','.','.','.','.', '.','.','.', '.', '.', '.', '.'};


[~,imgm_rawIdx] = ismember('imgm_raw',algSet.algNameSet);
[~,imgm_caoIdx] = ismember('imgm_cao',algSet.algNameSet);
[~,lne_imgm_rawIdx] = ismember('lne_imgm_raw', algSet.algNameSet);
[~,lne_imgm_caoIdx] = ismember('lne_imgm_cao', algSet.algNameSet);

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
    param.batchSize = batchSize;
    scrDenomCurrent = max(max(scrDenomMatInCnt(1:baseGraphCnt,1:baseGraphCnt)));
    baseMat_raw = rawMat(1:nodeCnt*baseGraphCnt,1:nodeCnt*baseGraphCnt);
    baseMat_cao = CAO(rawMat(1:nodeCnt*baseGraphCnt,1:nodeCnt*baseGraphCnt), nodeCnt, baseGraphCnt, target.config.iterRange,scrDenomCurrent, 'pair',1);
    param.n = nodeCnt; 

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
        
        %%%%%%%%%%% calculate the incremental matching with imgm_raw %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_rawIdx)
            if parak == 1
                prevMatching{imgm_rawIdx} = baseMat_raw;
            end
            matTmp{imgm_rawIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{imgm_rawIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{imgm_rawIdx};
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{imgm_rawIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{imgm_rawIdx},nodeCnt,param.N+batchSize,0);
            
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_rawIdx} = IMGM_local(affinity, simAP, matTmp{imgm_rawIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{imgm_rawIdx} = increMatching{imgm_rawIdx};
            
            acc{imgm_rawIdx} = cal_pair_graph_accuracy(increMatching{imgm_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{imgm_rawIdx} = cal_pair_graph_score(increMatching{imgm_rawIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{imgm_rawIdx} = cal_pair_graph_consistency(increMatching{imgm_rawIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, imgm_rawIdx, testk) = mean(acc{imgm_rawIdx}(:));
            scrAve(parak, imgm_rawIdx, testk) = mean(scr{imgm_rawIdx}(:));
            conPairAve(parak, imgm_rawIdx, testk) = mean(con{imgm_rawIdx}(:));
            timAve(parak, imgm_rawIdx, testk) = tEnd;
            countPairAve(parak, imgm_rawIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end

        %%%%%%%%%%% calculate the incremental matching with imgm_cao %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_caoIdx)
            if parak == 1
                prevMatching{imgm_caoIdx} = baseMat_cao;
            end
            matTmp{imgm_caoIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{imgm_caoIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{imgm_caoIdx};
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{imgm_caoIdx},affinity.GT(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)),nodeCnt,param.N+batchSize,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{imgm_caoIdx},nodeCnt,param.N+batchSize,0);
            
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_caoIdx} = IMGM_local(affinity, simAP, matTmp{imgm_caoIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{imgm_caoIdx} = increMatching{imgm_caoIdx};
            
            acc{imgm_caoIdx} = cal_pair_graph_accuracy(increMatching{imgm_caoIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{imgm_caoIdx} = cal_pair_graph_score(increMatching{imgm_caoIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{imgm_caoIdx} = cal_pair_graph_consistency(increMatching{imgm_caoIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, imgm_caoIdx, testk) = mean(acc{imgm_caoIdx}(:));
            scrAve(parak, imgm_caoIdx, testk) = mean(scr{imgm_caoIdx}(:));
            conPairAve(parak, imgm_caoIdx, testk) = mean(con{imgm_caoIdx}(:));
            timAve(parak, imgm_caoIdx, testk) = tEnd;
            countPairAve(parak, imgm_caoIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end

        %%%%%%%%%%%% calculate the incremental matching with lne_imgm_raw %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(lne_imgm_rawIdx)
            if parak == 1
                prevMatching{lne_imgm_rawIdx} = baseMat_raw;
            end
            matTmp{lne_imgm_rawIdx} = rawMat(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{lne_imgm_rawIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{lne_imgm_rawIdx};
            affScore = cal_pair_graph_inlier_score_local(affinity, matTmp{lne_imgm_rawIdx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{lne_imgm_rawIdx}, numPairMatch] = ANC_IMGM(affinity, affScore, matTmp{lne_imgm_rawIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{lne_imgm_rawIdx} = increMatching{lne_imgm_rawIdx};
            
            acc{lne_imgm_rawIdx} = cal_pair_graph_accuracy(increMatching{lne_imgm_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+1);
            scr{lne_imgm_rawIdx} = cal_pair_graph_score(increMatching{lne_imgm_rawIdx},affinity.GT,nodeCnt,param.N+1);
            con{lne_imgm_rawIdx} = cal_pair_graph_consistency(increMatching{lne_imgm_rawIdx},nodeCnt,param.N+1,0);
            
            accAve(parak, lne_imgm_rawIdx, testk) = mean(acc{lne_imgm_rawIdx}(:));
            scrAve(parak, lne_imgm_rawIdx, testk) = mean(scr{lne_imgm_rawIdx}(:));
            conPairAve(parak, lne_imgm_rawIdx, testk) = mean(con{lne_imgm_rawIdx}(:));
            timAve(parak, lne_imgm_rawIdx, testk) = tEnd;
            countPairAve(parak, lne_imgm_rawIdx, testk) = numPairMatch;
        end
        
        %%%%%%%%%%%% calculate the incremental matching with lne_imgm_cao %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(lne_imgm_caoIdx)
            if parak == 1
                prevMatching{lne_imgm_caoIdx} = baseMat_cao;
            end
            matTmp{lne_imgm_caoIdx} = rawMat(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{lne_imgm_caoIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{lne_imgm_caoIdx};
            affScore = cal_pair_graph_inlier_score_local(affinity, matTmp{lne_imgm_caoIdx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{lne_imgm_caoIdx}, numPairMatch] = ANC_IMGM(affinity, affScore, matTmp{lne_imgm_caoIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{lne_imgm_caoIdx} = increMatching{lne_imgm_caoIdx};
            
            acc{lne_imgm_caoIdx} = cal_pair_graph_accuracy(increMatching{lne_imgm_caoIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+1);
            scr{lne_imgm_caoIdx} = cal_pair_graph_score(increMatching{lne_imgm_caoIdx},affinity.GT,nodeCnt,param.N+1);
            con{lne_imgm_caoIdx} = cal_pair_graph_consistency(increMatching{lne_imgm_caoIdx},nodeCnt,param.N+1,0);
            
            accAve(parak, lne_imgm_caoIdx, testk) = mean(acc{lne_imgm_caoIdx}(:));
            scrAve(parak, lne_imgm_caoIdx, testk) = mean(scr{lne_imgm_caoIdx}(:));
            conPairAve(parak, lne_imgm_caoIdx, testk) = mean(con{lne_imgm_caoIdx}(:));
            timAve(parak, lne_imgm_caoIdx, testk) = tEnd;
            countPairAve(parak, lne_imgm_caoIdx, testk) = numPairMatch;
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
fprintf(algNamePreSpace);
for algk=1:algCnt
    if algSet.algEnable(algk)==0,continue;end
    fprintf([algSet.algNameSetDisplay{algk},algNameSepSpace]);
end
fprintf('\n');
fprintf('grh# itr#  ');
for algk=1:algCnt
    if algSet.algEnable(algk)==0,continue;end
    fprintf(' acc   scr   con   tim   pair');
end
fprintf('\n');
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
    fprintf(' %02d,  %02d ',viewCnt,testCnt);
    for algk=1:algCnt
        if algSet.algEnable(algk)==0,continue;end
        fprintf('| %.3f %.3f %.3f %.3f %4d',accAveFull(parak,algk),scrAveFull(parak,algk),conPairAveFull(parak,algk),timAveFull(parak,algk), countPairAveFull(parak, algk));
    end
    fprintf('\n');
end

save(savePath, 'target', 'algSet', 'accAveFull', 'scrAveFull', 'conPairAveFull', 'timAveFull', 'countPairAveFull');

legendOff = 0;
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
