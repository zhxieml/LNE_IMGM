%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;

target.config.graphMinCnt=20; 
target.config.graphMaxCnt=50; 
target.config.testCnt = 20; % v
target.config.maxNumSearch = 20;
target.config.batchSize = 1;
target.config.database = "synthetic"; % "willow", "synthetic"
load_target_data;

if strcmp(target.config.database, 'synthetic')
    savePath = sprintf('./exp/exp_online_%s_%s.mat', target.config.database, target.config.category);
elseif strcmp(target.config.database, 'willow')
    savePath = sprintf('./exp/exp_online_%s_%s_%s.mat', target.config.database, target.config.class, target.config.category);
end

% set algorithms
algNameSepSpace = '                    ';
algSet.algNameSet = {'cao_c_inc','cao_c_raw','imgm_d','imgm_r','lne_imgm'};
algSet.algEnable =  [    0,           0,         1,       1,        1];
algSet.algColor = { cao_cClr,cao_c_rawClr,imgm_dClr,imgm_rClr,lne_imgmClr};
algSet.algLineStyle = {'-',       '--',      '-',      '--',      '-'};
algSet.algMarker =    {'.',       '.',       '.',      '.',       '.'};

[~,cao_cIdx] = ismember('cao_c_inc',algSet.algNameSet);
[~,cao_c_rawIdx] = ismember('cao_c_raw',algSet.algNameSet);
[~,imgm_dIdx] = ismember('imgm_d',algSet.algNameSet);
[~,imgm_rIdx] = ismember('imgm_r',algSet.algNameSet);
[~,lne_imgmIdx] = ismember('lne_imgm',algSet.algNameSet);

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
    sigma = 0.3;
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


        %%%%%%%%%%%% calculate the incremental matching with lne_imgm %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(lne_imgmIdx)
            if parak == 1
                prevMatching{lne_imgmIdx} = baseMat;
            end
            matTmp{lne_imgmIdx} = rawMat(1:nodeCnt*(param.N+1),1:nodeCnt*(param.N+1));
            matTmp{lne_imgmIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{lne_imgmIdx};
            affScore = cal_pair_graph_inlier_score_local(affinity, matTmp{lne_imgmIdx}, nodeCnt, param.N+1, nodeCnt);
            tStart = tic;
            [increMatching{lne_imgmIdx}, numPairMatch] = LNE_IMGM(affinity, affScore, matTmp{lne_imgmIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{lne_imgmIdx} = increMatching{lne_imgmIdx};
            
            acc{lne_imgmIdx} = cal_pair_graph_accuracy(increMatching{lne_imgmIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+1);
            scr{lne_imgmIdx} = cal_pair_graph_score(increMatching{lne_imgmIdx},affinity.GT,nodeCnt,param.N+1);
            con{lne_imgmIdx} = cal_pair_graph_consistency(increMatching{lne_imgmIdx},nodeCnt,param.N+1,0);
            
            accAve(parak, lne_imgmIdx, testk) = mean(acc{lne_imgmIdx}(:));
            scrAve(parak, lne_imgmIdx, testk) = mean(scr{lne_imgmIdx}(:));
            conPairAve(parak, lne_imgmIdx, testk) = mean(con{lne_imgmIdx}(:));
            timAve(parak, lne_imgmIdx, testk) = tEnd;
            countPairAve(parak, lne_imgmIdx, testk) = numPairMatch;
        end

        %%%%%%%%%%% calculate the incremental matching with cao_c_inc %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            con{cao_cIdx} = cal_pair_graph_consistency(increMatching{cao_cIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, cao_cIdx, testk) = mean(acc{cao_cIdx}(:));
            scrAve(parak, cao_cIdx, testk) = mean(scr{cao_cIdx}(:));
            conPairAve(parak, cao_cIdx, testk) = mean(con{cao_cIdx}(:));
            timAve(parak, cao_cIdx, testk) = tEnd;
            countPairAve(parak, cao_cIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
        end
        
        %%%%%%%%%%%% calculate the incremental matching with cao_c_raw %%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(cao_c_rawIdx)
            tStart = tic;
            increMatching{cao_c_rawIdx} = CAO(rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize)), nodeCnt, param.N+batchSize, target.config.iterRange, scrDenomCurrent, 'exact', 1);
            tEnd = toc(tStart);
            
            acc{cao_c_rawIdx} = cal_pair_graph_accuracy(increMatching{cao_c_rawIdx}, affinity.GT, target.config.nOutlier, nodeCnt, param.N+batchSize);
            scr{cao_c_rawIdx} = cal_pair_graph_score(increMatching{cao_c_rawIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{cao_c_rawIdx} = cal_pair_graph_consistency(increMatching{cao_c_rawIdx},nodeCnt,param.N+batchSize,0);
            accAve(parak, cao_c_rawIdx, testk) = mean(acc{cao_c_rawIdx}(:));
            scrAve(parak, cao_c_rawIdx, testk) = mean(scr{cao_c_rawIdx}(:));
            conPairAve(parak, cao_c_rawIdx, testk) = mean(con{cao_c_rawIdx}(:));
            timAve(parak, cao_c_rawIdx, testk) = tEnd;
            countPairAve(parak, cao_c_rawIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
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
            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_dIdx} = IMGM_local(affinity, simAP, matTmp{imgm_dIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{imgm_dIdx} = increMatching{imgm_dIdx};
            
            acc{imgm_dIdx} = cal_pair_graph_accuracy(increMatching{imgm_dIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{imgm_dIdx} = cal_pair_graph_score(increMatching{imgm_dIdx},affinity.GT,nodeCnt,param.N+batchSize);
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
            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 3; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_rIdx} = IMGM_local(affinity, simAP, matTmp{imgm_rIdx}, target, param);
            tEnd = toc(tStart);
            prevMatching{imgm_rIdx} = increMatching{imgm_rIdx};
            
            acc{imgm_rIdx} = cal_pair_graph_accuracy(increMatching{imgm_rIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{imgm_rIdx} = cal_pair_graph_score(increMatching{imgm_rIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{imgm_rIdx} = cal_pair_graph_consistency(increMatching{imgm_rIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, imgm_rIdx, testk) = mean(acc{imgm_rIdx}(:));
            scrAve(parak, imgm_rIdx, testk) = mean(scr{imgm_rIdx}(:));
            conPairAve(parak, imgm_rIdx, testk) = mean(con{imgm_rIdx}(:));
            timAve(parak, imgm_rIdx, testk) = tEnd;
            countPairAve(parak, imgm_rIdx, testk) = param.N*batchSize + (batchSize+1)*batchSize/2;
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
