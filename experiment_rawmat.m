%% Initialization
clear *;clear -global *;close all;clc;
global affinity target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;

target.config.graphMinCnt=20; 
target.config.graphMaxCnt=50; 
target.config.testCnt = 50; % v
target.config.maxNumSearch = 20;
target.config.batchSize = 1;

target.config.database = "willow"; % "willow", "synthetic"
load_target_data;

% set algorithms
algNameSepSpace = '                    ';
algSet.algNameSet = {'imgm_raw','imgm_cao','anc_imgm_raw','anc_imgm_cao'};
algSet.algEnable =  [ 1,            1,             1,           1       ];
algSet.algColor =   {imgm_rawClr, imgm_caoClr, anc_imgm_rawClr, anc_imgm_caoClr};
algSet.algLineStyle = {'--','-','--','-','-','--','-','--','-', '--', '-', '--', '-'};
algSet.algMarker = {'.','.','.','.','.','.', '.','.','.', '.', '.', '.', '.'};


[~,imgm_rawIdx] = ismember('imgm_raw',algSet.algNameSet);
[~,imgm_caoIdx] = ismember('imgm_cao',algSet.algNameSet);
[~,anc_imgm_rawIdx] = ismember('anc_imgm_raw', algSet.algNameSet);
[~,anc_imgm_caoIdx] = ismember('anc_imgm_cao', algSet.algNameSet);

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
    param.batchSize = batchSize;
    scrDenomCurrent = max(max(scrDenomMatInCnt(1:baseGraphCnt,1:baseGraphCnt)));
    baseMat_raw = rawMat(1:nodeCnt*baseGraphCnt,1:nodeCnt*baseGraphCnt);
    baseMat_cao = CAO(rawMat(1:nodeCnt*baseGraphCnt,1:nodeCnt*baseGraphCnt), nodeCnt, baseGraphCnt, target.config.iterRange,scrDenomCurrent, 'pair',1);
    param.n = nodeCnt; 

    for parak = 1:paraCnt     
        param.N = baseGraphCnt + (parak-1)*batchSize; % 20

        %%%%%%%%%%% calculate the incremental matching with imgm_raw %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_rawIdx)
            if parak == 1
                prevMatching{imgm_rawIdx} = baseMat_raw;
            end
            matTmp{imgm_rawIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{imgm_rawIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{imgm_rawIdx};

            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            [increMatching{imgm_rawIdx}, numPairMatch] = IMGM_batch(affinity, target, matTmp{imgm_rawIdx}, nodeCnt, param.N, batchSize, 0, param);
                                                     % IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, param)
            tEnd = toc(tStart);
            prevMatching{imgm_rawIdx} = increMatching{imgm_rawIdx};
            
            acc{imgm_rawIdx} = cal_pair_graph_accuracy(increMatching{imgm_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{imgm_rawIdx} = cal_pair_graph_score(increMatching{imgm_rawIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{imgm_rawIdx} = cal_pair_graph_consistency(increMatching{imgm_rawIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, imgm_rawIdx, testk) = mean(acc{imgm_rawIdx}(:));
            scrAve(parak, imgm_rawIdx, testk) = mean(scr{imgm_rawIdx}(:));
            conPairAve(parak, imgm_rawIdx, testk) = mean(con{imgm_rawIdx}(:));
            timAve(parak, imgm_rawIdx, testk) = tEnd;
            countPairAve(parak, imgm_rawIdx, testk) = numPairMatch;
        end

        %%%%%%%%%%% calculate the incremental matching with imgm_cao %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_caoIdx)
            if parak == 1
                prevMatching{imgm_caoIdx} = baseMat_cao;
            end
            matTmp{imgm_caoIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{imgm_caoIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{imgm_caoIdx};

            param.iterMax = target.config.iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            [increMatching{imgm_caoIdx}, numPairMatch] = IMGM_batch(affinity, target, matTmp{imgm_caoIdx}, nodeCnt, param.N, batchSize, 0, param);
                                                     % IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, param)
            tEnd = toc(tStart);
            prevMatching{imgm_caoIdx} = increMatching{imgm_caoIdx};
            
            acc{imgm_caoIdx} = cal_pair_graph_accuracy(increMatching{imgm_caoIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{imgm_caoIdx} = cal_pair_graph_score(increMatching{imgm_caoIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{imgm_caoIdx} = cal_pair_graph_consistency(increMatching{imgm_caoIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, imgm_caoIdx, testk) = mean(acc{imgm_caoIdx}(:));
            scrAve(parak, imgm_caoIdx, testk) = mean(scr{imgm_caoIdx}(:));
            conPairAve(parak, imgm_caoIdx, testk) = mean(con{imgm_caoIdx}(:));
            timAve(parak, imgm_caoIdx, testk) = tEnd;
            countPairAve(parak, imgm_caoIdx, testk) = numPairMatch;
        end

        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_raw %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_rawIdx)
            % param for anc_imgm_rawIdx
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.iterMax = target.config.iterRange;
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            % previous matching: prevMatching{anc_imgm_rawIdx}
            if parak == 1
                prevMatching{anc_imgm_rawIdx} = baseMat_raw;
            end
            
            matTmp{anc_imgm_rawIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{anc_imgm_rawIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_rawIdx};

            tStart = tic;
            [increMatching{anc_imgm_rawIdx}, numPairMatch] = ANC_IMGM_batch(affinity, target, matTmp{anc_imgm_rawIdx}, nodeCnt, param.N, batchSize, 0, param);
                                                            %ANC_IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, param)
            tEnd = toc(tStart);
            prevMatching{anc_imgm_rawIdx} = increMatching{anc_imgm_rawIdx};
            
            acc{anc_imgm_rawIdx} = cal_pair_graph_accuracy(increMatching{anc_imgm_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{anc_imgm_rawIdx} = cal_pair_graph_score(increMatching{anc_imgm_rawIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{anc_imgm_rawIdx} = cal_pair_graph_consistency(increMatching{anc_imgm_rawIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, anc_imgm_rawIdx, testk) = mean(acc{anc_imgm_rawIdx}(:));
            scrAve(parak, anc_imgm_rawIdx, testk) = mean(scr{anc_imgm_rawIdx}(:));
            conPairAve(parak, anc_imgm_rawIdx, testk) = mean(con{anc_imgm_rawIdx}(:));
            timAve(parak, anc_imgm_rawIdx, testk) = tEnd;
            countPairAve(parak, anc_imgm_rawIdx, testk) = numPairMatch;
        end
        %%%%%%%%%%%% calculate the incremental matching with anc_imgm_raw %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(anc_imgm_caoIdx)
            % param for anc_imgm_caoIdx
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.iterMax = target.config.iterRange;
            param.subMethodParam.optType = 'exact';
            param.subMethodParam.useCstInlier = 1;
            param.bVerbose = 0;
            param.maxNumSearch = target.config.maxNumSearch;
            % previous matching: prevMatching{anc_imgm_caoIdx}
            if parak == 1
                prevMatching{anc_imgm_caoIdx} = baseMat_cao;
            end
            matTmp{anc_imgm_caoIdx} = rawMat(1:nodeCnt*(param.N+batchSize),1:nodeCnt*(param.N+batchSize));
            matTmp{anc_imgm_caoIdx}(1:nodeCnt*param.N,1:nodeCnt*param.N)=prevMatching{anc_imgm_caoIdx};

            tStart = tic;
            [increMatching{anc_imgm_caoIdx}, numPairMatch] = ANC_IMGM_batch(affinity, target, matTmp{anc_imgm_caoIdx}, nodeCnt, param.N, batchSize, 0, param);
                                                            %ANC_IMGM_batch(affinity, target, rawMat, nodeCnt, baseGraphCnt, batchSize, useAptOrder, param)
            tEnd = toc(tStart);
            prevMatching{anc_imgm_caoIdx} = increMatching{anc_imgm_caoIdx};
            
            acc{anc_imgm_caoIdx} = cal_pair_graph_accuracy(increMatching{anc_imgm_caoIdx},affinity.GT,target.config.nOutlier,nodeCnt,param.N+batchSize);
            scr{anc_imgm_caoIdx} = cal_pair_graph_score(increMatching{anc_imgm_caoIdx},affinity.GT,nodeCnt,param.N+batchSize);
            con{anc_imgm_caoIdx} = cal_pair_graph_consistency(increMatching{anc_imgm_caoIdx},nodeCnt,param.N+batchSize,0);
            
            accAve(parak, anc_imgm_caoIdx, testk) = mean(acc{anc_imgm_caoIdx}(:));
            scrAve(parak, anc_imgm_caoIdx, testk) = mean(scr{anc_imgm_caoIdx}(:));
            conPairAve(parak, anc_imgm_caoIdx, testk) = mean(con{anc_imgm_caoIdx}(:));
            timAve(parak, anc_imgm_caoIdx, testk) = tEnd;
            countPairAve(parak, anc_imgm_caoIdx, testk) = numPairMatch;
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
savePath = sprintf('exp_rawmat_%s_%s_batch_%d.mat', target.config.database, target.config.category, target.config.batchSize);
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
