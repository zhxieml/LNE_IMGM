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

varyMinGrhCnt=21; varyMaxGrhCnt=50; grhTestCnt = 5;% 
target.config.database = 'synthetic';% only synthetic test is allowed here
target.config.Sacle_2D = 0.05;
iterRange = 6;
graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
algNameSepSpace = '                    ';
algSet.algNameSet = {'cao_pc','cao_pc_raw','imgm_d','imgm_r', 'tbimgm_cao', 'tbimgm_cao_pc', 'tbimgm_cao_uc', 'tbimgm_qm', 'tbimgm_matchALS'};
algSet.algEnable = [1,1,1,1,1,1,0,0,0];
algSet.algColor = {cao_pcClr,cao_pc_rawClr,imgm_dClr,imgm_rCls, tbimgm_caoClr, tbimgm_cao_pcClr, tbimgm_cao_ucClr, tbimgm_qmClr, tbimgm_matchALSClr};
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


target.config.initConstWeight = .2; % initial weight for consitency regularizer, suggest 0.2-0.25
target.config.constStep = 1.1;% inflate parameter, suggest 1.1-1.2
target.config.constWeightMax = 1;
target.config.constIterImmune = 2; % in early iterations, not involve consistency, suggest 1-3
target.config.edgeAffinityWeight = 1;% in random graphs, only edge affinity is used, angle is meaningless
target.config.angleAffinityWeight = 1 - target.config.edgeAffinityWeight;
target.config.selectNodeMask = 1:1:nInlier+target.config.nOutlier;
target.config.selectGraphMask{1} = 1:graphMaxCnt;
paraCnt = graphMaxCnt-graphMinCnt+1;


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
X=cell(algCnt,1);
target.config.nodeCnt = length(target.config.selectNodeMask);
% target.config.graphCnt = min(max(graphRange),length(target.config.selectGraphMask{1}));
target.config.graphCnt = graphMaxCnt;
nodeCnt = target.config.nodeCnt;
graphCnt = target.config.graphCnt;

% paraCnt: iterate over graph #
% algCnt: iterate over algorithms
% testCnt: iterate over tests
timAve = zeros(paraCnt,algCnt,testCnt);timAveFull = zeros(paraCnt,algCnt);
accAve = zeros(paraCnt,algCnt,testCnt);accAveFull = zeros(paraCnt,algCnt);
scrAve = zeros(paraCnt,algCnt,testCnt);scrAveFull = zeros(paraCnt,algCnt);
conPairAve = zeros(paraCnt,algCnt,testCnt);conPairAveFull = zeros(paraCnt,algCnt);
accStd = zeros(paraCnt,algCnt,testCnt);accStdFull = zeros(paraCnt,algCnt);
scrStd = zeros(paraCnt,algCnt,testCnt);scrStdFull = zeros(paraCnt,algCnt);
conPairStd = zeros(paraCnt,algCnt,testCnt);conPairStdFull = zeros(paraCnt,algCnt);

% conMatPairGraph = cell(paraCnt,algCnt,testCnt);
% accMatPairGraph = cell(paraCnt,algCnt,testCnt);
% scrMatPairGraph = cell(paraCnt,algCnt,testCnt);
acc = cell{1, algCnt};
scr = cell{1, algCnt};
con = cell{1, algCnt};
prevMatching = cell{1, algCnt};
increMatching = cell{1, algCnt};
matTmp = cell{1, algCnt};

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
    sigma = 0;
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

    for parak = 1:paraCnt % IMGM count

        param.n = nodeCnt; 
        param.N = graphMinCnt + parak - 1;
        
        scrDenomCurrent = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
        baseMat = CAO(rawMat(1:end-nodeCnt*paraCnt,1:end-nodeCnt*paraCnt),nodeCnt, graphMinCnt, iterRange,scrDenomCurrent, 'pair',1);
        
        %%%%%%%%%%%% calculate the incremental matching with cao_pc_raw %%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(cao_pc_rawIdx)
            tStart = tic;
            increMatching{cao_pc_rawIdx} = CAO(rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak)), nodeCnt, graphMinCnt+parak , iterRange, scrDenomCurrent, 'pair',1);
            tEnd = toc(tStart);
            
            acc{cao_pc_rawIdx} = cal_pair_graph_accuracy(increMatching{cao_pc_rawIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{cao_pc_rawIdx} = cal_pair_graph_score(increMatching{cao_pc_rawIdx},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{cao_pc_rawIdx} = cal_pair_graph_consistency(increMatching{cao_pc_rawIdx},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, cao_pc_rawIdx, testk) = mean(acc{cao_pc_rawIdx}(:));
            scrAve(parak, cao_pc_rawIdx, testk) = mean(scr{cao_pc_rawIdx}(:));
            conPairAve(parak, cao_pc_rawIdx, testk) = mean(con{cao_pc_rawIdx}(:));
            timAve(parak, cao_pc_rawIdx, testk) = tEnd;
        end

        %%%%%%%%%%% calculate the incremental matching with cao_pc %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % previous matching: caoPrevMatching
        if algSet.algEnable(cao_pcIdx)
            if parak ==1
                prevMatching{cao_pcIdx} = baseMat;
            end
            matTmp{cao_pcIdx} = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
            matTmp{cao_pcIdx}(1:end-nodeCnt,1:end-nodeCnt)=prevMatching{cao_pcIdx};
            tStart = tic;
            increMatching{cao_pcIdx} = CAO(matTmp{cao_pcIdx}, nodeCnt, graphMinCnt+parak , iterRange, scrDenomCurrent, 'pair',1);
            tEnd = toc(tStart);
            prevMatching{cao_pcIdx} = increMatching{cao_pcIdx};
            
            acc{cao_pcIdx} = cal_pair_graph_accuracy(increMatching{cao_pcIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{cao_pcIdx} = cal_pair_graph_score(increMatching{cao_pcIdx},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{cao_pcIdx} = cal_pair_graph_consistency(increMatching{cao_pcIdx},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, cao_pcIdx, testk) = mean(acc{cao_pcIdx}(:));
            scrAve(parak, cao_pcIdx, testk) = mean(scr{cao_pcIdx}(:));
            conPairAve(parak, cao_pcIdx, testk) = mean(con{cao_pcIdx}(:));
            timAve(parak, cao_pcIdx, testk) = tEnd;
        end

        %%%%%%%%%%% calculate the incremental matching with imgm_d %%%%%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(imgm_dIdx)
            if parak == 1
                prevMatching{imgm_dIdx} = baseMat;
            end
            matTmp{imgm_dIdx} = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
            matTmp{imgm_dIdx}(1:end-nodeCnt,1:end-nodeCnt)=prevMatching{imgm_dIdx};
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{imgm_dIdx},affinity.GT(1:nodeCnt*(graphMinCnt+parak),1:nodeCnt*(graphMinCnt+parak)),nodeCnt,graphMinCnt+parak,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{imgm_dIdx},nodeCnt,graphMinCnt+parak,0);
            

            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
            param.iterMax = iterRange;
            param.visualization = 0;
            param.method = 1; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_dIdx} = IMGM_old(simAP, matTmp{imgm_dIdx}, param);
            tEnd = toc(tStart);
            prevMatching{imgm_dIdx} = increMatching{imgm_dIdx};
            
            acc{imgm_dIdx} = cal_pair_graph_accuracy(increMatching{imgm_dIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{imgm_dIdx} = cal_pair_graph_score(increMatching{imgm_dIdx},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{imgm_dIdx} = cal_pair_graph_consistency(increMatching{imgm_dIdx},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, imgm_dIdx, testk) = mean(acc{imgm_dIdx}(:));
            scrAve(parak, imgm_dIdx, testk) = mean(scr{imgm_dIdx}(:));
            conPairAve(parak, imgm_dIdx, testk) = mean(con{imgm_dIdx}(:));
            timAve(parak, imgm_dIdx, testk) = tEnd;
        end

        %%%%%%%%%%% calculate the incremental matching with imgm_r %%%%%%%%%%%%%%%%%%%%%%%        
        if algSet.algEnable(imgm_dIdx)
            if parak == 1
                prevMatching{imgm_dIdx} = baseMat;
            end
            matTmp{imgm_dIdx} = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
            matTmp{imgm_dIdx}(1:end-nodeCnt,1:end-nodeCnt)=prevMatching{imgm_dIdx};
            
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{imgm_dIdx},affinity.GT(1:nodeCnt*(graphMinCnt+parak),1:nodeCnt*(graphMinCnt+parak)),nodeCnt,graphMinCnt+parak,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{imgm_dIdx},nodeCnt,graphMinCnt+parak,0);
            

            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
            param.iterMax = iterRange;
            param.visualization = 0;
            param.method = 3; % isDPP = 1; isAP = 2; isRand = 3; isTIP = 4;
            tStart = tic;
            increMatching{imgm_dIdx} = IMGM_old(simAP, matTmp{imgm_dIdx}, param);
            tEnd = toc(tStart);
            prevMatching{imgm_dIdx} = increMatching{imgm_dIdx};
            
            acc{imgm_dIdx} = cal_pair_graph_accuracy(increMatching{imgm_dIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{imgm_dIdx} = cal_pair_graph_score(increMatching{imgm_dIdx},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{imgm_dIdx} = cal_pair_graph_consistency(increMatching{imgm_dIdx},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, imgm_dIdx, testk) = mean(acc{imgm_dIdx}(:));
            scrAve(parak, imgm_dIdx, testk) = mean(scr{imgm_dIdx}(:));
            conPairAve(parak, imgm_dIdx, testk) = mean(con{imgm_dIdx}(:));
            timAve(parak, imgm_dIdx, testk) = tEnd;
        end

        %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_pc %%%%%%%%%%%%%%%%%%%%%%%%%
        % param of tbimgm
        param.bVerbose = 0;
        param.maxNumSearch = 25;
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
            matTmp{tbimgm_cao_pcIdx} = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
            matTmp{tbimgm_cao_pcIdx}(1:end-nodeCnt,1:end-nodeCnt)=prevMatching{tbimgm_cao_pcIdx};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_cao_pcIdx},affinity.GT(1:nodeCnt*(graphMinCnt+parak),1:nodeCnt*(graphMinCnt+parak)),nodeCnt,graphMinCnt+parak,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_cao_pcIdx},nodeCnt,graphMinCnt+parak,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));

            tStart = tic;
            increMatching{tbimgm_cao_pcIdx} = TBIMGM(affinity, simAP, matTmp{tbimgm_cao_pcIdx}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao_pcIdx} = increMatching{tbimgm_cao_pcIdx};
            
            acc{tbimgm_cao_pcIdx} = cal_pair_graph_accuracy(increMatching{tbimgm_cao_pcIdx},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{tbimgm_cao_pcIdx} = cal_pair_graph_score(increMatching{tbimgm_cao_pcIdx},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{tbimgm_cao_pcIdx} = cal_pair_graph_consistency(increMatching{tbimgm_cao_pcIdx},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, tbimgm_cao_pcIdx, testk) = mean(acc{tbimgm_cao_pcIdx}(:));
            scrAve(parak, tbimgm_cao_pcIdx, testk) = mean(scr{tbimgm_cao_pcIdx}(:));
            conPairAve(parak, tbimgm_cao_pcIdx, testk) = mean(con{tbimgm_cao_pcIdx}(:));
            timAve(parak, tbimgm_cao_pcIdx, testk) = tEnd;
        end

        %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_cao)
            % param for tbimgm_cao
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.useCstDecay = 1;
            param.subMethodParam.cstDecay  = 0.7;
            param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = 5;
            param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'afnty';
            param.subMethodParam.useCstInlier = 1;
            % previous matching: prevMatching{tbimgm_cao}
            if parak == 1
                prevMatching{tbimgm_cao} = baseMat;
            end
            matTmp{tbimgm_cao} = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
            matTmp{tbimgm_cao}(1:end-nodeCnt,1:end-nodeCnt)=prevMatching{tbimgm_cao};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_cao},affinity.GT(1:nodeCnt*(graphMinCnt+parak),1:nodeCnt*(graphMinCnt+parak)),nodeCnt,graphMinCnt+parak,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_cao},nodeCnt,graphMinCnt+parak,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));

            tStart = tic;
            increMatching{tbimgm_cao} = TBIMGM(affinity, simAP, matTmp{tbimgm_cao}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao} = increMatching{tbimgm_cao};
            
            acc{tbimgm_cao} = cal_pair_graph_accuracy(increMatching{tbimgm_cao},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{tbimgm_cao} = cal_pair_graph_score(increMatching{tbimgm_cao},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{tbimgm_cao} = cal_pair_graph_consistency(increMatching{tbimgm_cao},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, tbimgm_cao, testk) = mean(acc{tbimgm_cao}(:));
            scrAve(parak, tbimgm_cao, testk) = mean(scr{tbimgm_cao}(:));
            conPairAve(parak, tbimgm_cao, testk) = mean(con{tbimgm_cao}(:));
            timAve(parak, tbimgm_cao, testk) = tEnd;
        end

        %%%%%%%%%%%% calculate the incremental matching with tbimgm_cao_uc %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_cao_uc)
            % param for tbimgm_cao_uc
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.useCstDecay = 1;
            param.subMethodParam.cstDecay  = 0.7;
            param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = 5;
            param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'unary';
            param.subMethodParam.useCstInlier = 1;
            % previous matching: prevMatching{tbimgm_cao_uc}
            if parak == 1
                prevMatching{tbimgm_cao_uc} = baseMat;
            end
            matTmp{tbimgm_cao_uc} = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
            matTmp{tbimgm_cao_uc}(1:end-nodeCnt,1:end-nodeCnt)=prevMatching{tbimgm_cao_uc};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_cao_uc},affinity.GT(1:nodeCnt*(graphMinCnt+parak),1:nodeCnt*(graphMinCnt+parak)),nodeCnt,graphMinCnt+parak,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_cao_uc},nodeCnt,graphMinCnt+parak,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
            param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));

            tStart = tic;
            increMatching{tbimgm_cao_uc} = TBIMGM(affinity, simAP, matTmp{tbimgm_cao_uc}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_cao_uc} = increMatching{tbimgm_cao_uc};
            
            acc{tbimgm_cao_uc} = cal_pair_graph_accuracy(increMatching{tbimgm_cao_uc},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{tbimgm_cao_uc} = cal_pair_graph_score(increMatching{tbimgm_cao_uc},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{tbimgm_cao_uc} = cal_pair_graph_consistency(increMatching{tbimgm_cao_uc},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, tbimgm_cao_uc, testk) = mean(acc{tbimgm_cao_uc}(:));
            scrAve(parak, tbimgm_cao_uc, testk) = mean(scr{tbimgm_cao_uc}(:));
            conPairAve(parak, tbimgm_cao_uc, testk) = mean(con{tbimgm_cao_uc}(:));
            timAve(parak, tbimgm_cao_uc, testk) = tEnd;
        end
        %%%%%%%%%%%% calculate the incremental matching with tbimgm_qm %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_qm)
            % param for tbimgm_qm
            param.subMethodParam.name = 'quickmatch';
            % previous matching: prevMatching{tbimgm_qm}
            if parak == 1
                prevMatching{tbimgm_qm} = baseMat;
            end
            matTmp{tbimgm_qm} = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
            matTmp{tbimgm_qm}(1:end-nodeCnt,1:end-nodeCnt)=prevMatching{tbimgm_qm};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_qm},affinity.GT(1:nodeCnt*(graphMinCnt+parak),1:nodeCnt*(graphMinCnt+parak)),nodeCnt,graphMinCnt+parak,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_qm},nodeCnt,graphMinCnt+parak,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;

            tStart = tic;
            increMatching{tbimgm_qm} = TBIMGM(affinity, simAP, matTmp{tbimgm_qm}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_qm} = increMatching{tbimgm_qm};
            
            acc{tbimgm_qm} = cal_pair_graph_accuracy(increMatching{tbimgm_qm},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{tbimgm_qm} = cal_pair_graph_score(increMatching{tbimgm_qm},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{tbimgm_qm} = cal_pair_graph_consistency(increMatching{tbimgm_qm},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, tbimgm_qm, testk) = mean(acc{tbimgm_qm}(:));
            scrAve(parak, tbimgm_qm, testk) = mean(scr{tbimgm_qm}(:));
            conPairAve(parak, tbimgm_qm, testk) = mean(con{tbimgm_qm}(:));
            timAve(parak, tbimgm_qm, testk) = tEnd;
        end     
        %%%%%%%%%%%% calculate the incremental matching with tbimgm_matchALS %%%%%%%%%%%%%%%%%%%%
        if algSet.algEnable(tbimgm_matchALS)
            % param for tbimgm_matchALS
            param.subMethodParam.name = 'quickmatch';
            % previous matching: prevMatching{tbimgm_matchALS}
            if parak == 1
                prevMatching{tbimgm_matchALS} = baseMat;
            end
            matTmp{tbimgm_matchALS} = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
            matTmp{tbimgm_matchALS}(1:end-nodeCnt,1:end-nodeCnt)=prevMatching{tbimgm_matchALS};
            scrDenomMatInCntTmp = cal_pair_graph_inlier_score(matTmp{tbimgm_matchALS},affinity.GT(1:nodeCnt*(graphMinCnt+parak),1:nodeCnt*(graphMinCnt+parak)),nodeCnt,graphMinCnt+parak,nodeCnt);
            conDenomMatInCntTmp = cal_pair_graph_consistency(matTmp{tbimgm_matchALS},nodeCnt,graphMinCnt+parak,0);
            
            simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;

            tStart = tic;
            increMatching{tbimgm_matchALS} = TBIMGM(affinity, simAP, matTmp{tbimgm_matchALS}, param);
            tEnd = toc(tStart);
            prevMatching{tbimgm_matchALS} = increMatching{tbimgm_matchALS};
            
            acc{tbimgm_matchALS} = cal_pair_graph_accuracy(increMatching{tbimgm_matchALS},affinity.GT,target.config.nOutlier,nodeCnt,graphMinCnt+parak);
            scr{tbimgm_matchALS} = cal_pair_graph_score(increMatching{tbimgm_matchALS},affinity.GT,nodeCnt,graphMinCnt+parak);
            con{tbimgm_matchALS} = cal_pair_graph_consistency(increMatching{tbimgm_matchALS},nodeCnt,graphMinCnt+parak,0);
            
            accAve(parak, tbimgm_matchALS, testk) = mean(acc{tbimgm_matchALS}(:));
            scrAve(parak, tbimgm_matchALS, testk) = mean(scr{tbimgm_matchALS}(:));
            conPairAve(parak, tbimgm_matchALS, testk) = mean(con{tbimgm_matchALS}(:));
            timAve(parak, tbimgm_matchALS, testk) = tEnd;
        end    

        fprintf('The accuracy is %4.3f\n', accAve(parak, :, testk));
        fprintf('The time is %4.3f\n', timAve(parak, :, testk));
        fprintf('Start from %d graphs, %d graphs incremented\n', graphMinCnt, parak);
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
    fprintf(' acc   scr   con   tim   ');
    fprintf(fidPerf,', acc,  score, consis, time');
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
    end
    fprintf(' %02d,  %02d ',viewCnt,iterRange);fprintf(fidPerf,' %02d,  %02d',viewCnt,iterRange);
    for algk=1:algCnt
        if algSet.algEnable(algk)==0,continue;end
        fprintf('| %.3f %.3f %.3f %.3f',accAveFull(parak,algk),scrAveFull(parak,algk),conPairAveFull(parak,algk),timAveFull(parak,algk));
        fprintf(fidPerf,', %.3f, %.3f, %.3f, %.3f',accAveFull(parak,algk),scrAveFull(parak,algk),conPairAveFull(parak,algk),timAveFull(parak,algk));% fprintf(cc,  score, consis
    end
    fprintf('\n');fprintf(fidPerf,'\n');
end

legendOff = 0;
xtag='view #';ytag=['accuracy (',target.config.database,' ',target.config.category,')'];
plotResult(legendOff,graphRange,accAveFull,accStdFull,algSet,xtag,ytag);