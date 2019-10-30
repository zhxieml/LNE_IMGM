%% Initialization
clear *;clear -global *;close all;clc; clear all;

init_path;

global affinity target
varyMinGrhCnt=4;varyMaxGrhCnt=32;grhTestCnt = 5;% 

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


target.config.database = 'synthetic';% only synthetic test is allowed here
target.config.Sacle_2D = 0.05;
iterRange = 6;
graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
if strcmp(target.config.testType,'massOutlier')% outlier test, see Fig.5(b) in the PAMI paper
    % cao cao_c cao_uc cao_pc are not used in massive outlier mode, because
    % no need to enforce consistency by post-step, thus we disable them:
    algNameSepSpace = '                ';
    algSet.algNameSet = {'mpm','rrwm','cao_s','cao_','cao_c_s','cao_c_','cao_uc_s','cao_uc_','cao_pc_s','cao_pc_','mOpt','mSync'};
    algSet.algEnable = [1,1,1,1,1,1,1,1,1,1,1,1];
    algSet.algColor = {mpmClr,rrwmClr,caoClr,caoClr,cao_cClr,cao_cClr,cao_ucClr,cao_ucClr,cao_pcClr,cao_pcClr,iccvClr,nipsClr};
    algSet.algLineStyle = {'--','--','--','-','--','-','--','-','--','-','-','-'};
    algSet.algMarker = {'.','.','.','.','.','.','.','.','.','.','.','.'};
    target.config.bGraphMatch = 0;% set to 1 use random graphs, otherwise use random points as set in the MPM code/paper
    target.config.category = 'outlier';% only outlier are supported here
    target.config.inCntType = 'exact';% set 'exact' for "more outlier case", e.g. Fig.5 and Fig.6
    nInlier = 6;target.config.nOutlier = 12;target.config.deform = .05;
    target.config.density = 1;target.config.complete = 1;
    graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
else
    algNameSepSpace = '                    ';
    algSet.algNameSet = {'cao','cao_raw','IMGM','IMGM_new'};
    algSet.algEnable = [1,1,1,1];
    algSet.algColor = {caoClr,cao_rawClr,IMGMClr,IMGM_newClr};
    algSet.algLineStyle = {'--',':','-','-.'};
    algSet.algMarker = {'.','.','.','.'};
%     algSet.algNameSet = {'rrwm','cao_','cao','cao_c_','cao_c','cao_uc_','cao_uc','cao_pc_','cao_pc','mOpt','mSync','IMGM_new'};
%     algSet.algEnable = [1,1,1,1,1,1,1,1,1,1,1,1];
%     algSet.algColor = {rrwmClr,caoClr,caoClr,cao_cClr,cao_cClr,cao_ucClr,cao_ucClr,cao_pcClr,cao_pcClr,iccvClr,nipsClr,};
%     algSet.algLineStyle = {'--','--','-','--','-','--','-','--','-','-','-'};
%     algSet.algMarker = {'.','.','.','.','.','.','.','.','.','.','.'};
    target.config.bGraphMatch = 1;
    target.config.inCntType = 'all';% set 'all' for "only a few outlier case", e.g. Fig.1&2&3&4
    target.config.category = 'deform';%'deform','outlier','density','complete'
    switch target.config.category
        case 'deform'% same setting with 5th row in Table 1 in the PAMI paper 
            nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.15;
            target.config.density = .9;target.config.complete = 1;
            graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
        case 'outlier'% same setting with 6th row in Table 1 in the PAMI paper 
            nInlier = 6;target.config.nOutlier = 4;target.config.deform = 0;
            target.config.density = 1;target.config.complete = 1;
            graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
        case 'density'% same setting with 7th row in Table 1 in the PAMI paper 
            nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.0;
            target.config.density = 0.5;target.config.complete = 1;
            graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
        case 'complete'% same setting with 8th row in Table 1 in the PAMI paper 
            nInlier = 10;target.config.nOutlier = 0;target.config.deform = 0.05;
            target.config.density = 1;target.config.complete = 0.1;     
    end
end
graphRange = graphMinCnt:4:graphMaxCnt;
target.config.initConstWeight = .2; % initial weight for consitency regularizer, suggest 0.2-0.25
target.config.constStep = 1.1;% inflate parameter, suggest 1.1-1.2
target.config.constWeightMax = 1;
target.config.constIterImmune = 2; % in early iterations, not involve consistency, suggest 1-3
target.config.edgeAffinityWeight = 1;% in random graphs, only edge affinity is used, angle is meaningless
target.config.angleAffinityWeight = 1 - target.config.edgeAffinityWeight;
target.config.selectNodeMask = 1:1:nInlier+target.config.nOutlier;
target.config.selectGraphMask{1} = 1:graphMaxCnt;
paraCnt=length(graphRange);
iterCnt = length(iterRange);

% [~,rrwmIdx] = ismember('rrwm',algSet.algNameSet);
% [~,mpmIdx] = ismember('mpm',algSet.algNameSet);
% [~,cao_Idx] = ismember('cao_',algSet.algNameSet);[~,cao_sIdx] = ismember('cao_s',algSet.algNameSet);
% [~,caoIdx] = ismember('cao',algSet.algNameSet);
% [~,cao_c_Idx] = ismember('cao_c_',algSet.algNameSet);[~,cao_c_sIdx] = ismember('cao_c_s',algSet.algNameSet);
% [~,cao_cIdx] = ismember('cao_c',algSet.algNameSet);
% [~,cao_uc_Idx] = ismember('cao_uc_',algSet.algNameSet);[~,cao_uc_sIdx] = ismember('cao_uc_s',algSet.algNameSet);
% [~,cao_ucIdx] = ismember('cao_uc',algSet.algNameSet);
% [~,cao_pc_Idx] = ismember('cao_pc_',algSet.algNameSet);[~,cao_pc_sIdx] = ismember('cao_pc_s',algSet.algNameSet);
% [~,cao_pcIdx] = ismember('cao_pc',algSet.algNameSet);
[~,caoIdx] = ismember('cao',algSet.algNameSet);
[~,rawCaoIdx] = ismember('cao_raw',algSet.algNameSet);
[~,IMGMIdx] = ismember('IMGM',algSet.algNameSet);
[~,IMGMNewIdx] = ismember('IMGM_new',algSet.algNameSet);

algSet.algNameSetDisplay{caoIdx} = 'cao';
algSet.algNameSetDisplay{rawCaoIdx} = 'raw_cao';
algSet.algNameSetDisplay{IMGMIdx} = 'IMGM';
algSet.algNameSetDisplay{IMGMNewIdx} = 'IMGM_new';

algCnt = length(algSet.algEnable); 
X=cell(algCnt,1);
target.config.nodeCnt = length(target.config.selectNodeMask);
% target.config.graphCnt = min(max(graphRange),length(target.config.selectGraphMask{1}));
target.config.graphCnt = 50;
nodeCnt = target.config.nodeCnt;
graphCnt = target.config.graphCnt;

% paraCnt: iterate over graph #
% algCnt: iterate over algorithms
% testCnt: iterate over tests
% timAve = zeros(paraCnt,algCnt,testCnt);timAveFull = zeros(paraCnt,algCnt);
% accAve = zeros(paraCnt,algCnt,testCnt);accAveFull = zeros(paraCnt,algCnt);
% scrAve = zeros(paraCnt,algCnt,testCnt);scrAveFull = zeros(paraCnt,algCnt);
% conPairAve = zeros(paraCnt,algCnt,testCnt);conPairAveFull = zeros(paraCnt,algCnt);
% accStd = zeros(paraCnt,algCnt,testCnt);accStdFull = zeros(paraCnt,algCnt);
% scrStd = zeros(paraCnt,algCnt,testCnt);scrStdFull = zeros(paraCnt,algCnt);
% conPairStd = zeros(paraCnt,algCnt,testCnt);conPairStdFull = zeros(paraCnt,algCnt);
% 
% conMatPairGraph = cell(paraCnt,algCnt,testCnt);
% accMatPairGraph = cell(paraCnt,algCnt,testCnt);
% scrMatPairGraph = cell(paraCnt,algCnt,testCnt);

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

test_one_increment = false;
test_multiple_increment = true;

algpar.algMethod = 'RRWM_old';

%% Run the test
for testk = 1:testCnt
    fprintf('Run test in round %d\n', testk);

    affinity = generateRandomAffinity(nInlier,testk); 
    affinity.GT = repmat(eye(nodeCnt,nodeCnt),graphCnt,graphCnt);

    % rrwm pairwise match, once for all graph pairs
    tStart = tic;
    rawMat = generatePairAssignment(algpar,nodeCnt,graphCnt,testk);
    rawTotalTime = toc(tStart);
    
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

    for parak = 1:paraCnt
        viewCnt = graphRange(parak);
        rawPairMatchTime = rawTotalTime*viewCnt*(viewCnt-1)/(graphCnt*(graphCnt-1));
        
        offsetStep = graphRange(end);%for synthetic test, no need inner loop for a given set of graphs of size viewCnt
        affinity.viewCnt = viewCnt;
        subParaCnt = length(1:offsetStep:graphCnt-viewCnt+1);
        
        for grhOffset = 1:offsetStep:graphCnt-viewCnt+1%this logic can be ignored in the demo code, it is used for real image data
            rawscopeX = (grhOffset-1)*nodeCnt+1:(grhOffset+viewCnt-1)*nodeCnt;
            rawScopeY = (grhOffset-1)*nodeCnt+1:(grhOffset+viewCnt-1)*nodeCnt;
            affinity.Xgt = affinity.GT(rawscopeX,rawScopeY);
            baseMat = rawMat(rawscopeX,rawScopeY);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             X{rrwmIdx} = baseMat;
%             timAve(parak,rrwmIdx,testk) = timAve(parak,rrwmIdx,testk) + rawPairMatchTime;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            scrDenom = max(max(scrDenomMatInCnt(1:viewCnt,1:viewCnt)));%used to normalize the affinity score in the objective
            tStart = tic;
            % the following function returns both consistency and affinity, here only consistency needed
            % it also does not discriminate between inliers and outliers
            singleGraphConstList = cal_single_graph_consistency_score(baseMat,nodeCnt,viewCnt);
            [C,refConstGraph] = max(singleGraphConstList);
            % given the reference graph r, first compute the unary
            % graph-wise consistency, then rank them into cstGrhList
            cstGrhList = rankGrhByConsistencyToRefGrh(baseMat,refConstGraph,nodeCnt,viewCnt);
            updGrhList = [cstGrhList,refConstGraph];%for consistency rank
            refGraph = updGrhList(end);
            cstCalTime = toc(tStart);
                                   
            param.n = 10; param.N = graphCnt - 1;
            param.iterMax = iterRange;
            param.visualization = 0;
            param.method = 1;
            
            scrDenomCurrent = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % param for IMGM_new
            param.bVerbose = 0;
            param.maxNumSearch = 25;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            param.subMethodParam.name = 'DPMC';
            param.subMethodParam.name = 'CAO';
            param.subMethodParam.useCstDecay = 1;
            param.subMethodParam.cstDecay  = 0.7;
            param.subMethodParam.useWeightedDecay  = 0;
            param.subMethodParam.iterMax = 5;
            param.subMethodParam.scrDenom = scrDenomCurrent;
            param.subMethodParam.optType = 'pair';
            param.subMethodParam.useCstInlier = 1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if test_one_increment
                % tStart = tic;
                currentMat = CAO(rawMat(1:end-nodeCnt,1:end-nodeCnt),nodeCnt, param.N, iterRange,scrDenomCurrent, 'pair',1);
                % timsCost = toc(tStart)
                
                affScore = scrDenomMatInCnt;
                rawMatTmp = rawMat;
                rawMatTmp(1:end-nodeCnt,1:end-nodeCnt)=currentMat;
                scrDenomMatInCntTmp = cal_pair_graph_inlier_score(rawMatTmp,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
                conDenomMatInCntTmp = cal_pair_graph_consistency(rawMatTmp,nodeCnt,graphCnt,0);
                %%%%%%%%%%%%%%%testing different similarity%%%%%%%%%%
                sigma = 0;
                simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tStart = tic;
                Xmatching = IMGM_old(simAP, rawMatTmp, param);
                timsCost = toc(tStart)
                accIMGM = cal_pair_graph_accuracy(Xmatching,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
                scrIMGM = cal_pair_graph_score(Xmatching,affinity.GT,nodeCnt,graphCnt);
                conIMGM = cal_pair_graph_consistency(Xmatching,nodeCnt,graphCnt,0);
                results = ['The results of IMGM on 50 graphs, accuracy:',num2str(mean(accIMGM(:))),', score:',num2str(mean(scrIMGM(:))),', consistency:',num2str(mean(conIMGM(:)))];
                disp(results);
%                 
%                 scrDenomCurrent = max(max(scrDenomMatInCnt(1:end,1:end)));
%                 tStart = tic;
%                 Xoriginal = CAO(rawMatTmp, nodeCnt, graphCnt, iterRange, scrDenomCurrent, 'pair',1);
%                 timsCost = toc(tStart)
%                 accCAO = cal_pair_graph_accuracy(Xoriginal,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
%                 scrCAO = cal_pair_graph_score(Xoriginal,affinity.GT,nodeCnt,graphCnt);
%                 conCAO = cal_pair_graph_consistency(Xoriginal,nodeCnt,graphCnt,0);
%                 results = ['The results of CAO-R on 50 graphs, accuracy:',num2str(mean(accCAO(:))),', score:',num2str(mean(scrCAO(:))),', consistency:',num2str(mean(conCAO(:)))];
%                 disp(results);
%                 
%                 tStart = tic;
%                 Xoriginal = CAO(rawMat, nodeCnt, graphCnt, iterRange, scrDenomCurrent, 'pair',1);
%                 timsCost = toc(tStart)
%                 accCAO = cal_pair_graph_accuracy(Xoriginal,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
%                 scrCAO = cal_pair_graph_score(Xoriginal,affinity.GT,nodeCnt,graphCnt);
%                 conCAO = cal_pair_graph_consistency(Xoriginal,nodeCnt,graphCnt,0);
%                 results = ['The results of CAO on 50 graphs, accuracy:',num2str(mean(accCAO(:))),', score:',num2str(mean(scrCAO(:))),', consistency:',num2str(mean(conCAO(:)))];
%                 disp(results);
                
                %%%%%% IMGM_new
%                 hyperGraph = Prim(simAP);
%                 hyperGraph(param.N:end, :) = 0;
%                 hyperGraph(:, param.N:end) = 0;
%                 
%                 tStart = tic;
%                 Xoriginal = IMGM_single_step(affinity, simAP, rawMatTmp, hyperGraph, param);
%                 timsCost = toc(tStart)
%                 accCAO = cal_pair_graph_accuracy(Xoriginal,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
%                 scrCAO = cal_pair_graph_score(Xoriginal,affinity.GT,nodeCnt,graphCnt);
%                 conCAO = cal_pair_graph_consistency(Xoriginal,nodeCnt,graphCnt,0);
%                 results = ['The results of IMGM_new on 50 graphs, accuracy:',num2str(mean(accCAO(:))),', score:',num2str(mean(scrCAO(:))),', consistency:',num2str(mean(conCAO(:)))];
%                 disp(results);
                hyperGraph = zeros(affinity.graphCnt, affinity.graphCnt);
                hyperGraph(1:param.N, 1:param.N) = Prim(simAP(1:param.N, 1:param.N));

                tStart = tic;
                Xoriginal = IMGM_single_step(affinity, simAP, rawMatTmp, param);
                timsCost = toc(tStart)
                accCAO = cal_pair_graph_accuracy(Xoriginal,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
                scrCAO = cal_pair_graph_score(Xoriginal,affinity.GT,nodeCnt,graphCnt);
                conCAO = cal_pair_graph_consistency(Xoriginal,nodeCnt,graphCnt,0);
                results = ['The results of IMGM_new on 50 graphs, accuracy:',num2str(mean(accCAO(:))),', score:',num2str(mean(scrCAO(:))),', consistency:',num2str(mean(conCAO(:)))];
                disp(results);
            end
            
            %%%%%%%%%%%%%%%%%%%%% multiple incremental tests %%%%%%%%%%%%%%%%%%%%%
            if test_multiple_increment
                IMGMcount = 15;
                scrResult = zeros(4,IMGMcount);
                accResult = zeros(4, IMGMcount);
                conResult = zeros(4, IMGMcount);
                timeResult = zeros(4, IMGMcount);
                baseGraphCnt = graphCnt - IMGMcount;

                for i = 1:IMGMcount
                    % calculate the base matching with CAO for the 1st iter
                    if i == 1
                        baseMat = CAO(rawMat(1:end-nodeCnt*IMGMcount,1:end-nodeCnt*IMGMcount),nodeCnt, baseGraphCnt, iterRange,scrDenomCurrent, 'pair',1);
                    end
                    % calcualte the raw matching with CAO
                    tStart = tic;
                    caoMat = CAO(rawMat, nodeCnt, graphCnt, iterRange, scrDenomCurrent, 'pair',1);
                    tEnd = toc(tStart);
                    
                    
                    % calculate the incremental matching with CAO
                    % previous matching: caoPrevMathcing
                    if i ==1
                        caoPrevMatching = baseMat;
                    end
                    caoMatTmp = rawMat(1:end-nodeCnt*(IMGMcount - i),1:end-nodeCnt*(IMGMcount - i));
                    caoMatTmp(1:end-nodeCnt,1:end-nodeCnt)=caoPrevMatching;
                    tStart = tic;
                    caoIncreMatching = CAO(caoMatTmp, nodeCnt, baseGraphCnt+i , iterRange, scrDenomCurrent, 'pair',1);
                    tEnd = toc(tStart);
                    caoPrevMatching = caoIncreMatching;
                    
                    accIncreCAO = cal_pair_graph_accuracy(caoIncreMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+i);
                    scrIncreCAO = cal_pair_graph_score(caoIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
                    conIncreCAO = cal_pair_graph_consistency(caoIncreMatching,nodeCnt,baseGraphCnt+i,0);
                    
                    accResult(2,i)=mean(accIncreCAO(:));
                    scrResult(2,i)=mean(scrIncreCAO(:));
                    conResult(2,i)=mean(conIncreCAO(:));
                    timeResult(2,i)=tEnd;
                    
                    % calculate the incremental matching with raw CAO
                    tStart = tic;
                    caoRawMatching = CAO(rawMat(1:end-nodeCnt*(IMGMcount - i),1:end-nodeCnt*(IMGMcount - i)), nodeCnt, baseGraphCnt+i , iterRange, scrDenomCurrent, 'pair',1);
                    tEnd = toc(tStart);
                    
                    accRawCAO = cal_pair_graph_accuracy(caoRawMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+i);
                    scrRawCAO = cal_pair_graph_score(caoRawMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
                    conRawCAO = cal_pair_graph_consistency(caoRawMatching,nodeCnt,baseGraphCnt+i,0);
                    
                    accResult(1,i)=mean(accRawCAO(:));
                    scrResult(1,i)=mean(scrRawCAO(:));
                    conResult(1,i)=mean(conRawCAO(:));
                    timeResult(1,i)=tEnd;
                    
                    % calculate the incremental matching with IMGM_old
                    % previous matching: imgmPrevMatching
                    if i == 1
                        imgmOldPrevMatching = baseMat;
                    end
                    imgmOldMatTmp = rawMat(1:end-nodeCnt*(IMGMcount - i),1:end-nodeCnt*(IMGMcount - i));
                    imgmOldMatTmp(1:end-nodeCnt,1:end-nodeCnt)=imgmOldPrevMatching;
                    
                    scrDenomMatInCntTmp = cal_pair_graph_inlier_score(imgmOldMatTmp,affinity.GT(1:nodeCnt*(baseGraphCnt+i),1:nodeCnt*(baseGraphCnt+i)),nodeCnt,baseGraphCnt+i,nodeCnt);
                    conDenomMatInCntTmp = cal_pair_graph_consistency(imgmOldMatTmp,nodeCnt,baseGraphCnt+i,0);
                    
                    %%%%%%%%%%%%%%%%%% AP inital arguments %%%%%%%%%%%%%%%%%%
                    sigma = 0;
                    simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
                    param.N = baseGraphCnt + i - 1;
                    
                    param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    tStart = tic;
                    imgmOldIncreMatching = IMGM_old(simAP, imgmOldMatTmp, param);
                    tEnd = toc(tStart);
                    imgmOldPrevMatching = imgmOldIncreMatching;
                    
                    accImgmU = cal_pair_graph_accuracy(imgmOldIncreMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+i);
                    scrImgmU = cal_pair_graph_score(imgmOldIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
                    conImgmU = cal_pair_graph_consistency(imgmOldIncreMatching,nodeCnt,baseGraphCnt+i,0);
                    
                    accResult(3,i)=mean(accImgmU(:));
                    scrResult(3,i)=mean(scrImgmU(:));
                    conResult(3,i)=mean(conImgmU(:));
                    timeResult(3,i)=tEnd;
                    
                    % calculate the incremental matching with IMGM
                    if i == 1
                        imgmPrevMatching = baseMat;
                    end
                    imgmMatTmp = rawMat(1:end-nodeCnt*(IMGMcount - i),1:end-nodeCnt*(IMGMcount - i));
                    imgmMatTmp(1:end-nodeCnt,1:end-nodeCnt)=imgmPrevMatching;
                    
                    tStart = tic;
                    imgmIncreMatching = IMGM_single_step(affinity, simAP, imgmMatTmp, param);
                    tEnd = toc(tStart);
                    imgmPrevMatching = imgmIncreMatching;
                    
                    accImgm = cal_pair_graph_accuracy(imgmIncreMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+i);
                    scrImgm = cal_pair_graph_score(imgmIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
                    conImgm = cal_pair_graph_consistency(imgmIncreMatching,nodeCnt,baseGraphCnt+i,0);
                    
                    accResult(4,i)=mean(accImgm(:));
                    scrResult(4,i)=mean(scrImgm(:));
                    conResult(4,i)=mean(conImgm(:));
                    timeResult(4,i)=tEnd;
                    
%                     a = 1;
                    
                    fprintf('The accuracy is %4.3f\n', accResult(:, i));
                    fprintf('The time is %4.3f\n', timeResult(:, i));
                    fprintf('Start from %d graphs, %d graphs incremented\n', baseGraphCnt, i);
                end
%                 figure
%                 plot(1:IMGMcount, accResult(1, :), '-o');

                plotResult_new(0, 1:IMGMcount, timeResult', algSet, 'increment #', 'time');
                plotResult_new(0, 1:IMGMcount, accResult', algSet, 'increment #', 'acc');
                plotResult_new(0, 1:IMGMcount, scrResult', algSet, 'increment #', 'src');
                plotResult_new(0, 1:IMGMcount, conResult', algSet, 'increment #', 'con');

            end
            %%%%%%%%%%%%%%%%%%%%% multiple incremental tests %%%%%%%%%%%%%%%%%%%%%
            
        end % for grhOffSet
    end % for paraCnt
end % for testCnt
