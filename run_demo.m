%    ****************************************** Notice **********************************
%    This demo contains and tests the algorithmic code of two following papers (and their conference versions):  
%    [1] alg1/2/3: Multi-Graph Matching via Affinity Optimization with Graduated Consistency Regularization
%    IEEE Transactions on Pattern Analysis and Machine Intelligence, 2015 
%    Conference version: Graduated Consistency-Regularized Optimization for Multi-Graph Matching, ECCV 2014  
% 
%    [2] matchOpt (mOpt): Consistency-Driven Alternating Optimization for Multigraph Matching: A Unified Approach
%    IEEE Transactions on Image Processing, 2015
%    Conference version: Joint Optimization for Consistent Multiple Graph Matching, ICCV 2013
% 
%    In this demo code, it reproduces (randomly) the bottom row of Fig.3 and Fig.5(b)
%    on synthetic dataset, hence it is self-contained and no dependency on other data sources
%    Code written by Junchi Yan, Shanghai Jiao Tong University, 2013-2015
%    Note the demo code is extracted from the raw code of the TPAMI and TIP implementation for unit test. 
%    Please cite the above papers, if you would like to use the code in your research
%   
%    Other relevant multi-graph/pairwise matching work compared in this demo code:
%    [3] matchSync (mSync): Pachauri et al. Solving the multi-way matching problem by permutation synchronization, NIPS 2013
%    [4] mpm: Cho et al. Finding matches in a haystack: a max-pooling strategy for graph matching in the presence of outliers, CVPR 2014
%    [5] matchLift (mLift): Chen et al. Near-optimal joint object matching via convex relaxation, ICML 2014 (not tested in this demo code -- refer to the raw code from authors)
%    **************************************************************************************
%    Before runing this code, set the following three places
%    1) set target.config.testType = formal/massOutlier, formal for Fig.3(fghi)
%    2) set varyMinGrhCnt=4;varyMaxGrhCnt=32;grhTestCnt = 5; as you want
%    3) set target.config.category, e.g. outlier/deform/density/complete for different settings
%    *************************************************************************************
clear *;clear -global *;close all;clc; clear all;
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

% varyMinGrhCnt=4;varyMaxGrhCnt=32;grhTestCnt = 20;

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
    algSet.algNameSet = {'rrwm','cao_','cao','cao_c_','cao_c','cao_uc_','cao_uc','cao_pc_','cao_pc','mOpt','mSync'};
    algSet.algEnable = [1,1,1,1,1,1,1,1,1,1,1];
    algSet.algColor = {rrwmClr,caoClr,caoClr,cao_cClr,cao_cClr,cao_ucClr,cao_ucClr,cao_pcClr,cao_pcClr,iccvClr,nipsClr};
    algSet.algLineStyle = {'--','--','-','--','-','--','-','--','-','-','-'};
    algSet.algMarker = {'.','.','.','.','.','.','.','.','.','.','.'};
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

[~,rrwmIdx] = ismember('rrwm',algSet.algNameSet);
[~,mpmIdx] = ismember('mpm',algSet.algNameSet);
[~,cao_Idx] = ismember('cao_',algSet.algNameSet);[~,cao_sIdx] = ismember('cao_s',algSet.algNameSet);
[~,caoIdx] = ismember('cao',algSet.algNameSet);
[~,cao_c_Idx] = ismember('cao_c_',algSet.algNameSet);[~,cao_c_sIdx] = ismember('cao_c_s',algSet.algNameSet);
[~,cao_cIdx] = ismember('cao_c',algSet.algNameSet);
[~,cao_uc_Idx] = ismember('cao_uc_',algSet.algNameSet);[~,cao_uc_sIdx] = ismember('cao_uc_s',algSet.algNameSet);
[~,cao_ucIdx] = ismember('cao_uc',algSet.algNameSet);
[~,cao_pc_Idx] = ismember('cao_pc_',algSet.algNameSet);[~,cao_pc_sIdx] = ismember('cao_pc_s',algSet.algNameSet);
[~,cao_pcIdx] = ismember('cao_pc',algSet.algNameSet);
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
timAve = zeros(paraCnt,algCnt,testCnt);timAveFull = zeros(paraCnt,algCnt);
accAve = zeros(paraCnt,algCnt,testCnt);accAveFull = zeros(paraCnt,algCnt);
scrAve = zeros(paraCnt,algCnt,testCnt);scrAveFull = zeros(paraCnt,algCnt);
conPairAve = zeros(paraCnt,algCnt,testCnt);conPairAveFull = zeros(paraCnt,algCnt);
accStd = zeros(paraCnt,algCnt,testCnt);accStdFull = zeros(paraCnt,algCnt);
scrStd = zeros(paraCnt,algCnt,testCnt);scrStdFull = zeros(paraCnt,algCnt);
conPairStd = zeros(paraCnt,algCnt,testCnt);conPairStdFull = zeros(paraCnt,algCnt);

conMatPairGraph = cell(paraCnt,algCnt,testCnt);
accMatPairGraph = cell(paraCnt,algCnt,testCnt);
scrMatPairGraph = cell(paraCnt,algCnt,testCnt);

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
% now start to test
estErr = zeros(testCnt,1);
for testk = 1:testCnt
    affinity = generateRandomAffinity(nInlier,testk);
    affinity.GT = repmat(eye(nodeCnt,nodeCnt),graphCnt,graphCnt);%just use identity matrix as grund truth matchings
    
    % rrwm pairwise match, once for all graph pairs
    tStart = tic;% compute all pairwise matchings at one time -> rawMat 
    rawMat = generatePairAssignment(algpar,nodeCnt,graphCnt,testk);% generate matchings by pairwise matching solver
    rawTotalTime = toc(tStart);
    
    % mpmMethod, once for all graph pairs
    [~,mpmIdx] = ismember('mpm',algSet.algNameSet);
    if mpmIdx>0&&algSet.algEnable(mpmIdx)
        tStart = tic;
        mpmMat = generatePairAssignment(mpmAlgPar,nodeCnt,graphCnt,testk);
        mpmTotalTime = toc(tStart);
    end
            
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
    for parak=1:paraCnt % test on various sizes of graph sets, 4,8,12,...,32
        viewCnt=graphRange(parak);
        rawPairMatchTime = rawTotalTime*viewCnt*(viewCnt-1)/(graphCnt*(graphCnt-1));
        offsetStep = graphRange(end);%for synthetic test, no need inner loop for a given set of graphs of size viewCnt
        affinity.viewCnt = viewCnt;
        subParaCnt = length(1:offsetStep:graphCnt-viewCnt+1);
        for grhOffset = 1:offsetStep:graphCnt-viewCnt+1%this logic can be ignored in the demo code, it is used for real image data
            rawscopeX = (grhOffset-1)*nodeCnt+1:(grhOffset+viewCnt-1)*nodeCnt;
            rawScopeY = (grhOffset-1)*nodeCnt+1:(grhOffset+viewCnt-1)*nodeCnt;
            affinity.Xgt = affinity.GT(rawscopeX,rawScopeY);
            baseMat = rawMat(rawscopeX,rawScopeY);
            
            X{rrwmIdx} = baseMat;
            timAve(parak,rrwmIdx,testk) = timAve(parak,rrwmIdx,testk) + rawPairMatchTime;
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
            
            % mpm method
            [~,mpmIdx] = ismember('mpm',algSet.algNameSet);
            if mpmIdx>0&&algSet.algEnable(mpmIdx)
                X{mpmIdx} = mpmMat(rawscopeX,rawScopeY);
                mpmPairMatchTime = mpmTotalTime*viewCnt*(viewCnt-1)/(graphCnt*(graphCnt-1));
                timAve(parak,mpmIdx,testk) = timAve(parak,mpmIdx,testk) + mpmPairMatchTime;
            end
            % nipsMethod
            [~,nipsIdx] = ismember('mSync',algSet.algNameSet);
            if nipsIdx>0&&algSet.algEnable(nipsIdx)
                tStart = tic;
                X{nipsIdx} = SynchronizePermute(baseMat,nodeCnt,viewCnt,'sync');
                timAve(parak,nipsIdx,testk) = timAve(parak,nipsIdx,testk) + toc(tStart) + rawPairMatchTime;
            end
            % iccvMethod
            [~,iccvIdx] = ismember('mOpt',algSet.algNameSet);
            if iccvIdx>0&&algSet.algEnable(iccvIdx) 
                tStart = tic;
                algpar.bPathSelect = 1;% a flag parameter used in the T-IP paper, see more details in the paper abut path selection discussion
                X{iccvIdx} = ConsistMultiMatch(updGrhList,nodeCnt,viewCnt,algpar,baseMat);
                timAve(parak,iccvIdx,testk) = timAve(parak,iccvIdx,testk) + toc(tStart) + rawPairMatchTime + cstCalTime;
            end
            % now the suite of algorithms proposed in the PAMI paper
            if cao_Idx>0&&algSet.algEnable(cao_Idx)% use consistency-driven inlier elicination, by default for few oulier case
                tStart = tic;% see Alg.1 in the pami paper, CAO: composition based affinity optimization 
                X{cao_Idx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'afnty',1);
                timeCost = toc(tStart) + rawPairMatchTime;timAve(parak,cao_Idx,testk) = timAve(parak,cao_Idx,testk) + timeCost; 
                if caoIdx&&algSet.algEnable(caoIdx),tStart = tic;
                    X{caoIdx} = SynchronizePermute(X{cao_Idx},nodeCnt,viewCnt,'sync');
                    timAve(parak,caoIdx,testk) = timAve(parak,caoIdx,testk) + toc(tStart) + timeCost;
                end
            end
            if cao_sIdx>0&&algSet.algEnable(cao_sIdx)% use affinity-driven inlier elicination
                tStart = tic;% see Alg.1 in the pami paper, CAO: composition based affinity optimization 
                X{cao_sIdx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'afnty',0);
                timeCost = toc(tStart) + rawPairMatchTime;timAve(parak,cao_sIdx,testk) = timAve(parak,cao_sIdx,testk) + timeCost; 
            end
            if cao_uc_Idx&&algSet.algEnable(cao_uc_Idx),tStart = tic;
                X{cao_uc_Idx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'unary',1);
                timAve(parak,cao_uc_Idx,testk) = timAve(parak,cao_uc_Idx,testk) + ...
                    toc(tStart) + rawTotalTime*viewCnt*(viewCnt-1)/(graphCnt*(graphCnt-1));
                if cao_ucIdx&&algSet.algEnable(cao_ucIdx)
                    X{cao_ucIdx} = SynchronizePermute(X{cao_uc_Idx},nodeCnt,viewCnt,'sync');
                    timAve(parak,cao_ucIdx,testk) = timAve(parak,cao_ucIdx,testk) + toc(tStart) + ...
                    rawTotalTime*viewCnt*(viewCnt-1)/(graphCnt*(graphCnt-1));
                end
            end
            if cao_uc_sIdx&&algSet.algEnable(cao_uc_sIdx),tStart = tic;
                X{cao_uc_sIdx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'unary',0);
                timAve(parak,cao_uc_sIdx,testk) = timAve(parak,cao_uc_sIdx,testk) + ...
                    toc(tStart) + rawTotalTime*viewCnt*(viewCnt-1)/(graphCnt*(graphCnt-1));
            end
            if cao_pc_Idx&&algSet.algEnable(cao_pc_Idx),tStart = tic;
                X{cao_pc_Idx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'pair',1);
                timeCost = toc(tStart) + rawPairMatchTime;
                timAve(parak,cao_pc_Idx,testk) = timAve(parak,cao_pc_Idx,testk) + timeCost;
                if cao_pcIdx&&algSet.algEnable(cao_pcIdx),tStart = tic;
                    X{cao_pcIdx} = SynchronizePermute(X{cao_pc_Idx},nodeCnt,viewCnt,'sync');
                    timAve(parak,cao_pcIdx,testk) = timAve(parak,cao_pcIdx,testk) + timeCost + toc(tStart);
                end
            end
            if cao_pc_sIdx&&algSet.algEnable(cao_pc_sIdx),tStart = tic;
                X{cao_pc_sIdx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'pair',0);
                timeCost = toc(tStart) + rawPairMatchTime;
                timAve(parak,cao_pc_sIdx,testk) = timAve(parak,cao_pc_sIdx,testk) + timeCost;
            end
            if cao_c_Idx&&algSet.algEnable(cao_c_Idx),tStart = tic;
                X{cao_c_Idx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'exact',1);
                timeCost = toc(tStart) + rawPairMatchTime;
                timAve(parak,cao_c_Idx,testk) = timAve(parak,cao_c_Idx,testk) + timeCost;
                if cao_cIdx&&algSet.algEnable(cao_cIdx),tStart = tic;
                    X{cao_cIdx} = SynchronizePermute(X{cao_c_Idx},nodeCnt,viewCnt,'sync');
                    timAve(parak,cao_cIdx,testk) = timAve(parak,cao_cIdx,testk) + toc(tStart) + timeCost;
                end
            end
            if cao_c_sIdx&&algSet.algEnable(cao_c_sIdx),tStart = tic;
                X{cao_c_sIdx} = CAO(baseMat,nodeCnt,viewCnt,iterRange,scrDenom,'exact',0);
                timeCost = toc(tStart) + rawPairMatchTime;
                timAve(parak,cao_c_sIdx,testk) = timAve(parak,cao_c_sIdx,testk) + timeCost;
            end
           %% IMGM
            if 1
                param.n = 10; param.N = graphCnt - 1;
                param.iterMax = iterRange;
                param.visualization = 0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                param.method = 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                scrDenomCurrent = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
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
                
                scrDenomCurrent = max(max(scrDenomMatInCnt(1:end,1:end)));
                tStart = tic;
                Xoriginal = CAO(rawMatTmp, nodeCnt, graphCnt, iterRange, scrDenomCurrent, 'pair',1);
                timsCost = toc(tStart)
                accCAO = cal_pair_graph_accuracy(Xoriginal,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
                scrCAO = cal_pair_graph_score(Xoriginal,affinity.GT,nodeCnt,graphCnt);
                conCAO = cal_pair_graph_consistency(Xoriginal,nodeCnt,graphCnt,0);
                results = ['The results of CAO-R on 50 graphs, accuracy:',num2str(mean(accCAO(:))),', score:',num2str(mean(scrCAO(:))),', consistency:',num2str(mean(conCAO(:)))];
                disp(results);
                
                tStart = tic;
                Xoriginal = CAO(rawMat, nodeCnt, graphCnt, iterRange, scrDenomCurrent, 'pair',1);
                timsCost = toc(tStart)
                accCAO = cal_pair_graph_accuracy(Xoriginal,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
                scrCAO = cal_pair_graph_score(Xoriginal,affinity.GT,nodeCnt,graphCnt);
                conCAO = cal_pair_graph_consistency(Xoriginal,nodeCnt,graphCnt,0);
                results = ['The results of CAO on 50 graphs, accuracy:',num2str(mean(accCAO(:))),', score:',num2str(mean(scrCAO(:))),', consistency:',num2str(mean(conCAO(:)))];
                disp(results);
                
%                 accRAW = cal_pair_graph_accuracy(rawMat,affinity.GT, target.config.nOutlier, nodeCnt, 32);
%                 scrRAW = cal_pair_graph_score(rawMat,affinity.GT,nodeCnt,32);
%                 conRAW = cal_pair_graph_consistency(rawMat,nodeCnt,32,0);
%                 results = ['The results of Raw on 32 graphs, accuracy:',num2str(mean(accRAW(:))),', score:',num2str(mean(scrRAW(:))),', consistency:',num2str(mean(conRAW(:)))];
%                 disp(results);
                
%                 tStart = tic;
%                 Xsync = SynchronizePermute(rawMat,nodeCnt,32,'sync');
%                 timsCost = toc(tStart)
%                 accSYN = cal_pair_graph_accuracy(Xsync,affinity.GT,target.config.nOutlier,nodeCnt,32);
%                 scrSYN = cal_pair_graph_score(Xsync,affinity.GT,nodeCnt,32);
%                 conSYN = cal_pair_graph_consistency(Xsync,nodeCnt,32,0);
%                 results = ['The results of SYNC on 32 graphs, accuracy:',num2str(mean(accSYN(:))),', score:',num2str(mean(scrSYN(:))),', consistency:',num2str(mean(conSYN(:)))];
%                 disp(results);
            end
            
            %% multiple incremental tests
            if 1 
                IMGMcount = 10;
                scrResult = zeros(3,IMGMcount);
                accResult = zeros(3, IMGMcount);
                conResult = zeros(3, IMGMcount);
                timeResult = zeros(3, IMGMcount);
                baseGraphCnt = graphCnt - IMGMcount;
                param.n = 10;
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
                    
                    % calculate the incremental matching with IMGM
                    % previous matching: imgmPrevMatching
                    if i == 1
                        imgmPrevMatching = baseMat;
                    end
                    imgmMatTmp = rawMat(1:end-nodeCnt*(IMGMcount - i),1:end-nodeCnt*(IMGMcount - i));
                    imgmMatTmp(1:end-nodeCnt,1:end-nodeCnt)=imgmPrevMatching;
                    
                    scrDenomMatInCntTmp = cal_pair_graph_inlier_score(imgmMatTmp,affinity.GT(1:nodeCnt*(baseGraphCnt+i),1:nodeCnt*(baseGraphCnt+i)),nodeCnt,baseGraphCnt+i,nodeCnt);
                    conDenomMatInCntTmp = cal_pair_graph_consistency(imgmMatTmp,nodeCnt,baseGraphCnt+i,0);
                    
                    %%%%%%%%%%%%%%%%%% AP inital arguments %%%%%%%%%%%%%%%%%%
                    sigma = 0;
                    simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
                    param.N = baseGraphCnt + i - 1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    tStart = tic;
                    imgmIncreMatching = IMGM(simAP, imgmMatTmp, param);
                    tEnd = toc(tStart);
                    imgmPrevMatching = imgmIncreMatching;
                    
                    accImgmU = cal_pair_graph_accuracy(imgmIncreMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+i);
                    scrImgmU = cal_pair_graph_score(imgmIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
                    conImgmU = cal_pair_graph_consistency(imgmIncreMatching,nodeCnt,baseGraphCnt+i,0);
                    
                    accResult(3,i)=mean(accImgmU(:));
                    scrResult(3,i)=mean(scrImgmU(:));
                    conResult(3,i)=mean(conImgmU(:));
                    timeResult(3,i)=tEnd;
                    
                    
                    
                    
                    a = 1;
                end
            end
            a = 1;
           %%
            % now compute the performance
            for algk = 1:algCnt
                if algSet.algEnable(algk)==0||isempty(X{algk}),continue;end
                txt = algSet.algNameSet(algk);
                % compute the accuracy, affinity score, consistency for each pair of graphs
                accMatPairGraph{parak,algk,testk} = cal_pair_graph_accuracy(X{algk},affinity.Xgt,target.config.nOutlier,nodeCnt,viewCnt);
                scrMatPairGraph{parak,algk,testk} = cal_pair_graph_score(X{algk},affinity.Xgt,nodeCnt,viewCnt);
                conMatPairGraph{parak,algk,testk} = cal_pair_graph_consistency(X{algk},nodeCnt,viewCnt,0);
                % the overall accuracy and affinity score, note the diagnal
                % of xxxMatPairGraph is dismissed because they are meaningless
                accAve(parak,algk,testk) = accAve(parak,algk,testk)+mean(accMatPairGraph{parak,algk,testk}(logical(triu(ones(size(accMatPairGraph{parak,algk,testk})),1))));
                accStd(parak,algk,testk) = accStd(parak,algk,testk)+std(accMatPairGraph{parak,algk,testk}(logical(triu(ones(size(accMatPairGraph{parak,algk,testk})),1))));
                scrAve(parak,algk,testk) = scrAve(parak,algk,testk)+mean(scrMatPairGraph{parak,algk,testk}(logical(triu(ones(size(accMatPairGraph{parak,algk,testk})),1))));%mean(scrtmp(:,parak,algk));
                scrStd(parak,algk,testk) = scrStd(parak,algk,testk)+std(scrMatPairGraph{parak,algk,testk}(logical(triu(ones(size(accMatPairGraph{parak,algk,testk})),1))));%std(scrtmp(:,parak,algk));
                conPairAve(parak,algk,testk) = conPairAve(parak,algk,testk)+mean(conMatPairGraph{parak,algk,testk}(logical(triu(ones(size(accMatPairGraph{parak,algk,testk})),1))));
                conPairStd(parak,algk,testk) = conPairStd(parak,algk,testk)+std(conMatPairGraph{parak,algk,testk}(logical(triu(ones(size(accMatPairGraph{parak,algk,testk})),1))));
            end %for algk
        end % for grhOffSet
        accAve(parak,:,testk) = accAve(parak,:,testk)/subParaCnt;
        scrAve(parak,:,testk) = scrAve(parak,:,testk)/subParaCnt;
        accStd(parak,:,testk) = accStd(parak,:,testk)/subParaCnt;
        scrStd(parak,:,testk) = scrStd(parak,:,testk)/subParaCnt;
        conPairAve(parak,:,testk) = conPairAve(parak,:,testk)/subParaCnt;
        conPairStd(parak,:,testk) =conPairStd(parak,:,testk)/subParaCnt;
        timAve(parak,:,testk) = timAve(parak,:,testk)/subParaCnt;
        inlierAcc(parak,testk) = inlierAcc(parak,testk)/subParaCnt;
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

% print out the performance for different algorithms

algNamePreSpace = '                    ';
fprintf('--------------------------------------------------------------test %02d performance-------------------------------------------------------------------\n',testk);
fprintf(fidPerf,'test%02d\n',testk);
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
for parak=1:paraCnt%graph #
    viewCnt=graphRange(parak);
    fprintf(' %02d   %02d ',viewCnt,iterRange);
    fprintf(fidPerf,' %02d,  %02d',viewCnt,iterRange);
    for algk=1:algCnt
        if algSet.algEnable(algk)==0,continue;end
        fprintf('| %.3f %.3f %.3f %.3f',accAve(parak,algk,testk),scrAve(parak,algk,testk),conPairAve(parak,algk,testk),timAve(parak,algk,testk));
        fprintf(fidPerf,', %.3f, %.3f, %.3f, %.3f',accAve(parak,algk,testk),scrAve(parak,algk,testk),conPairAve(parak,algk,testk),timAve(parak,algk,testk));% fprintf(cc,  score, consis
    end
    fprintf('\n');
    fprintf(fidPerf,'\n');
end
end % for testCnt
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
    viewCnt=graphRange(parak);
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