%% Initialization
clear *;clear -global *;close all;clc; clear all;

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


target.config.database = 'synthetic';% only synthetic test is allowed here
target.config.Sacle_2D = 0.05;
iterRange = 6;
graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;

algNameSepSpace = '                    ';
algSet.algColor = {caoClr,cao_rawClr,IMGMClr,IMGM_newClr};
algSet.algLineStyle = {'--',':','-','-.'};
algSet.algMarker = {'.','.','.','.'};
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



[~,caoIdx] = ismember('cao',algSet.algNameSet);
[~,rawCaoIdx] = ismember('cao_raw',algSet.algNameSet);
[~,IMGMIdx] = ismember('IMGM_old',algSet.algNameSet);
[~,IMGMNewIdx] = ismember('IMGM_cao',algSet.algNameSet);

algSet.algNameSetDisplay{caoIdx} = 'cao';
algSet.algNameSetDisplay{rawCaoIdx} = 'raw_cao';
algSet.algNameSetDisplay{IMGMIdx} = 'IMGM_old';
algSet.algNameSetDisplay{IMGMNewIdx} = 'IMGM_cao';

algCnt = length(algSet.algEnable); 
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
% fidPerf = fopen('results.csv','w');
% fprintf(fidPerf, 'testType,bGraphMatch,unaryFeat,edgeFeat,inCntType,database,category,testCnt,iter#,total node#,outlier#,complete,density,deform,graph#,alg#,scale,edgeWeight,initConstWeight,consStep,constWeightMax,iterImmune\n');
% fprintf(fidPerf, '%s,%d,%d,%d,%s,%s,%s,%d,%d,%d,%d,%.2f,%.2f,%.2f,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%d\n',...
%     target.config.testType,target.config.bGraphMatch,target.config.bUnaryEnable,target.config.bEdgeEnable,target.config.inCntType,target.config.database,target.config.category,...
%     testCnt,max(iterRange),nodeCnt,target.config.nOutlier,target.config.complete,target.config.density,target.config.deform,graphCnt,sum(algSet.algEnable),target.config.Sacle_2D,target.config.edgeAffinityWeight,...
%     target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
% fprintf('testType=%s, bGraphMatch=%d, unaryEnable=%d, edgeEnable=%d, inCntType=%s, database=%s, category=%s, iter#=%d, test#=%d, node#=%d, outlier#=%d,complete=%.2f, density=%.2f, deform=%.2f, graph#=%d, alg#=%d, edgeWeight=%.2f, scale=%.2f, initW=%.2f, stepW=%.2f, maxW=%.2f,iterImmune=%d\n',...
%     target.config.testType,target.config.bGraphMatch,target.config.bUnaryEnable,target.config.bEdgeEnable,target.config.inCntType,target.config.database,target.config.category,max(iterRange),testCnt,...
%     nodeCnt,target.config.nOutlier,target.config.complete,target.config.density,target.config.deform,graphCnt,sum(algSet.algEnable),...
%     target.config.edgeAffinityWeight,target.config.Sacle_2D,target.config.initConstWeight,target.config.constStep,target.config.constWeightMax,target.config.constIterImmune);
% fprintf('\n');fprintf(fidPerf,'\n');
% inlierAcc = zeros(paraCnt,testCnt);
% estErr = zeros(testCnt,1);


%% Run the test
%%%%%%%%%%%%%%%%%%  experiment 1 different incremental algorithms %%%%%%%%%%%%%%%%%%%%%
global affinity target
varyMinGrhCnt=21; varyMaxGrhCnt=50; testCnt = 5;% 
graphRange = graphMinCnt:1:graphMaxCnt;
paraCnt=length(graphRange);
for testk = 1:testCnt
    fprintf('Run test in round %d\n', testk);

    affinity = generateRandomAffinity(nInlier,testk); 
    affinity.GT = repmat(eye(nodeCnt,nodeCnt),graphCnt,graphCnt);

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

    for parak = 1:paraCnt

        param.n = nodeCnt; 
        param.N = baseGraphCnt + parak - 1;
        param.iterMax = iterRange;
        param.visualization = 0;
        param.method = 1;
        
        scrDenomCurrent = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
        
        % calculate the base matching with CAO for the 1st iter
        if parak == 1
            baseMat = CAO(rawMat(1:end-nodeCnt*paraCnt,1:end-nodeCnt*paraCnt),nodeCnt, baseGraphCnt, iterRange,scrDenomCurrent, 'pair',1);
        end
        %%%%%%%%%%%% calculate the incremental matching with raw CAO %%%%%%%%%%%%%%%%%%%%%
        tStart = tic;
        caoRawMatching = CAO(rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak)), nodeCnt, baseGraphCnt+parak , iterRange, scrDenomCurrent, 'pair',1);
        tEnd = toc(tStart);
        
        accRawCAO = cal_pair_graph_accuracy(caoRawMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+parak);
        scrRawCAO = cal_pair_graph_score(caoRawMatching,affinity.GT,nodeCnt,baseGraphCnt+parak);
        conRawCAO = cal_pair_graph_consistency(caoRawMatching,nodeCnt,baseGraphCnt+parak,0);
        
        accAve(parak, 1, testk) = mean(accRawCAO(:));
        scrAve(parak, 1, testk) = mean(scrRawCAO(:));
        conPairAve(parak, 1, testk) = mean(conRawCAO(:));
        timAve(parak, 1, testk) = tEnd;

        %%%%%%%%%%% calculate the incremental matching with CAO %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % previous matching: caoPrevMatching
        if parak ==1
            caoPrevMatching = baseMat;
        end
        caoMatTmp = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
        caoMatTmp(1:end-nodeCnt,1:end-nodeCnt)=caoPrevMatching;
        tStart = tic;
        caoIncreMatching = CAO(caoMatTmp, nodeCnt, baseGraphCnt+parak , iterRange, scrDenomCurrent, 'pair',1);
        tEnd = toc(tStart);
        caoPrevMatching = caoIncreMatching;
        
        accIncreCAO = cal_pair_graph_accuracy(caoIncreMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+parak);
        scrIncreCAO = cal_pair_graph_score(caoIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+parak);
        conIncreCAO = cal_pair_graph_consistency(caoIncreMatching,nodeCnt,baseGraphCnt+parak,0);
        
        accAve(parak, 2, testk) = mean(accIncreCAO(:));
        scrAve(parak, 2, testk)=mean(scrIncreCAO(:));
        conPairAve(parak, 2, testk)=mean(conIncreCAO(:));
        timAve(parak, 2, testk) = tEnd;
                
        %%%%%%%%%%% calculate the incremental matching with IMGM_old %%%%%%%%%%%%%%%%%%%%%%%
        % previous matching: imgmOldPrevMatching
        if parak == 1
            imgmOldPrevMatching = baseMat;
        end
        imgmOldMatTmp = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
        imgmOldMatTmp(1:end-nodeCnt,1:end-nodeCnt)=imgmOldPrevMatching;
        
        scrDenomMatInCntTmp = cal_pair_graph_inlier_score(imgmOldMatTmp,affinity.GT(1:nodeCnt*(baseGraphCnt+parak),1:nodeCnt*(baseGraphCnt+parak)),nodeCnt,baseGraphCnt+parak,nodeCnt);
        conDenomMatInCntTmp = cal_pair_graph_consistency(imgmOldMatTmp,nodeCnt,baseGraphCnt+parak,0);
        
        % AP inital arguments
        sigma = 0;
        simAP = (1-sigma)*scrDenomMatInCntTmp + sigma*conDenomMatInCntTmp;
        param.subMethodParam.scrDenom = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));

        
        tStart = tic;
        imgmOldIncreMatching = IMGM_old(simAP, imgmOldMatTmp, param);
        tEnd = toc(tStart);
        imgmOldPrevMatching = imgmOldIncreMatching;
        
        accImgmU = cal_pair_graph_accuracy(imgmOldIncreMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+parak);
        scrImgmU = cal_pair_graph_score(imgmOldIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+parak);
        conImgmU = cal_pair_graph_consistency(imgmOldIncreMatching,nodeCnt,baseGraphCnt+parak,0);
        
        accAve(parak, 3, testk)=mean(accImgmU(:));
        scrAve(parak, 3, testk)=mean(scrImgmU(:));
        conPairAve(parak, 3, testk)=mean(conImgmU(:));
        timAve(parak, 3, testk)=tEnd;
        

        %%%%%%%%%%%% calculate the incremental matching with IMGM_cao %%%%%%%%%%%%%%%%%%%%%%%%%
        % param of IMGM_cao
        param.bVerbose = 0;
        param.maxNumSearch = 25;
        % param for IMGM_cao
        param.subMethodParam.name = 'CAO';
        param.subMethodParam.useCstDecay = 1;
        param.subMethodParam.cstDecay  = 0.7;
        param.subMethodParam.useWeightedDecay  = 0;
        param.subMethodParam.iterMax = 5;
        param.subMethodParam.scrDenom = scrDenomCurrent;
        param.subMethodParam.optType = 'pair';
        param.subMethodParam.useCstInlier = 1;
        % previous matching: imgmPrevMatching
        if parak == 1
            imgmPrevMatching = baseMat;
        end
        imgmMatTmp = rawMat(1:end-nodeCnt*(paraCnt - parak),1:end-nodeCnt*(paraCnt - parak));
        imgmMatTmp(1:end-nodeCnt,1:end-nodeCnt)=imgmPrevMatching;
        tStart = tic;
        imgmIncreMatching = IMGM_single_step(affinity, simAP, imgmMatTmp, param);
        tEnd = toc(tStart);
        imgmPrevMatching = imgmIncreMatching;
        
        accImgm = cal_pair_graph_accuracy(imgmIncreMatching,affinity.GT,target.config.nOutlier,nodeCnt,baseGraphCnt+parak);
        scrImgm = cal_pair_graph_score(imgmIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+parak);
        conImgm = cal_pair_graph_consistency(imgmIncreMatching,nodeCnt,baseGraphCnt+parak,0);
        
        accAve(parak, 4, testk)=mean(accImgm(:));
        scrAve(parak, 4, testk)=mean(scrImgm(:));
        conPairAve(parak, 4, testk)=mean(conImgm(:));
        timAve(parak, 4, testk)=tEnd;

        fprintf('The accuracy is %4.3f\n', accAve(parak, :, testk));
        fprintf('The time is %4.3f\n', timAve(parak, :, testk));
        fprintf('Start from %d graphs, %d graphs incremented\n', baseGraphCnt, parak);
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