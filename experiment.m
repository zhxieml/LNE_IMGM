%% Initialization
clear *;clear -global *;close all;clc; clear all;

global affinity;
global target;

nodeCnt = 10;
graphCnt = 50;
testCnt = 5;
nInlier = 10;
paraCnt = 8;
varyMinGrhCnt=4;varyMaxGrhCnt=32;grhTestCnt = 5;% 

target.config.database = 'synthetic';% only synthetic test is allowed here
target.config.testType = 'formal';% massOutlier

target.config.category = 'deform';%'deform','outlier','density','complete'
switch target.config.category
    case 'deform'% same setting with 5th row in Table 1 in the PAMI paper 
        nInlier = 10;target.config.nOutlier = 0;
        target.config.deform = 0.15; % Gaussian noise 
        target.config.density = .9; % density parameter 老
        target.config.complete = 1;
        graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
    case 'outlier'% same setting with 6th row in Table 1 in the PAMI paper 
        nInlier = 6;target.config.nOutlier = 4;
        target.config.deform = 0; % Gaussian noise 
        target.config.density = 1; % density parameter 老
        target.config.complete = 1;
        graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
    case 'density'% same setting with 7th row in Table 1 in the PAMI paper 
        nInlier = 10;target.config.nOutlier = 0;
        target.config.deform = 0.0; % Gaussian noise 
        target.config.density = 0.5; % density parameter 老
        target.config.complete = 1;
        graphMinCnt = varyMinGrhCnt;graphMaxCnt = varyMaxGrhCnt;testCnt = grhTestCnt;
    case 'complete'% same setting with 8th row in Table 1 in the PAMI paper 
        nInlier = 10;target.config.nOutlier = 0;
        target.config.deform = 0.05; % Gaussian noise 
        target.config.density = 1; % density parameter 老
        target.config.complete = 0.1;     
end

iterRange = 6;
graphRange = graphMinCnt:4:graphMaxCnt;

target.config.graphCnt = graphCnt;
target.config.nodeCnt = nodeCnt;
target.config.inCntType = 'exact'; % set 'exact' for "more outlier case", e.g. Fig.5 and Fig.6
target.config.Sacle_2D = 0.05; % 'Sacle' is a typo

% used for CAO
target.config.initConstWeight = .2; % initial weight for consitency regularizer, suggest 0.2-0.25
target.config.constStep = 1.1;% inflate parameter, suggest 1.1-1.2
target.config.constWeightMax = 1;
target.config.constIterImmune = 2; % in early iterations, not involve consistency, suggest 1-3

target.config.edgeAffinityWeight = 1;% in random graphs, only edge affinity is used, angle is meaningless
target.config.angleAffinityWeight = 1 - target.config.edgeAffinityWeight;
target.config.bGraphMatch = 0;% set to 1 use random graphs, otherwise use random points as set in the MPM code/paper
target.config.category = 'outlier';% only outlier are supported here

setObsoleteVariables;% some old parameters are used for debug and other tests, less relevant to the algorithm
algpar = setPairwiseSolver();

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

% Choose to run which part of the test
test_one_increment = false;
test_multiple_increment = true;

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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            param.method = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            scrDenomCurrent = max(max(scrDenomMatInCnt(1:param.N,1:param.N)));
            
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
            end
            
            %%%%%%%%%%%%%%%%%%%%% multiple incremental tests %%%%%%%%%%%%%%%%%%%%%
            if test_multiple_increment
                IMGMcount = 10;
                scrResult = zeros(3,IMGMcount);
                accResult = zeros(3, IMGMcount);
                conResult = zeros(3, IMGMcount);
                timeResult = zeros(3, IMGMcount);
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
                    imgmIncreMatching = IMGM_old(simAP, imgmMatTmp, param);
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
            %%%%%%%%%%%%%%%%%%%%% multiple incremental tests %%%%%%%%%%%%%%%%%%%%%


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
        end
    end
end