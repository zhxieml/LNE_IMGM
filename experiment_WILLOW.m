%% Initialization
clear *;clear -global *;close all;clc; clear all;

init_path;

global affinity target

setPlotColor;

dataDir = 'C:\Users\xzh87\OneDrive\Reference\Lab\code\WILLOW-ObjectClass_dataset';
affinityDir=[dataDir '/affinity'];
resDir = './res';

testCnt = 5;
class = 'Duck';
graphCnt = 50;
inCnt = 10;
outCnt = 0;
nodeCnt = inCnt + outCnt;

incrementCnt = 20;
scrResult = zeros(testCnt, 4,incrementCnt);
accResult = zeros(testCnt, 4, incrementCnt);
conResult = zeros(testCnt, 4, incrementCnt);
timeResult = zeros(testCnt, 4, incrementCnt);
baseGraphCnt = graphCnt - incrementCnt;

% parameter for CAO
iterRange = 6;
target.config.inCnt = inCnt;
target.config.testType = 'formal';
target.config.constIterImmune = 2;% immune from consistency regularization, only use affinity in earl
target.config.constStep = 1.1;% the inflate parameter, e.g. 1.05-1.1
target.config.initConstWeight = 0.2;% initial consistency regularizer weight, e.g 0.2-0.25
target.config.constWeightMax = 1;% the upperbound, always set to 1
target.config.graphCnt = graphCnt;
target.config.nodeCnt = nodeCnt;

% load the data
loadPath = fullfile(affinityDir, class, sprintf('in%02d_out%02d.mat', inCnt, outCnt));
load(loadPath, 'rawMat', 'affinity');

% parameter for IMGM_old
affinity.BiDir = 1;
affinity.edgeAffinityWeight = 1;
affinity.angleAffinityWeight = 0;
affinity.nodeCnt = nodeCnt;
affinity.graphCnt = graphCnt;

affinity.EG = cell(graphCnt, 1);
[r,c] = find(ones(nodeCnt, nodeCnt));
for i = 1:graphCnt
    affinity.EG{i} = [r, c]';
    affinity.nP{i} = nodeCnt;
    affinity.nE{i} = nodeCnt * nodeCnt;
end

target.pairwiseMask = cell(1);
target.pairwiseMask{1} = ones(graphCnt*nodeCnt,graphCnt*nodeCnt);

%% Run the test
for testk = 1:testCnt
    fprintf('Run test in round %d\n', testk);
    
    scrDenomMatInCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,inCnt);
    scrDenomMatInCntGT = cal_pair_graph_inlier_score(affinity.GT,affinity.GT,nodeCnt,graphCnt,inCnt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    scrDenomCurrent = max(max(scrDenomMatInCnt(1:end,1:end)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     param.n = 10; param.N = graphCnt - 1;
%     param.iterMax = iterRange;
%     param.visualization = 0;
%     param.method = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % param for IMGM_new
    param.bVerbose = 0;
    param.maxNumSearch = 25;
    param.n = 10;
    param.iterMax = iterRange;
    param.visualization = 0;
    param.method = 1;
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
    
    for i = 1:incrementCnt
        if i == 1
            baseMat = CAO(rawMat(1:end-nodeCnt*incrementCnt,1:end-nodeCnt*incrementCnt), nodeCnt, baseGraphCnt, iterRange, scrDenomCurrent, 'pair',1);
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
        caoMatTmp = rawMat(1:end-nodeCnt*(incrementCnt - i),1:end-nodeCnt*(incrementCnt - i));
        caoMatTmp(1:end-nodeCnt,1:end-nodeCnt)=caoPrevMatching;
        tStart = tic;
        caoIncreMatching = CAO(caoMatTmp, nodeCnt, baseGraphCnt+i, iterRange, scrDenomCurrent, 'pair',1);
        tEnd = toc(tStart);
        caoPrevMatching = caoIncreMatching;
        
        accIncreCAO = cal_pair_graph_accuracy(caoIncreMatching,affinity.GT,outCnt,nodeCnt,baseGraphCnt+i);
        scrIncreCAO = cal_pair_graph_score(caoIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
        conIncreCAO = cal_pair_graph_consistency(caoIncreMatching,nodeCnt,baseGraphCnt+i,0);
        
        accResult(testk,1,i)=mean(accIncreCAO(:));
        scrResult(testk,1,i)=mean(scrIncreCAO(:));
        conResult(testk,1,i)=mean(conIncreCAO(:));
        timeResult(testk,1,i)=tEnd;
        
        % calculate the incremental matching with raw CAO
        tStart = tic;
        caoRawMatching = CAO(rawMat(1:end-nodeCnt*(incrementCnt - i),1:end-nodeCnt*(incrementCnt - i)), nodeCnt, baseGraphCnt+i , iterRange, scrDenomCurrent, 'pair',1);
        tEnd = toc(tStart);
        
        accRawCAO = cal_pair_graph_accuracy(caoRawMatching,affinity.GT,outCnt,nodeCnt,baseGraphCnt+i);
        scrRawCAO = cal_pair_graph_score(caoRawMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
        conRawCAO = cal_pair_graph_consistency(caoRawMatching,nodeCnt,baseGraphCnt+i,0);
        
        accResult(testk,2,i)=mean(accRawCAO(:));
        scrResult(testk,2,i)=mean(scrRawCAO(:));
        conResult(testk,2,i)=mean(conRawCAO(:));
        timeResult(testk,2,i)=tEnd;
        
        % calculate the incremental matching with IMGM_old
        % previous matching: imgmPrevMatching
        if i == 1
            imgmOldPrevMatching = baseMat;
        end
        imgmOldMatTmp = rawMat(1:end-nodeCnt*(incrementCnt - i),1:end-nodeCnt*(incrementCnt - i));
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
        
        accImgmU = cal_pair_graph_accuracy(imgmOldIncreMatching,affinity.GT,outCnt,nodeCnt,baseGraphCnt+i);
        scrImgmU = cal_pair_graph_score(imgmOldIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
        conImgmU = cal_pair_graph_consistency(imgmOldIncreMatching,nodeCnt,baseGraphCnt+i,0);
        
        accResult(testk,3,i)=mean(accImgmU(:));
        scrResult(testk,3,i)=mean(scrImgmU(:));
        conResult(testk,3,i)=mean(conImgmU(:));
        timeResult(testk,3,i)=tEnd;
        
        % calculate the incremental matching with IMGM
        if i == 1
            imgmPrevMatching = baseMat;
        end
        imgmMatTmp = rawMat(1:end-nodeCnt*(incrementCnt - i),1:end-nodeCnt*(incrementCnt - i));
        imgmMatTmp(1:end-nodeCnt,1:end-nodeCnt)=imgmPrevMatching;
        
        tStart = tic;
        imgmIncreMatching = IMGM_single_step(affinity, simAP, imgmMatTmp, param);
        tEnd = toc(tStart);
        imgmPrevMatching = imgmIncreMatching;
        
        accImgm = cal_pair_graph_accuracy(imgmIncreMatching,affinity.GT,outCnt,nodeCnt,baseGraphCnt+i);
        scrImgm = cal_pair_graph_score(imgmIncreMatching,affinity.GT,nodeCnt,baseGraphCnt+i);
        conImgm = cal_pair_graph_consistency(imgmIncreMatching,nodeCnt,baseGraphCnt+i,0);
        
        accResult(testk,4,i)=mean(accImgm(:));
        scrResult(testk,4,i)=mean(scrImgm(:));
        conResult(testk,4,i)=mean(conImgm(:));
        timeResult(testk,4,i)=tEnd;
                
        fprintf('The accuracy is %4.3f\n', accResult(testk, :, i));
        fprintf('The time is %4.3f\n', timeResult(testk, :, i));
        fprintf('Start from %d graphs, %d graphs incremented\n', baseGraphCnt, i);
    end
    
    plotResult_new(0, 1:incrementCnt, timeResult(testk, :, :)', algSet, 'increment #', 'time');
    plotResult_new(0, 1:incrementCnt, accResult(testk, :, :)', algSet, 'increment #', 'acc');
    plotResult_new(0, 1:incrementCnt, scrResult(testk, :, :)', algSet, 'increment #', 'src');
    plotResult_new(0, 1:incrementCnt, conResult(testk, :, :)', algSet, 'increment #', 'con');
    
    pause;

end


