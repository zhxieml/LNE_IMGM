close all;
set_conf_img;
%% Params for test algorithms
set_param_GM;
set_alg;

global affinity target;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% target.config.nodeCnt = length(target.config.selectNodeMask);
% target.config.graphCnt = min(max(graphRange),length(target.config.selectGraphMask{1}));
% nodeCnt = target.config.nodeCnt;
% graphCnt = target.config.graphCnt;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cls = conf.class;

disp('== initialize an object model...');
listOfModelDataFile = dir(fullfile(conf.imgDir, cls, '*.png'));
listOfFeatureDataFile = dir(fullfile(pwd,'data', 'feature', cls, '*.mat'));


%% Initial visualization
% hAlg1 = figure('Name',sprintf('%s', alg(1).strName), 'NumberTitle', 'off'); 
% hFeat = figure('Name','model initialization', 'NumberTitle', 'off'); 

graphCnt = size(listOfFeatureDataFile, 1);
nodeCnt = 10;
target.config.graphCnt = graphCnt;
target.config.nodeCnt = nodeCnt;
target.config.inCnt = 10;
target.config.nOutlier = 0;
target.config.testType = 'formal';
iterRange = 6;

rawMat = zeros(nodeCnt*graphCnt, nodeCnt*graphCnt);
affinity.GT = zeros(nodeCnt*graphCnt, nodeCnt*graphCnt);

if 0
    for xview = 1:graphCnt-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        filePathName1 = fullfile(conf.imgDir, cls, listOfModelDataFile(xview).name );

        anno_filePathName1 = fullfile(conf.annoDir, cls, [ listOfModelDataFile(xview).name(1:end-4) '.mat' ]);
        load(anno_filePathName1, 'pts_coord');%% load the annotated data    

        % extract features
        viewInfo1 = extract_localfeatures_mcho( filePathName1, fparam, 'verbose', true );

        idx_sel1 = [];
        for i=1:size(pts_coord,2)
            pt1 = pts_coord(1:2,i);
            %idx1 = nearestneighbour(pt1, viewInfo1.frame(1:2,:), 'Radius', distClose);
            idx1 = nearestneighbour(pt1, viewInfo1.frame(1:2,:), 'NumberOfNeighbours', 1);
            idx_sel1 = [ idx_sel1 idx1 ];
        end
        %[ x_sel1 idx_sel1 ] = sel_pts_by_meanshift( viewInfo1.frame(1:2,:), 25);
        %idx_sel1
        viewInfo1.type = viewInfo1.type( idx_sel1 );
        viewInfo1.frame = viewInfo1.frame( :, idx_sel1);
        viewInfo1.desc = viewInfo1.desc( :, idx_sel1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     feature_filePathName1 = fullfile(listOfFeatureDataFile(xview).folder, listOfFeatureDataFile(xview).name);
    %     load(feature_filePathName1, 'cdata');
    %     viewInfo1 = cdata.view;
    %     clear cdata;

        for yview = xview+1:graphCnt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            filePathName2 = fullfile(conf.imgDir, cls, listOfModelDataFile(yview).name ); 

            anno_filePathName2 = fullfile(conf.annoDir, cls, [ listOfModelDataFile(yview).name(1:end-4) '.mat' ]);
            load(anno_filePathName2, 'pts_coord');%% load the annotated data

            % extract features
            viewInfo2 = extract_localfeatures_mcho( filePathName2, fparam, 'verbose', true );

            idx_sel2 = [];
            for j=1:size(pts_coord,2)
                pt2 = pts_coord(1:2,j);
                %idx1 = nearestneighbour(pt1, viewInfo1.frame(1:2,:), 'Radius', distClose);
                idx2 = nearestneighbour(pt2, viewInfo2.frame(1:2,:), 'NumberOfNeighbours', 1);
                idx_sel2 = [ idx_sel2 idx2 ];
            end
            %[ x_sel1 idx_sel1 ] = sel_pts_by_meanshift( viewInfo1.frame(1:2,:), 25);
            %idx_sel1
            viewInfo2.type = viewInfo2.type( idx_sel2 );
            viewInfo2.frame = viewInfo2.frame( :, idx_sel2);
            viewInfo2.desc = viewInfo2.desc( :, idx_sel2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         feature_filePathName2 = fullfile(listOfFeatureDataFile(yview).folder, listOfFeatureDataFile(yview).name);
    %         load(feature_filePathName2, 'cdata');
    %         viewInfo2 = cdata.view;
    %         clear cdata;

             %% make it into cdata
            cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;
            cdata.view(1) = viewInfo1;  cdata.view(2) = viewInfo2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cdata2.fparam = fparam;  cdata2.aparam = aparam;  cdata2.mparam = mparam;
            cdata2.view(1) = viewInfo2;  cdata2.view(2) = viewInfo1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% visualize it
    %         imgInput = appendimages( cdata.view(1).img, cdata.view(2).img );
    %         imgInput = double(imgInput)./255;
    %         offset_x = size(cdata.view(1).img,2);
    %         set(0, 'CurrentFigure', hFeat); clf;
    %         %figure('Name','model initialization', 'NumberTitle', 'off'); 
    %         imshow(rgb2gray(imgInput)); hold on;
    %         visualizeFeatures(cdata.view(1).frame, 'style', 'point', 'colorcode', 'g');
    %         visualizeFeatures(cdata.view(2).frame, 'style', 'point', 'colorcode', 'g', 'offset', [offset_x 0] );

            %% perform matching
            res = wrapper_matchingBasic( alg(1), cdata );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            res2 = wrapper_matchingBasic( alg(1), cdata2 );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     fprintf('%10s - score:%.5f (%.5f), %.3fsec \n', alg(1).strName, res.score_raw, res.score_GM, res.time);   
            fprintf('%10s - score:%.5f (%.5f), xxxsec \n', alg(1).strName, res.score_raw, res.score_GM);      

            %% generate raw matrix
    %         mask = target.pairwiseMask{testk}(xscope, yscope);
            affinity.K{xview, yview} = res.affinityMatrix;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            affinity.K{yview, xview} = res2.affinityMatrix;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            pairMatch = zeros(nodeCnt, nodeCnt);
            match = res.cand_matchlist(:,find(res.X));
            for pairIndex = 1:size(match, 2)
                node1 = match(1, pairIndex);
                node2 = match(2, pairIndex);
                pairMatch(node1, node2) = 1;
            end

            xscope = (xview-1)*nodeCnt+1:xview*nodeCnt;
            yscope = (yview-1)*nodeCnt+1:yview*nodeCnt;
            rawMat(xscope, yscope) = pairMatch;

            %% generate ground truth affinity.GT
            conf.affinity.GTDir = 'C:/Users/xzh87/OneDrive/Reference/Lab/code/IMGM/data/ground_truth';
            xname = listOfModelDataFile(xview).name(end-7:end-4);
            yname = listOfModelDataFile(yview).name(end-7:end-4);
            GTfile = fullfile(conf.affinity.GTDir, cls, [cls '_' xname '_' yname '.mat'] ); 
            load(GTfile, 'groundTruth');
            affinity.GT(xscope, yscope) = groundTruth.assign;

            %% visualize matching
    %         set(0, 'CurrentFigure', hAlg1); clf;
    %         %figure('Name',sprintf('%s', alg(1).strName), 'NumberTitle', 'off'); 
    %         imshow(rgb2gray(imgInput)); hold on;
    %         colorCode = makeColorCode(nnz(res.X));
    % 
    %         visualizeMatches( cdata.view(1).frame, cdata.view(2).frame, res.cand_matchlist(:,find(res.X)), ...
    %             'offset', [size(cdata.view(1).img,2) 0], 'style', 'frame', 'colorcode', colorCode, 'weight', res.X_raw(find(res.X)) ); 

    %         pause;
        end
    end

    rawMat = rawMat + rawMat' + eye(nodeCnt*graphCnt,nodeCnt*graphCnt);
    
%     save('rawdata.mat');
    save('rawMat.mat', 'rawMat');
    save('affinity.mat', 'affinity');
else
    load('rawMat.mat');
    load('affinity.mat');
end

%% Some param added according to the old IMGM
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

%% IMGM
IMGM_param.n = 10; IMGM_param.N = graphCnt - 1;
IMGM_param.iterMax = iterRange;
IMGM_param.visualization = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IMGM_param.method = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scrDenomMattarget.config.inCnt = cal_pair_graph_inlier_score(rawMat,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
scrDenomCurrent = max(max(scrDenomMattarget.config.inCnt(1:31,1:31)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some question %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentMat = CAO(rawMat(1:end-nodeCnt,1:end-nodeCnt), nodeCnt, graphCnt-1, iterRange,scrDenomCurrent, 'pair',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


affScore = scrDenomMattarget.config.inCnt;

rawMatTmp = rawMat;
rawMatTmp(1:end-nodeCnt,1:end-nodeCnt)=currentMat;
scrDenomMattarget.config.inCntTmp = cal_pair_graph_inlier_score(rawMatTmp,affinity.GT,nodeCnt,graphCnt,target.config.inCnt);
conDenomMattarget.config.inCntTmp = cal_pair_graph_consistency(rawMatTmp,nodeCnt,graphCnt,0);
%%%%%%%%%%%%%%%testing different similarity%%%%%%%%%%
sigma = 0;
simAP = (1-sigma)*scrDenomMattarget.config.inCntTmp + sigma*conDenomMattarget.config.inCntTmp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart = tic;
Xmatching = IMGM_old(simAP, rawMatTmp, IMGM_param);
timsCost = toc(tStart)
accIMGM = cal_pair_graph_accuracy(Xmatching,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
scrIMGM = cal_pair_graph_score(Xmatching,affinity.GT,nodeCnt,graphCnt);
conIMGM = cal_pair_graph_consistency(Xmatching,nodeCnt,graphCnt,0);
results = ['The results of IMGM, accuracy:',num2str(mean(accIMGM(:))),', score:',num2str(mean(scrIMGM(:))),', consistency:',num2str(mean(conIMGM(:)))];
disp(results);

scrDenomCurrent = max(max(scrDenomMattarget.config.inCnt(1:end,1:end)));
tStart = tic;
Xoriginal = CAO(rawMatTmp, nodeCnt, graphCnt, iterRange, scrDenomCurrent, 'pair',1);
timsCost = toc(tStart)
accCAO = cal_pair_graph_accuracy(Xoriginal,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
scrCAO = cal_pair_graph_score(Xoriginal,affinity.GT,nodeCnt,graphCnt);
conCAO = cal_pair_graph_consistency(Xoriginal,nodeCnt,graphCnt,0);
results = ['The results of CAO-R, accuracy:',num2str(mean(accCAO(:))),', score:',num2str(mean(scrCAO(:))),', consistency:',num2str(mean(conCAO(:)))];
disp(results);

tStart = tic;
Xoriginal = CAO(rawMat, nodeCnt, graphCnt, iterRange, scrDenomCurrent, 'pair',1);
timsCost = toc(tStart)
accCAO = cal_pair_graph_accuracy(Xoriginal,affinity.GT,target.config.nOutlier,nodeCnt,graphCnt);
scrCAO = cal_pair_graph_score(Xoriginal,affinity.GT,nodeCnt,graphCnt);
conCAO = cal_pair_graph_consistency(Xoriginal,nodeCnt,graphCnt,0);
results = ['The results of CAO, accuracy:',num2str(mean(accCAO(:))),', score:',num2str(mean(scrCAO(:))),', consistency:',num2str(mean(conCAO(:)))];
disp(results);