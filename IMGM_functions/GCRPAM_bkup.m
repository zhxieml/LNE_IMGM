function [P,corr,nodeCon,inCnt] = CAO(rawMat,nodeCnt,graphCnt,iterMax,scrDenom,optType)
% Input:
% rawMat is a matrix of size nodeCnt*graphCnt \times nodeCnt*graphCnt.
% rawMat contains the pairwise matchings for all graph pairs computed by a
% pairwise matching solver, such as RRWM.
% The values of optType include (refer to PAMI paper): 'unary' (CAO-UC),
% 'score' (CAO), 'pair' (CAO-PC), 'exact' (CAO-C).
% scrDenom is used to normalize affinity score, in order to be comparable
% to consistency metrics in [0,1].
% if no iteration needed, just return the raw and initial matchings
if iterMax==0,P = rawMat;return;end
global affinity
global nodeConsistencyMask
global target
global matchMat
useCstInlier = target.config.useCstInlier;%set to 1 if use consistency to ,否则用affinity来甄别
boostSampleRate= target.config.boostSampleRate;
inCnt = target.config.inCnt;% the number of inliers, can be specified to different values by manual or automatically 
corr = zeros(iterMax,12);
% 下面两个标签只有在score max时才打开
nodeConsistencyMask =  zeros(nodeCnt,graphCnt);
if strcmp(target.config.testType,'massOutlier'),massOutlierMode = 1;else massOutlierMode = 0;end
calPairCorrFlag = 0;%overall graph level pairwise consistency/score correlation
calAnchCorrFlag = target.config.calAnchCorrFlag;%anchor graph level pairwise consistency/score correlation
constIterImune = target.config.constIterImune;%前constIterImune次专注于提升score
curIdx = 1;
lastIdx = 2;
constStep = target.config.constStep;
initConstWeight = target.config.initConstWeight;
constWeightMax = target.config.constWeightMax;
if massOutlierMode
    if useCstInlier    
        nodeConsistencyMask = cal_node_consistency_mask(rawMat,nodeCnt,graphCnt,inCnt);
    else
        nodeConsistencyMask = cal_node_affinity_mask(rawMat,nodeCnt,graphCnt,inCnt);
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

iter = 0;
err = ones(iterMax,1);
I = eye(nodeCnt*graphCnt,nodeCnt*graphCnt);
matchMat = rawMat;
% permMat = generatePermMatFromMatchMat(matchMat,nodeCnt,graphCnt);
constWeight = initConstWeight;
%这里已经考虑了inlier eliciting的mechanism
unaryListConsistency{curIdx} = cal_single_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode);
unaryListConsistency{lastIdx} = unaryListConsistency{curIdx};
%这里已经考虑了inlier eliciting的mechanism
pairListConsistency{curIdx} = cal_pair_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode);
pairListConsistency{lastIdx} = pairListConsistency{curIdx};

[X,Y] = meshgrid(1:graphCnt,1:graphCnt);
X = X(:);Y = Y(:);
while iter<iterMax
    tempMat = zeros(nodeCnt*graphCnt,nodeCnt*graphCnt); 
anchCorrCnt = 0;anchCorrSum=zeros(10,1);
for vk = 1:graphCnt^2% optimize to one-layer loop
    xview = X(vk); yview = Y(vk);
    if mod(vk-1,graphCnt) == 0,xscope = (xview-1)*nodeCnt+1:xview*nodeCnt;end
    if X(vk)>= Y(vk),continue;end
    yscope = (yview-1)*nodeCnt+1:yview*nodeCnt;
    if strcmp(optType,'score')%只考虑score anchCorr = [scr-acc,cst-acc,scr]';
        [tempMat(xscope,yscope),anchCorr] = find1stOrderPathByScoreGraphPairConsistency(...
                    xview,yview,unaryListConsistency{curIdx},pairListConsistency{curIdx},nodeCnt,graphCnt,constWeight,scrDenom,'score',massOutlierMode,boostSampleRate,calAnchCorrFlag);
        if sum(isnan(anchCorr))==0,anchCorrSum = anchCorrSum+anchCorr;anchCorrCnt = anchCorrCnt+1;end
    elseif strcmp(optType,'cstcy')%只考虑consistency
        tempMat(xscope,yscope) = find1stOrderPathByScoreGraphPairConsistency(...
                    xview,yview,unaryListConsistency{curIdx},pairListConsistency{curIdx},nodeCnt,graphCnt,constWeight,scrDenom,'cstcy',massOutlierMode,boostSampleRate,calAnchCorrFlag);
    else%exact模式，score和consistency均考虑，分两个阶段
        if iter>=constIterImune%兼顾score和consistency 超过一定迭代次数才启动consistency
            tempMat(xscope,yscope) = find1stOrderPathByScoreGraphPairConsistency(...
                    xview,yview,unaryListConsistency{curIdx},pairListConsistency{curIdx},nodeCnt,graphCnt,constWeight,scrDenom,optType,massOutlierMode,boostSampleRate,calAnchCorrFlag);
        else%imune 阶段靠score
            tempMat(xscope,yscope) = find1stOrderPathByScoreGraphPairConsistency(...
                    xview,yview,unaryListConsistency{curIdx},pairListConsistency{curIdx},nodeCnt,graphCnt,constWeight,scrDenom,'score',massOutlierMode,boostSampleRate,calAnchCorrFlag);
        end
    end
end%for vk
%     if exist('matchMatOld','var'),matchMatOldOld = matchMatOld;end
    matchMatOld = matchMat;
    matchMat = tempMat+tempMat'+I;%matchMat被更新了一次
    if iter>=constIterImune
        constWeight = min([constWeightMax,constStep*constWeight]);
    else
%         if iter == constIterImune-1 
%             if strcmp(target.config.inCntType,'est')%更新inCnt
%             [~,inCnt] = cal_single_node_consistency(matchMat,nodeCnt,graphCnt);
             % 更新和计算node-wise consistency只在纯score的时候进行
            if massOutlierMode
                if useCstInlier    
                    nodeConsistencyMask = cal_node_consistency_mask(matchMat,nodeCnt,graphCnt,inCnt);
                else
                    nodeConsistencyMask = cal_node_affinity_mask(matchMat,nodeCnt,graphCnt,inCnt);
                end
            end
%         end
    end
    iter = iter + 1;
    err(iter) = sum(abs(matchMatOld(:)-matchMat(:)));
    % 更新权重,逐步放大
    
        
    % 更新consistency
    unaryListConsistency{lastIdx} = unaryListConsistency{curIdx};
    unaryListConsistency{curIdx} = cal_single_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode);%0.001秒的运行时间
%     debugunaryListConsistency(:,iter) = unaryListConsistency{curIdx};
    pairListConsistency{lastIdx} = pairListConsistency{curIdx};
    pairListConsistency{curIdx} = cal_pair_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode);
   
%     constMean(iter) = mean(pairListConsistency{curIdx}(logical(triu(ones(size(pairListConsistency{curIdx})),1))));
%     disp(['itr = ',num2str(iter),' c=',num2str(constMean(iter)),'err=',num2str(err(iter))]);
    if err(iter)==0
%         disp(['itr = ',num2str(iter),' c=',num2str(constMean(iter)),', err=',num2str(err(iter))]);
%         if constMean(iter) < 1%收敛到不动点，且该不动点非full consistent
%             disp(err);
%             disp(constMean);
%             dlmwrite('C:\pamilog1029.txt',[constMean(1:iter),err(1:iter)],'delimiter','\t','precision', 6,'-append');
%         end
        break;
    end
%     pairListConstDiff =  pairListConsistency{curIdx} - pairListConsistency{lastIdx};
%     a = length(find(pairListConstDiff<0));
%     b = length(find(pairListConstDiff>0));
%     if a >0, disp(['itr = ',num2str(iter),' exists dec samples ',num2str(a),'/',num2str(length(pairListConstDiff(:))),' inc samples ',num2str(b),'/',num2str(length(pairListConstDiff(:)))]);end
    % compute each pair of graphs for accuracy, score, consistency
    if calPairCorrFlag>0||calAnchCorrFlag>0
        acc = cal_pair_graph_accuracy(matchMat,affinity.Xgt,target.config.nOutlier,nodeCnt,graphCnt);
        scr = cal_pair_graph_score(matchMat,affinity.Xgt,nodeCnt,graphCnt);%要raw
%         scr = scr/scrDenom;
        cst = cal_pair_graph_consistency(matchMat,nodeCnt,graphCnt,massOutlierMode);
        acc = triu(acc,1)-tril(ones(graphCnt,graphCnt),0);scr = triu(scr,1)-tril(ones(graphCnt,graphCnt),0);cst = triu(cst,1)-tril(ones(graphCnt,graphCnt),0);
        acc((acc<0))=[];scr((scr<0))=[];cst((cst<0))=[];
        acc = acc(:);scr = scr(:);cst = cst(:);
        if calPairCorrFlag>0
            tmp = corrcoef(acc,scr);corr(iter:end,1)=tmp(2);
            tmp = corrcoef(acc,cst);corr(iter:end,2)=tmp(2);
            corr(iter:end,3)=mean(scr);
            corr(iter:end,4)=mean(acc);
            corr(iter:end,5)=mean(cst);
        end
        if calAnchCorrFlag>0%corr=[scr-acc,cst-acc,gcst-acc,pcst-acc,scr,acc,cst]
            corr(iter:end,1)=anchCorrSum(1)/anchCorrCnt;%scr-acc corr
            corr(iter:end,2)=anchCorrSum(2)/anchCorrCnt;%gcst-acc corr
            corr(iter:end,3)=anchCorrSum(3)/anchCorrCnt;%pcst-acc corr
            corr(iter:end,4)=anchCorrSum(4)/anchCorrCnt;%ecst-acc corr
            corr(iter:end,5)=anchCorrSum(5)/anchCorrCnt;%gquality-acc corr
            corr(iter:end,6)=anchCorrSum(6)/anchCorrCnt;%pquality-acc corr
            corr(iter:end,7)=anchCorrSum(7)/anchCorrCnt;%equality-acc corr
            corr(iter:end,8)=anchCorrSum(8)/anchCorrCnt;%gcst-cst err
            corr(iter:end,9)=anchCorrSum(9)/anchCorrCnt;%pcst-cst err
            corr(iter:end,10)=mean(scr);
            corr(iter:end,11)=mean(acc);
            corr(iter:end,12)=mean(cst);
        end
    end
end%for iter
% disp(['itr = ',num2str(iter)]);
% if iter == iterMax && err(iter)>0%到了最大迭代次数还未收敛
%     errOld = sum(abs(matchMatOldOld(:)-matchMat(:)));
%     disp(['itr = ',num2str(iter),' c=',num2str(constMean(iter)),', err=',num2str(err(iter))]);
%     dlmwrite('C:\pamilog1029.txt',errOld);
%     dlmwrite('C:\pamilog1029.txt',[constMean(1:iter),err(1:iter)],'delimiter','\t','precision', 6,'-append');
% end

nodeCon = nodeConsistencyMask;
P = matchMat;
%anchCorrMat = [scores,exaCon,trueAcc,graphCon,pairCon];
function [P anchCorr]= find1stOrderPathByScoreGraphPairConsistency(xview,yview,graphCon,pairListConsistency,...
    nodeCnt,graphCnt,constWeight,scrDenom,metricType,massOutlierMode,boostSampleRate,calAnchCorrFlag)
global affinity nodeConsistencyMask
global matchMat %permMat
    if calAnchCorrFlag>0
        trueAcc = zeros(graphCnt,1);
    end
    pairCon = zeros(graphCnt,1);
    xscope = (xview-1)*nodeCnt+1:xview*nodeCnt;
    yscope = (yview-1)*nodeCnt+1:yview*nodeCnt;
%     p0 = mat2vec(matchMat(xscope,yscope));
    Y = zeros(nodeCnt*nodeCnt,graphCnt);
%     permY = zeros(nodeCnt,graphCnt);
    for anchor=1:graphCnt
        ascope = (anchor-1)*nodeCnt+1:anchor*nodeCnt;
        P1 = matchMat(xscope,ascope);P2 = matchMat(ascope,yscope);
        Y(:,anchor) = mat2vec(P1*P2);%Xxy=Xxa*Xay
        if strcmp(metricType,'pair')||calAnchCorrFlag>0%pair也要动态计算,但代价较低，注意pairCon不能先uniqu再算，因为严格跟中间的k相关
            pairCon(anchor) = sqrt(pairListConsistency(xview,anchor)*pairListConsistency(anchor,yview));
        end
    end
    [~,m1,n1] = unique(Y','rows','first');%b1 = A(m1) and A = b1(n1):
    uniLen = length(m1);%totLen = length(n1);disp([num2str(totLen),' ',num2str(uniLen)]);
    uniCon = zeros(uniLen,1);uniScr = zeros(uniLen,1);
    if calAnchCorrFlag>0,uniAcc = zeros(uniLen,1);end
    % 只计算unique的变量
    for i=1:uniLen
        a = m1(i);
        p = sparse(Y(:,a));
        if massOutlierMode,b=repmat(nodeConsistencyMask(:,a)',nodeCnt,1);p = p.*b(:);end
        if ~strcmp(metricType,'cstcy')
            uniScr(i) = p'*affinity.K{xview,yview}*p/scrDenom;%score是必须的,cons是可选的
        end
%         tic
        %uniCon是算法2中的exact consistency
        if strcmp(metricType,'exact')||strcmp(metricType,'cstcy')||calAnchCorrFlag>0,uniCon(i) = cal_single_pair_consistency(matchMat,vec2mat(Y(:,a),nodeCnt,nodeCnt),xview,yview,nodeCnt,graphCnt,massOutlierMode);end
%         toc
%         tic%发现fast的版本反而比原始版本慢一些
%         if strcmp(metricType,'exact')||strcmp(metricType,'cstcy')||calAnchCorrFlag>0
%             [r,~,~] = find(vec2mat(Y(:,a),nodeCnt,nodeCnt));
%             uniCon(i) = fast_cal_single_pair_consistency(permMat((xview-1)*nodeCnt+1:xview*nodeCnt,:),permMat(:,yview),r,xview,yview,nodeCnt,graphCnt);
%         end
%         toc
        if calAnchCorrFlag>0,uniAcc(i) = cal_acc(vec2mat(p,nodeCnt,nodeCnt),target.config.nOutlier,affinity.GT(xscope,yscope));end
    end
    %pairCon graphCon没用unique 所以不用恢复
    scores = uniScr(n1);%从unique里恢复全部score
    if strcmp(metricType,'exact')||strcmp(metricType,'cstcy')||calAnchCorrFlag>0
        exaCon = uniCon(n1);%从unique里恢复全部exaCon
    end
    if calAnchCorrFlag>0
        trueAcc = uniAcc(n1);
        graphCstQuality = (1-constWeight)*scores + constWeight*graphCon;
        pairCstQuality = (1-constWeight)*scores + constWeight*pairCon;
        exaCstQuality = (1-constWeight)*scores + constWeight*exaCon;
    end
    
    switch metricType
        case 'unary'% CAO-UC
            quality = (1-constWeight)*scores + constWeight*graphCon;
        case 'pair'% CAO-PC
            quality = (1-constWeight)*scores + constWeight*pairCon;
        case 'exact'% CAO-C
            quality = (1-constWeight)*scores + constWeight*exaCon;
        case 'score'% CAO
            quality = scores;
        case 'cstcy'% used in Fig.1 in PAMI paper for CAO^{cst}
            quality = exaCon;
    end
    %对样本随机采样
%     quality(random('uniform',0,1,[graphCnt 1])>boostSampleRate) = 0;
    quality2 = quality; 
    quality = zeros(size(quality));
    if boostSampleRate>0 && boostSampleRate<=1
        quality(1:ceil(length(quality)*boostSampleRate)) = quality2(1:ceil(length(quality)*boostSampleRate));
    else
        quality = quality2;
    end
% % % % % % %     
    [~, idx] = max(quality);
    P = vec2mat(Y(:,idx(1)),nodeCnt,nodeCnt);
    anchCorr = zeros(10,1);
    if calAnchCorrFlag>0
        tmp = corrcoef(scores,trueAcc);anchCorr(1)=tmp(1,2);
        tmp = corrcoef(graphCon,trueAcc);anchCorr(2)=tmp(1,2);
        tmp = corrcoef(pairCon,trueAcc);anchCorr(3)=tmp(1,2);
        tmp = corrcoef(exaCon,trueAcc);anchCorr(4)=tmp(1,2);
        tmp = corrcoef(graphCstQuality,trueAcc);anchCorr(5)=tmp(1,2);
        tmp = corrcoef(pairCstQuality,trueAcc);anchCorr(6)=tmp(1,2);
        tmp = corrcoef(exaCstQuality,trueAcc);anchCorr(7)=tmp(1,2);
        anchCorr(8)=mean(abs(exaCon-graphCon));
        anchCorr(9)=mean(abs(exaCon-pairCon));
        anchCorr(10)=mean(scores);
    end