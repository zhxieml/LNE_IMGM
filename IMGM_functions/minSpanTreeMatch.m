function [minSpanTreeMat treeAdjRes] = minSpanTreeMatch(scoreAdj,rawMat,nodeCnt,graphCnt)
global treeAdj visitedNodes seq
bPlot = 0;
I = eye(nodeCnt);
minSpanTreeMat = zeros(nodeCnt*graphCnt,nodeCnt*graphCnt);
treeAdj = buildMinSpanTree(scoreAdj);
% if bPlot
%     totalFigNum = 5;
%     subplot(1,totalFigNum,1);shM(scoreAdj);title('scoreAdj');
%     subplot(1,totalFigNum,2);shM(rawMat);title('rawMat');
% end
% for each pair of start and end points, compute the path on the span tree
for startp=1:graphCnt
    for endp = startp+1:graphCnt
        seq = [];
        visitedNodes = zeros(graphCnt,1);
        seq(1) = startp;
        visitedNodes(startp) = 1;
        findPathP2P(startp,endp,nodeCnt,graphCnt);%存在seq里面
        P = I;
        for s=1:length(seq)-1
            xview = seq(s);
            yview = seq(s+1);
            P = P*rawMat((xview-1)*nodeCnt+1:xview*nodeCnt,(yview-1)*nodeCnt+1:yview*nodeCnt);
        end
        minSpanTreeMat((startp-1)*nodeCnt+1:startp*nodeCnt,(endp-1)*nodeCnt+1:endp*nodeCnt) = P;
    end
end
minSpanTreeMat = minSpanTreeMat + minSpanTreeMat' + eye(nodeCnt*graphCnt);
if bPlot
subplot(1,totalFigNum,3);shM(minSpanTreeMat);title('minSpanTreeMat');
spanMat = eye(nodeCnt*graphCnt);
for i=1:nodeCnt
    for j =1:nodeCnt
        if treeAdj(i,j)
            spanMat((i-1)*graphCnt+1:i*graphCnt,(j-1)*graphCnt+1:j*graphCnt) = rawMat((i-1)*graphCnt+1:i*graphCnt,(j-1)*graphCnt+1:j*graphCnt);
        end
    end
end
subplot(1,totalFigNum,4);shM(spanMat);title('spanBaseMat');
subplot(1,totalFigNum,5);shM(treeAdj);title('treeAdj');
end
treeAdjRes = treeAdj;
% ok = 1;

function success = findPathP2P(currentp,endp,nodeCnt,graphCnt)%为节点r找相连节点c
global treeAdj visitedNodes seq %minSpanTreeMat
success = 0;
seqLen = length(seq);
for c=1:graphCnt
    if visitedNodes(c)>0, continue;end
    if treeAdj(currentp,c)>0
%         fprintf('%d ',c);
        seq(end+1) = c;
        visitedNodes(c) = 1;
        moveFlag = 1;
        if c == endp
            success = 1; 
%             fprintf(' success (%d) ',c);
            return;
        else
            success = success + findPathP2P(c,endp,nodeCnt,graphCnt);
        end
    end
end
if success == 0
    if seqLen>1
        seq = seq(1:seqLen-1);
    end
end