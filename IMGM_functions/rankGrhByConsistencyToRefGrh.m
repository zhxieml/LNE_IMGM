function constGrhList = rankGrhByConsistencyToRefGrh(X,refGraph,nodeCnt,graphCnt)
otherGrhList = [1:refGraph-1,refGraph+1:graphCnt];
errX = zeros(graphCnt-1,1);
constX = zeros(graphCnt-1,1);
rscope = (refGraph-1)*nodeCnt+1:refGraph*nodeCnt;
for gk = 1:graphCnt-1
    x = otherGrhList(gk);
%     fprintf('%d=x,%d=r\n',x,refGraph);
    xscope = (x-1)*nodeCnt+1:x*nodeCnt;
    Xxr = X(xscope,rscope);
    for k=1:graphCnt
        if k==refGraph || k==x, continue;end
        kscope = (k-1)*nodeCnt+1:k*nodeCnt;
        % Xxr=Xxk*Xkr
        Xxkr = X(xscope,kscope)*X(kscope,rscope);
        errX(gk) = errX(gk) + sum(abs(Xxr(:) - Xxkr(:)))/(2*nodeCnt);%除以2是因为排列阵相减误差会被放大2倍
    end
    constX(gk) = 1-sum(errX(gk))/graphCnt;
end

% [B,IX]=sort(constX,'descend');
[B,IX]=sort(constX,'ascend');
constGrhList = otherGrhList(IX);
