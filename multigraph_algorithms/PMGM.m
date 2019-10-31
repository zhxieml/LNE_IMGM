% permutation sync for graph matching
% see more details in the NIP13 paper:
% Pachauri et al. Solving the multi-way matching problem by permutation synchronization, NIPS 2013
function [X] = PMGM(rawMat,nodeCnt,graphCnt)
n = nodeCnt;
% N = graphCnt;
N = max([nodeCnt,graphCnt]);% for PMGM, graphCnt must be larger than nodeCnt
cutoffeig_num = n;

[u,~,~] = svd(rawMat);%[U,S,V] = SVD(X) X = U*S*V'
cu = u(:,1:cutoffeig_num);%use the first cutoffeig_num eigen vectors

B={};
for i = 1:n:N*n%divide each left-singular vectors into N pieces£ºn*N x cutoffeig_num -> n x cutoffeig_num
    B{(i+(n-1))/n} = cu(i:(i+n)-1,:);%n x cutoffeig_num
end
% Recovered permutation
recPM=cell(N,1);
recDiscMat = zeros(nodeCnt*graphCnt,nodeCnt*graphCnt);
I = eye(n);
for i = 1:N
    tmp =  B{i}*B{1}';%use B{1} as the reference, see more details in the algrithmic chart of the paper
    [tmpU,~,tmpV]=svd(tmp);
    recPM{i} = tmpU*tmpV';%closest orthogonal matrix to tmp
    % discretize
    [c,~] = munkres(-recPM{i}');
    recPM{i} = I(:,c);
% % % % % % % % % % % % % % % % % % % 
end
for i = 1:N
    for j = i+1:N
        recDiscMat((i-1)*n+1:i*n,(j-1)*n+1:j*n) = recPM{i}*recPM{j}';       
    end
end
X = recDiscMat + recDiscMat' + eye(nodeCnt*graphCnt);