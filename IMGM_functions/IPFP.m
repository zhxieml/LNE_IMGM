%% Integer Projected Fixed Point Method
% by Jungmin Lee in CVL
%
% input:  affinity matrix M
%         edge incidence matrix E12 ( nNode of graph1 by nNode of graph2 )
%
function X = IPFP( M, nP1, nP2, varargin )

E12 = ones(nP1,nP2);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);
%% Step 1: initialize
n1 = size(group1,2);
n2 = size(group2,2);
E12 = ones(n1,n2);
nSize = n1*n2;

for i = 1:n1
    tmp1 = find(group1(:,i));
    match(tmp1,1) = i;
end
for i = 1:n2
    tmp2 = find(group2(:,i));
    match(tmp2,2) = i;
end

for i = 1:length(match)
    full(i,1) = (match(i,2)-1)*n1+match(i,1);
end

%% get groups for bistochastic normalization
% [idx1 ID1] = make_groups(group1);
% [idx2 ID2] = make_groups(group2);
% 
% if ID1(end) < ID2(end)
%     [idx1 ID1 idx2 ID2 dumVal dumSize] = make_groups_slack(idx1, ID1, idx2, ID2);
%     dumDim = 1;
% elseif ID1(end) > ID2(end)
%     [idx2 ID2 idx1 ID1 dumVal dumSize] = make_groups_slack(idx2, ID2, idx1, ID1);
%     dumDim = 2;
% else
%     dumDim = 0; dumVal = 0; dumSize = 0;
% end
% idx1 = idx1-1; idx2 = idx2-1;

tolC = 1e-3;
nMatch = length(M);
X = ones(nMatch,1)/nMatch;
B = X;
X_answer = X;

S = X'*M*X;
k = 0;

%% Step 2
while 1
    B1 = M*X;
    B_temp = zeros(nSize,1);
    for i = 1:length(full)
        B_temp(full(i,1),1) = B1(i,1);
    end
    B_mat = postHungarian( E12, B_temp);
    B_mat_temp = reshape(B_mat, nSize, 1);
    for i = 1:length(full)
        B(i,1) = B_mat_temp(full(i,1),1);
    end
%     B1 = exp(1000*B);
%     X_slack = [B1; dumVal*ones(dumSize,1)];
%     X_slack = mexBistocNormalize_match_slack(X_slack, int32(idx1), int32(ID1), int32(idx2), int32(ID2), tolC, dumDim, dumVal, int32(1000));
%     B_soft = X_slack(1:nMatch);
%     B = greedyMapping(B_soft, group1, group2);
    
    sumB = sum(B); % normalization of sum 1
    if sumB>0, B = B./sumB; end
    
    C = X'*M*(B-X);
    D = (B-X)'*M*(B-X);
    
    %% Step 3
    if D >= 0
        X2 = B;
    else
        r = min(-C/D, 1);
        X2 = X+r*(B-X);
    end
    
    %% Step 4
    if B'*M*B >= S
        S = B'*M*B;
        X_answer = B;
    end
    
    %% Step 5
    if sum(X2 - X) == 0 || k > 100
        if k == 100
            disp('Exception! in IPFP!')
        end
        break
    end
    
    %% Step 6
    k = k+1;
    X = X2;
end

%% Return the solution
%X = reshape(X_answer, n1, n2);
X = X_answer;
from1to2 = sum(X,2);
from2to1 = sum(X,1);

if ~isempty(find(from1to2 > 1)) || ~isempty(find(from2to1 > 1))
	%X = reshape(B, n1, n2);
    X = B;
end

end


