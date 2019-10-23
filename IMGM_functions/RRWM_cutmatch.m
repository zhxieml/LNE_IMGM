function [ X ] = RRWM_cutmatch( M, nP1, nP2, algpar)
%% Default Parameters
param = struct( ...
    'c', 0.2, ...                   % prob. for walk or reweighted jump?
    'amp_max', 30, ...              % maximum value for amplification procedure
    'thresConvergence', 1e-25, ...  % convergence threshold for random walks
    'tolC', 1e-3 ...                % convergence threshold for the Sinkhorn method
    );
E12 = ones(nP1,nP2);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);

%% parameter structure -> parameter value
strField = fieldnames(param);
for i = 1:length(strField), eval([strField{i} '=param.' strField{i} ';']); end
% % % % % % % % % % % % % % % % % % % %
% get groups for bistochastic normalization
[idx1 ID1] = make_groups(group1);
[idx2 ID2] = make_groups(group2);
%
if ID1(end) < ID2(end)
    [idx1 ID1 idx2 ID2 dumVal dumSize] = make_groups_slack(idx1, ID1, idx2, ID2);
    dumDim = 1;
elseif ID1(end) > ID2(end)
    [idx2 ID2 idx1 ID1 dumVal dumSize] = make_groups_slack(idx2, ID2, idx1, ID1);
    dumDim = 2;
else
    dumDim = 0; dumVal = 0; dumSize = 0;
end
idx1 = idx1-1; idx2 = idx2-1;

% note that this matrix is column-wise stochastic
d = sum(M, 1); % degree : col sum
maxD = max(d);
Mo = M ./ maxD; % nomalize by the max degree

% initialize answer
nMatch = length(M);
prev_score = ones(nMatch,1)/nMatch; % buffer for the previous score
prev_score2 = prev_score;         % buffer for the two iteration ahead
prev_assign = ones(nMatch,1)/nMatch; % buffer for Sinkhorn result of prev_score

bCont = 1;  iter_i = 0;
mode = 2;

%% start main iteration
while bCont && iter_i < algpar.iterMax1
    
    iter_i = iter_i + 1;
    
    %% random walking with reweighted jumps
    cur_score = Mo * ( c*prev_score + (1-c)*prev_assign );
    
    sumCurScore = sum(cur_score); % normalization of sum 1
    if sumCurScore>0, cur_score = cur_score./sumCurScore; end
    
    %% update reweighted jumps
    cur_assign = cur_score;
    % attenuate small values and amplify large values
    amp_value = amp_max/ max(cur_assign);  % compute amplification factor
    cur_assign = exp( amp_value*cur_assign );  % amplification
    
    % Sinkhorn method of iterative bistocastic normalizations
    X_slack = [cur_assign; dumVal*ones(dumSize,1)];
    if mode == 1
        X_slack = mexBistocNormalize_match_slack(X_slack, int32(idx1), int32(ID1), int32(idx2), int32(ID2), tolC, dumDim, dumVal, int32(1000));
    elseif mode == 2
        X_slack = reshape(X_slack,[nP1,nP2]);
        X_slack = BregmanBiStochZeroSlack(X_slack,100);
        X_slack = X_slack(:);
    end
    cur_assign = X_slack(1:nMatch);
    
    sumCurAssign = sum(cur_assign); % normalization of sum 1
    if sumCurAssign>0, cur_assign = cur_assign./sumCurAssign; end
    
    %% Check the convergence of random walks
    if 1
        diff1 = sum((cur_score-prev_score).^2);
        diff2 = sum((cur_score-prev_score2).^2); % to prevent oscillations
        diff_min = min(diff1, diff2);
        if diff_min < thresConvergence
            bCont = 0;
        end
    else
        normed_cur_score = cur_score/norm(cur_score);
        if norm(M*normed_cur_score - la*normed_cur_score,1) < thresConvergence2
            bCont = 0;
        end
        la = normed_cur_score'*M*normed_cur_score;
    end
    
    prev_score2 = prev_score;
    prev_score = cur_score;
    prev_assign = cur_assign;
    
end % end of main iteration

X = cur_score;
X = reshape(X,[nP1,nP2]);
% X = BregmanBiStochZero(X,100);
X = X(:);
end