%% Implementaion of PageRank matching algortihm
function [ X score ] = RWPRM_v5( problem, AlgConfig )
% use asymmetrically normalized transition matrix
M = problem.M;
E12= problem.E12;
nSize = length(M);
if nSize ~= nnz(E12) 
    error('E12 and M conflict in size!');
end
[L12(:,1) L12(:,2)] = find(E12);% transfrom E12 to list
%eliminate conflicting elements
for i = 1:nSize
    for j = 1:nSize
        if L12(i,1) == L12(j,1) ||...
                L12(i,2) == L12(j,2) % one-to-one
            M(i,j) = 0;
            M(j,i) = 0;
        end
    end
end

isDisplay = 0;

%c_i = [ 1.0 0.8*ones(1,50) ];
%c = 0.2;
beta = 30;
c = 0.2; % no reweighting
iterMax = 300;
tolC=1e-3;%1e-3;
%thresConvergence=1e-25;

%% augmented transition matrix
options.disp = 0;
M = sparse(M);

% note that this matrix is column-wise stochastic
d = sum(M, 1); % degree : col sum
maxD = max(d);

Mo = M; 
Mo = Mo ./ maxD;

% initialize answer
answer = zeros(nSize,1);
prev_score = ones(nSize,1)/nSize;
prev_score2 = prev_score;
prev_assign = ones(nSize,1)/nSize;
%prev_assign = zeros(nSize,1);

% for convergence check of power iteration
thresConvergence = length(prev_score)*norm(M,1)*eps;
la = prev_score'*M*prev_score;

nDiff = 1;
iter_i = 0;
%show_script_for_RWPRM;
Xt = zeros(size(E12));

if isDisplay
    figure(1);
end

while nDiff > 0 && iter_i < iterMax
    
    iter_i = iter_i + 1;
    
    if isDisplay
        subplot(231);
        imagesc(reshape(prev_score, n1, n2));
        title(sprintf('(%d) prev score',iter_i));
    end
    
    cur_score = Mo * ( c*prev_score + (1-c)*prev_assign );
    %cur_score = c*Mo*prev_score + (1-c)*prev_assign;
    
    sumCurScore = sum(cur_score); % normalization
    if sumCurScore>0, cur_score = cur_score./sumCurScore; end
    
    if isDisplay
        subplot(232);
        imagesc(reshape(cur_score, n1, n2));
        title(sprintf('(%d) updated score',iter_i));
    end
    %cur_assign = cur_score;
    %cur_assign = exp( (5*2^(iter_i-1)*n1*n2)*cur_score );

    cur_assign = cur_score;
    %cur_assign = cur_assign - min(cur_assign);
    amp_value = beta/ max(cur_assign);%log(10000) 
    %amp_value = 20 * nSize;%log(10000) 
    cur_assign = exp( amp_value*cur_assign );      

    if isDisplay
        subplot(233);
        imagesc(reshape(cur_assign, n1, n2));
        title(sprintf('(%d) exponentiation',iter_i));
    end
    
    Xt(find(E12)) = cur_assign(:);
    [X,Xslack]=bistocNormalize_slack(Xt,tolC);
    cur_assign = X(find(E12));
    
    sumCurAssign = sum(cur_assign); % normalization
    if sumCurAssign>0, cur_assign = cur_assign./sumCurAssign; end
    
    %X
    if isDisplay
        subplot(234);
        imagesc(reshape(cur_assign, n1, n2));
        %bar3(reshape(cur_score, n1, n2),'detachted');
        title(sprintf('(%d) softassign',iter_i));
    end
    
    if isDisplay
        %seeds = zeros(nSize,1);
        %seeds(find(answer)) = cur_assign(find(answer));
        subplot(236);
        %bar3(problem.T12,'detachted');
        imagesc(problem.T12);
        title(sprintf('ground truth',iter_i));
        pause;
    end

    %nDiff = 1;%sum(xor(prev_answer, answer));        
    %prevent oscillations
%    diff1 = sum((cur_score-prev_score).^2);
%    diff2 = sum((cur_score-prev_score2).^2);
%    diff_min = min(diff1, diff2);
%    if diff_min < thresConvergence
%        nDiff = 0;
%    end
    normed_cur_score = cur_score/norm(cur_score);
    if norm(M*normed_cur_score - la*normed_cur_score,1) < thresConvergence
        nDiff = 0;
    end
    
    %prev_score2 = prev_score;
    la = normed_cur_score'*M*normed_cur_score;
    prev_score = cur_score;
    prev_assign = cur_assign;
    
end
%fprintf(' %d',iter_i);
score = cur_score;

%% discretize the solution    
% discretization
Xt(find(E12)) = cur_score(:);
%Xt(find(E12)) = cur_assign(:);
%     X2=computeXorthonormal(X2);
Xd=discretisationMatching_hungarian(Xt, E12);
%Xd=computeDiscreteCorrespondancesGreedy(Xt, E12);
X = Xd(find(E12));

%answer = X(:);
%X = reshape(answer, n1, n2);

end


function [X,Xslack]=bistocNormalize_slack(X,tolC)
[n1,n2]=size(X);
if n1~=n2
    Xslack=X;
    if n1>n2
        Xslack(:,n2+1:n1)=1;
    else
        Xslack(n1+1:n2,:)=1;
    end
    Xslack = bistocNormalize(Xslack,tolC,1000);
    X=Xslack(1:n1,1:n2);
else
    Xslack=X;
    X = bistocNormalize(X,tolC,1000);
end
end

function display_result(X,idxFig)
if nargin>1
    fig = figure(idxFig);
    set(fig, 'MenuBar', 'none');
    set(fig, 'Toolbar', 'none');
end

%imagesc(X);
bar3(X,'detachted');
drawnow;
%caxis([0,1]);
%pause(0.01);
%disp2('bounds(X(:))');
%colormap(gray);

end