%% Implementaion of PageRank matching algortihm
function [ X score ] = RWPRM( problem )
% use asymmetrically normalized transition matrix
M = problem.M;
[list(:,1) list(:,2)] = find(problem.E12);
answer = zeros(size(M,1),1);

iterMax = 10;
c = 0.7;
bUnbiased = 1;

%% augmented transition matrix
options.disp = 0;
M = sparse(M);
nSize = size(M,1);

Ma = M'; % transpose for readability

d = sum(Ma, 2); % degree : row sum
maxD = max(d);
%Ma(:,nSize+1) = maxD - d;
Ma = Ma ./ maxD;

Ma = Ma * c;
%sum(Ma,1)
%pause;

Ma(:,nSize+1) = [ d*(1-c)/maxD ];
%sum(Brw,1)
%size(Brw)
%pause;

Ma = Ma'; % transpose for readability

seeds = ones(nSize,1)./nSize;

% initialize answer
answer = zeros(nSize,1);
prev_score = ones(nSize,1);

nDiff = 1;
iter_i = 0;
%show_script_for_RWPRM;

while nDiff > 0 && iter_i < iterMax
    %fprintf('.');
    iter_i = iter_i + 1;
    
    seeds = zeros(nSize,1);
    seeds(find(answer)) = prev_score(find(answer));
    if bUnbiased
        seeds = Ma(1:nSize,1:nSize) * seeds;  seeds = seeds(1:nSize);
    end
    seedNorm = sum(seeds);
    if seedNorm>0 
        seeds = seeds./sum(seeds);
    end
    Ma(:,nSize+1) = [ seeds; 0 ];
    [eigenVectorMat eigenValue] = eigs(Ma, 1, 'lm', options);
%     eigenVectorMat
    %eigenValue
    eigenVector = eigenVectorMat(:,1);
    eigenVector = abs(eigenVector);
    %pause;
%     pageRank = ( (speye(nSize) - c*Brw) \ seeds ).* (1-c);
%     if bUnbiased
%         pageRank = pageRank - seeds .* (1-c); % bias subtraction
%         pageRank = pageRank./sum(pageRank);
%     end
    cur_score = eigenVector(1:nSize);
    
    prev_answer = answer;
    prev_score = cur_score;
    %% Get index from pageRank
    answer = zeros(size(M,1),1);
    while sum(cur_score) > 0 %&& length(answer) < 2
        [maxValue maxIdx] = max(cur_score);
        if maxValue == 0
            break;
        end
        answer(maxIdx) = 1;
        for i = 1:length(list)
            if list(maxIdx,1) == list(i,1) || list(maxIdx,2) == list(i,2) % one-to-one
                cur_score(i) = 0;
            end
        end
    end
    
    nDiff = sum(xor(prev_answer, answer));
    
    %show_script_for_RWPRM;
end
%fprintf('\n');

score = eigenVector;

X = reshape(answer, size(problem.E12,1), size(problem.E12,2));