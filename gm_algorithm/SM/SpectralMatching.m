%% Implementaion of Leordeanu Marius' spectral matching algortihm
% by Jungmin Lee in CVL
function X = SpectralMatching( M, E12 )

n1 = size(E12,1); % size of graph 1
n2 = size(E12,2); % size of graph 2

[L12(:,1) L12(:,2)] = find(E12);

%% Get eigenVector with largest eigenValue
options.disp = 0;
M = sparse(M);
[eigenVectorMat eigenValue] = eigs(M, 1, 'lm', options);
eigenVector = eigenVectorMat(:,1);
eigenVector = abs(eigenVector);

%% Get index from principle eigenVector
% X = zeros(length(M),1);
% while sum(eigenVector) > 0
%     [maxValue maxIdx] = max(eigenVector);
%     if maxValue == 0
%         break;
%     end
%     X(maxIdx) = 1;
%     eigenVector(find(L12(:,1) == L12(maxIdx,1))) = 0;
%     eigenVector(find(L12(:,2) == L12(maxIdx,2))) = 0;
% end
% 
% X = reshape(X, size(E12,1), size(E12,2));

X = reshape(eigenVector, n1, n2);
Xd=discretisationMatching_hungarian(X, E12);
X = Xd;
