function [ X score ] = SM( affinityMatrix, group1, group2 )

%% Get eigenVector with largest eigenValue
options.disp = 0;
M = sparse(affinityMatrix);

%v = ones(nSize,1)/nSize;
%[eigenValue eigenVector] = powerit(M, v);
[eigenVectorMat eigenValue] = eigs(M, 1, 'lm', options);
eigenVector = eigenVectorMat(:,1);
eigenVector = abs(eigenVector);

% Discretize using the overlap matrix 
X = greedyMapping(eigenVector, group1, group2);
score = eigenVector;
