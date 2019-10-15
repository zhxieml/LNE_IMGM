function [ X score ] = SM( problem, AlgConfig )

if AlgConfig.bNormalize
    M = problem.Mnorm;
else
    M = problem.M;
end
%nSize = length(M);
%E12 = problem.E12;
% n1 = size(E12,1); % size of graph 1
% n2 = size(E12,2); % size of graph 2

%% Get eigenVector with largest eigenValue
options.disp = 0;
M = sparse(M);

%v = ones(nSize,1)/nSize;
%[eigenValue eigenVector] = powerit(M, v);
[eigenVectorMat eigenValue] = eigs(M, 1, 'lm', options);
eigenVector = eigenVectorMat(:,1);
eigenVector = abs(eigenVector);

% if 1
%     % discretization
%     Xt = zeros(size(E12));
%     Xt(find(E12)) = eigenVector(:);
%     Xd=discretisationMatching_hungarian(Xt, E12);
%     %Xd=computeDiscreteCorrespondancesGreedy(Xt, E12);
%     X = Xd(find(E12));
% else
%     [L12(:,1) L12(:,2)] = find(E12);
%     %[L12(:,2) L12(:,1)] = find(E12');    % transfrom E12 to list
%     L12(1:15,:)
%     X=greedyMapping(eigenVector, L12);
%     
% end

% Discretize using the overlap matrix 
X = greedyMapping(eigenVector, problem.overlapMatrix);
score = eigenVector;

if isa(AlgConfig.postProcess, 'function_handle')
    X = feval(AlgConfig.postProcess, problem, X);
end
