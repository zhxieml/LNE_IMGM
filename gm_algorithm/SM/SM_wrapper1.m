function [ X score ] = SM_wrapper1( problem, AlgConfig )

if AlgConfig.bNormalize
    M = problem.Mnorm;
else
    M = problem.affinityMatrix;
end
%conf = problem.conflictMatrix1 | problem.conflictMatrix2;
%M = M.*(~conf);

%[ X score ] = SM( M, problem.group1, problem.group2);
[ score ] = SM( M );

X = greedyMapping(abs(score), problem.group1, problem.group2);

%if isa(AlgConfig.postProcess, 'function_handle')
%    X = feval(AlgConfig.postProcess, problem, X);
%end

