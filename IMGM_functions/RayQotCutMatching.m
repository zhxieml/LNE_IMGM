% function RayQotCutMatching calculates the eigenvector corresponding to
% the largest eigenvalue in terms of Rayleigh quotient
% Input:
%     Weight: weight matrix encoding the edge weights
%     X: matching within the graph
%     lambda_1, lambda_2: refer to paper
%     Type: 1. Rayleigh quotient
%           2. normalized cut according to Shi (2002)
%           3. normalized according to Jordan & Weiss (2002)
function y = RayQotCutMatching(Weight, X, lambda_1, lambda_2, Type)
    tol = 1e-6;
    degs = -lambda_1*sum(Weight,2)+lambda_2*sum(X,2);
    
    D = sparse(1:size(Weight,1),1:size(Weight,2),degs);
    W = -lambda_1*Weight + lambda_2*X;
    L = D - W;
    
    switch Type
        case 2
            degs(degs==0)=eps;
            D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
        
            % calculate normalized Laplacian
            L = D * L;
        case 3
            % avoid dividing by zero
            degs(degs == 0) = eps;
            % calculate D^(-1/2)
            D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
        
            % calculate normalized Laplacian
            L = D * L * D;
    end
    
    diff   = eps;
    L = (L+L')/2;
    [U, V] = eigs(-L, 2, 'SA');
    V = diag(V);
    if abs(0-V(1))<tol
        y = U(:,2);
    else
        y = U(:,1);
    end
    a =1;   
end