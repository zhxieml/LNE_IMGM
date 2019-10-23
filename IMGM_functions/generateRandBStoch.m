function P = generateRandBStoch(n, weight)
    P = zeros(n,n);
    I = eye(n);
    for i = 1:length(weight)
        randInd = randperm(n);
        randMat = I(randInd,:);
        P = P + weight(i)*randMat;
    end
end