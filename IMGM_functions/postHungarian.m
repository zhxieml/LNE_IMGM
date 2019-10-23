function Xd = postHungarian( E12, X)

[n1 n2] = size(E12);
Xd = discretisationMatching_hungarian(reshape(X,n1,n2),E12);