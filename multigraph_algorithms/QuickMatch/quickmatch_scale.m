function Scale = quickmatch_scale(data, n)
% """
% calculate the $d_{ik}$ for each point
% :param data: [m*n, p] data, m is number of images
%             for data[x][p], i = x // n, k = x % n
% :param n: number of points per image
% :return: a $d_{ik}$ array of shape[m*n,]
% """
    nPoint = size(data, 1);
    Scale = zeros(1, nPoint);
    for ii = 0:n:nPoint-1
        dsq_i = euclidean_dist_matrix(data(ii+1:ii+n, :), data(ii+1:ii+n, :));
        dsq_i = dsq_i + 1e5*eye(n);
        k = 1:n;
        Scale(ii+k) =  min(dsq_i, [], 1);
    end
end