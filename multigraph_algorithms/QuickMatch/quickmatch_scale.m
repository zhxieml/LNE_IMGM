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
    for ii = 1:n:nPoint
        dsq_i = euclidean_dist_matrix(data(ii:ii+n), data(ii:ii+n));
        for k = 1:n
            min_d = 1e5;
            for y = 1:n
                if k == y
                    continue;
                end
                min_d = min(min_d, dsq_i(k, y));
            end
            Scale(ii+k) = min_d;
        end
    end

end