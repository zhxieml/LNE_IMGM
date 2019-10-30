function tree_density = quickmatch_density(data, scales, param) 
    % rho=0.25;
    % kernel=None;
    % kernel_support_radius=None;
    % flag_low_memory=False;
    % chunk_size=None;
    % flag_recursive=False;
    rho = param.rho;
    kernel = param.kernel;
    kernel_support_radius = param.kernel_support_radius;
    flag_low_memory = param.flag_low_memory;
    chunk_size = param.chunk_size;
    flag_recursive = param.flag_recursive;

    Sigma = rho * scales;
    if (~flag_low_memory) || flag_recursive
        dsq = data;
        dsq = dsq / Sigma;
        if ~kernel
            p = dsq;
        else
            if kernel_support_radius
                p = zeros(size(dsq));
                flag = dsq < kernel_support_radius;
                p(flag) = kernel(-dsq(flag));
            else
                p = kernel(-dsq);
            end
        end
        tree_density = sum(p,1);
    else
        % here data is [p, n] raw data
        npoint = size(data, 1);
        tree_density = zeros(1, npoint);
        for ii = 1:chunk_size:npoint
            jj = min(ii+chunk_size, npoint);
            dsq = euclidean_dist_matrix(data(:, ii:jj), data);
            param.flag_recursive = true;
            tree_density(ii:jj) = quickmatch_density(dsq, scales, param);
        end
    end
end