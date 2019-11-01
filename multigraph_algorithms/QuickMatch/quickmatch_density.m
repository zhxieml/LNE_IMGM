function tree_density = quickmatch_density(data, scales, param) 
    rho = 0.25;
    kernel = [];
    kernel_support_radius = false;
    flag_low_memory = false;
    chunk_size = 0;
    flag_recursive = false;
    if isfield(param, 'rho')
        rho = param.rho;
    end
    if isfield(param, 'kernel')
        kernel = param.kernel;
    end
    if isfield(param, 'kernel_support_radius')
        kernel_support_radius = param.kernel_support_radius;
    end
    if isfield(param, 'flag_low_memory')
        flag_low_memory = param.flag_low_memory;
    end
    if isfield(param, 'chunk_size')
        chunk_size = param.chunk_size;
    end
    if isfield(param, 'flag_recursive')
        flag_recursive = param.flag_recursive;
    end

    Sigma = rho * scales;
    if (~flag_low_memory) || flag_recursive
        dsq = data;
        dsq = dsq / Sigma;
        if isempty(kernel)
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