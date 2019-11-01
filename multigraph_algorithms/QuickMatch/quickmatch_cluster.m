function clusters = quickmatch_cluster(data, n, param)
    % """

    % :param rho:
    % :param rho_edge:
    % :param data: [m*n, p]
    % :param n: number of points per image
    % :param kernel:
    % :param flag_low_memory:
    % :return:
    % """
    Scale = quickmatch_scale(data, n);
    dsq = euclidean_dist_matrix(data, data);
    tree_density = quickmatch_density(dsq, Scale, param);
    [tree_parents, tree_distance] = quickmatch_tree(tree_density, dsq, n);
    clusters = quickmatch_breaktree_merge(tree_parents, tree_distance, Scale, param.rho_edge);
    clusters(cellfun(@isempty,clusters))=[];
end

