function clusters = quickmatch_breaktree_merge(tree_parent, tree_distance, scales, rho_edge)
    npoint = size(tree_parent,1);
    clusters_indicator = 1:npoint;
    clusters = num2cell(1:npoint);
    match_dis = scales;

    [~, idx_sorted] = sort(tree_distance, 1);
    for i = idx_sorted
        p = tree_parent(i);
        c1 = clusters_indicator(i);
        c2 = clusters_indicator(p);
        if i ~= p && c1 ~= c2
            match_dis_c1_c2 = min(match_dis(c1), match_dis(c2));
            if tree_distance(i) <= rho_edge * match_dis_c1_c2
                clusters_indicator(clusters{c2}) = c1;
                clusters{c1} = [clusters{c1}, clusters{c2}];
                clusters{c2} = [];
                match_dis(c1) = match_dis_c1_c2;
            end
        end
    end
end