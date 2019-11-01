function [tree_parent, tree_distance] = quickmatch_tree(density, dsq, n)
    % :param density:
    % :param dsq: distance matrix
    % :param n: number of points per image
    % :return:
    n = int32(n);
    nPoint = size(density, 1);
    tree_parent = zeros(1, nPoint,'int32');
    tree_distance = zeros(1, nPoint);

    for ii = 1:nPoint
        Parent = ii;
        for jj = 1:nPoint
            % parent in the tree is given by the closest point
            % with higher density from another image
            if density(jj) <= density(ii) || idivide(ii, n) == idivide(jj, n)
                % density is lower or from the same image
                continue;
            end
            if Parent == ii || dsq(ii, jj) < dsq(ii, Parent)
                Parent = jj;
            end
        end
        tree_parent(ii) = Parent;
        tree_distance(ii) = dsq(ii, Parent);
    end