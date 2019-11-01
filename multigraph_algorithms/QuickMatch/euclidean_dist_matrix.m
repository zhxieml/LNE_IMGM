function dsq = euclidean_dist_matrix(a, b)
    % """
	% compute Euclidean Distance Matrix
	% :param a: [n1, p]
	% :param b: [n2, p]
	% :return: dsq [n1, n2]
	% """
	n1 = size(a, 1);
	n2 = size(b, 1);
	a_square = sum(a.^2, 2);
	b_square = sum(b.^2, 2);
	dsq = a_square*ones(1, n2) + ones(n1, 1)*b_square' - 2*a*b';
	dsq =  sqrt(dsq);
end