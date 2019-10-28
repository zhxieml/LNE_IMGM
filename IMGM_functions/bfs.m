function visited = bfs(graph, root, maxNum)
    np = length(graph);
    len = np+1;
    queue = ones(1, len, 'uint32');
    front = 1;
    tail = 2;
    queue(front) = root;
    visited = zeros(1, np, 'logical'); 
    n = 0;
    while front ~= tail && n < maxNum
        visited(queue(front)) = true;
        n = n+1;
        unVisited = graph(queue(front), :) & (~visited);
        front = mod(front, len)+1;
        if nnz(unVisited) ~= 0 && n < maxNum
            for ii = find(unVisited)
                queue(tail) = ii;
                tail = mod(tail, len)+1;
            end
        end
    end
end