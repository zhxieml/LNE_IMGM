function visited = dfs(graph, root, maxNum)
    np = length(graph);
    len = np+1;
    stack = ones(1, len, 'uint32');
    top = 1;
    stack(top) = root;

    visited = zeros(1, np, 'logical'); 
    n = 0;
    
    while top && n < maxNum
       visited(stack(top)) = true;
       n = n + 1;
       
       unVisited = graph(stack(top), :) & (~visited);
       top = top - 1;
       
       for ii = find(unVisited)
           top = top + 1;
           stack(top) = ii;
       end
    end
end