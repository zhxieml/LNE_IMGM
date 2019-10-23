function treeAdj = buildMinSpanTree(scoreAdj)
% [treeAdj, Cost] =  UndirectedMaximumSpanningTree(max(scoreAdj(:))-scoreAdj);
[treeAdj, Cost] =  UndirectedMaximumSpanningTree(scoreAdj);
