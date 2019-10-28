function hyperGraph = Prim(affScore)
    graphCnt = size(affScore, 1);
    graphCnt = int32(graphCnt);
    hyperGraph = zeros(graphCnt, graphCnt, 'logical');
    searched = zeros(1, graphCnt, 'logical');
    [~, I] = sort(affScore(:), 'descend');
    I = int32(I);
    
    a = mod(I - 1, graphCnt) + 1;
    b = idivide(I - 1, graphCnt) + 1;
    
    searched(a(1)) = true;
    for i = 1:length(I)
        if all(searched)
            break;
        end
        if searched(a(i)) && ~searched(b(i))
            searched(b(i)) = true;
            hyperGraph(a(i), b(i)) = 1;
            hyperGraph(b(i), a(i)) = 1;
        end
    end
end