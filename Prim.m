function hyperGraph = Prim(affScore)
    graphCnt = size(affScore, 1);
    graphCnt = int32(graphCnt);
    hyperGraph = zeros(graphCnt, graphCnt, 'logical');
    
    searched = zeros(1, graphCnt, 'logical');
    searched(1) = true;
    
    [~, I] = sort(affScore(:), 'descend');
    I = int32(I);
    
    a = mod(I - 1, graphCnt) + 1;
    b = idivide(I - 1, graphCnt) + 1;
    
    for i = 1:length(I)
        if size(unique(searched), 2) == 1
            break;
        end
        
        if searched(a(i)) && ~searched(b(i))
            searched(b(i)) = true;
            hyperGraph(a(i), b(i)) = 1;
            hyperGraph(b(i), a(i)) = 1;
        end
    end
    
%     if searched(a) && ~searched(b)
%         searched(b) = true;
%         hyperGraph(a, b) = 1;
%         hyperGraph(b, a) = 1;
%     end
    
%     I = int32(I);
%     
%     for ii = I
%         a = mod(ii-1, graphCnt) + 1;
%         b = idivide(ii-1, graphCnt) + 1;
%         if searched(a) && visited(b) && ~searched(b)
%             searched(b) = 1;
%             hyperGraph(a, b) = 1;
%             hyperGraph(b, a) = 1;
%         end
%     end
end