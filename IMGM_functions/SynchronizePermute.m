function XX = SynchronizePermute(X,nodeCnt,graphCnt,mode)
if nargin == 3, mode='auto';end;
switch mode
%     case 'msts'
%         scoreAdj = generateSuperGraphScoreAdj(X,nodeCnt,graphCnt);
%         XX = minSpanTreeMatch(scoreAdj,X,nodeCnt,graphCnt);
%     case 'mstc'
%         cstAdj = cal_pair_graph_consistency(X,nodeCnt,graphCnt);
%         XX = minSpanTreeMatch(cstAdj,X,nodeCnt,graphCnt);
    case 'sync'
        if nodeCnt<=graphCnt
        XX = PMGM(X,nodeCnt,graphCnt);
        else
            cstAdj = cal_pair_graph_consistency(X,nodeCnt,graphCnt);
            XX = minSpanTreeMatch(cstAdj,X,nodeCnt,graphCnt);
        end
%     case 'rdm'
%         ref = 1;
%         XX = ConsistentReferenceMultiMatch(ref,nodeCnt,graphCnt,'',X);
    case 'auto'
        cstAdj = cal_pair_graph_consistency(X,nodeCnt,graphCnt);
        cst = mean(cstAdj);
        if cst>0.3
            if nodeCnt<=graphCnt
                XX = PMGM(X,nodeCnt,graphCnt);
            else
                cstAdj = cal_pair_graph_consistency(X,nodeCnt,graphCnt);
                XX = minSpanTreeMatch(cstAdj,X,nodeCnt,graphCnt);
            end
        else
            scoreAdj = generateSuperGraphScoreAdj(X,nodeCnt,graphCnt);
            XX = minSpanTreeMatch(scoreAdj,X,nodeCnt,graphCnt);
        end
end
    