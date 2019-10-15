function test_IMGM
init_path;
set_alg;
set_conf_img;

cls = 'Duck';
listOfModelDataFile = dir(fullfile(conf.dataDir, 'feature', cls, '*.mat'));
gtDir = fullfile(conf.dataDir, 'ground_truth', cls);
for i = 1:length(listOfModelDataFile)
    path = fullfile(listOfModelDataFile(i).folder, listOfModelDataFile(i).name);
    cdatai = load(path);
    viewInfo(i) = cdatai.cdata.view;
end
tic
[ matchScore, hyperGraph, clusterMember, interClusterGraph, assignMat] = IMGM(alg(1), 110, viewInfo);
t = toc;
fprintf('time consumed: %.3fsec\n', t);
% Evaluation

nGraph = size(hyperGraph, 1);
np = size(viewInfo(1).desc, 2);
stk = zeros(1, nGraph, 'uint32');
visited = zeros(1, nGraph, 'logical');
top = 1;
stk(top) = 1;
acc = 0;

while top >= 1
    visited(stk(top)) = true;
    adjSet = hyperGraph(stk(top), :) & (~visited);
    adj = find(adjSet, 1);
    if isempty(adj)
        top = top - 1;
    else
        % find X(stk(top) -> adj)
        if stk(top) < adj
            name1 = viewInfo(stk(top)).fileName(end-7:end-4);
            name2 = viewInfo(adj).fileName(end-7:end-4);
            gtFileName = sprintf("%s_%s_%s.mat", cls, name1, name2);
            gt = load(fullfile(gtDir, gtFileName));
            gtX = gt.groundTruth.assign;
        else
            name1 = viewInfo(adj).fileName(end-7:end-4);
            name2 = viewInfo(stk(top)).fileName(end-7:end-4);
            gtFileName = sprintf("%s_%s_%s.mat",cls, name1, name2);
            gt = load(fullfile(gtDir, gtFileName));
            gtX = gt.groundTruth.assign';
        end
        X = assignMat{adj, stk(top)};
        % calculate accuracy
        acc = acc + sum(X.*gtX(:));
        
        % update stack
        top = top + 1;
        stk(top) = adj;
    end
end
acc = acc / (np*(nGraph-1));

fprintf("total accuracy = %.3f\n", acc);
end