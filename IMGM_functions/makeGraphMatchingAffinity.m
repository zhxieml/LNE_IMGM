function affinity = makeGraphMatchingAffinity(parak,target,bfullconnect,bedgesparse)
if length(target)==1%虚拟随机实验
    Data = loadRandomData(target,bfullconnect,bedgesparse);%仿真实验没有点集，target存的是参数
else
    Data = loadTargetData(target,bfullconnect,bedgesparse);%真实实验有点集，都在target里存着
end
Sacle_2D = 0.15;
% [KP, KQ, G, H, EG] = calNodeEdgeAffinity(Data, Sacle_2D);
affinity = calNodeEdgeAffinity(Data, Sacle_2D);
% affinity.weight = [1 1 1];
if length(target)==1%仿真实验
    affinity.asgT{1}.X = eye(target.nInlier(min([parak,length(target.nInlier)]))+target.nOutlier(min([parak,length(target.nOutlier)])));
    affinity.nOutlier = target.nOutlier(min([parak,length(target.nOutlier)]));
else
    affinity.asgT{1}.X = eye(length(target{1}));
    affinity.nOutlier = 0;%真实实验目前没有outlier
end
% affinity.asgT{2} = affinity.asgT{1};
% affinity.asgT{3} = affinity.asgT{1};
affinity.adj{1} = Data{1}.adjMatrix;
affinity.adj{2} = Data{2}.adjMatrix;
% affinity.adj{3} = Data{3}.adjMatrix;
%假设outlier每个图数目一样

