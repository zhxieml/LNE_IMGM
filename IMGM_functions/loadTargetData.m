function Data = loadTargetData(target,bfullconnect,bedgesparse)

Data=cell(2,1);
if length(target)>1
for viewk=1:2
    Data{viewk}.nP = size(target{viewk}.point,1);%点的数目
    Data{viewk}.edge = zeros(Data{viewk}.nP,Data{viewk}.nP);%边
    Data{viewk}.point = target{viewk}.point;%点的坐标
end
for r = 1: Data{1}.nP
     for c = r+1: Data{1}.nP
        Data{1}.edge(r,c) = sqrt((target{1}.point(r,1)-target{1}.point(c,1))^2+(target{1}.point(r,2)-target{1}.point(c,2))^2);
        Data{2}.edge(r,c) = sqrt((target{2}.point(r,1)-target{2}.point(c,1))^2+(target{2}.point(r,2)-target{2}.point(c,2))^2);
%         Data{3}.edge(r,c) = sqrt((target{3}.point(r,1)-target{3}.point(c,1))^2+(target{3}.point(r,2)-target{3}.point(c,2))^2);
     end
end
Data{1}.edge = Data{1}.edge/max(Data{1}.edge(:));
Data{2}.edge = Data{2}.edge/max(Data{2}.edge(:));
% Data{3}.edge = Data{3}.edge/mx(Data{3}.edge(:));

Data{1}.edge = Data{1}.edge + Data{1}.edge';
Data{2}.edge = Data{2}.edge + Data{2}.edge';
% Data{3}.edge = Data{3}.edge + Data{3}.edge';

if bfullconnect
    maskDelaunay = ones(Data{1}.nP, Data{1}.nP);
else
% 在这里稀疏化edge矩阵
for viewk=1:2
Data{viewk}.tri = delaunay(Data{viewk}.point(:,1),Data{viewk}.point(:,2));
% Data.tri{viewk} = Data.tri{viewk}(1:sideDownSample:end,:);
Data{viewk}.adjMatrix = zeros( Data{viewk}.nP, Data{viewk}.nP);
triNum=size(Data{viewk}.tri,1);
for i=1:triNum
    Data{viewk}.adjMatrix(Data{viewk}.tri(i,1),Data{viewk}.tri(i,2))=1;
    Data{viewk}.adjMatrix(Data{viewk}.tri(i,2),Data{viewk}.tri(i,1))=1;
    Data{viewk}.adjMatrix(Data{viewk}.tri(i,2),Data{viewk}.tri(i,3))=1;
    Data{viewk}.adjMatrix(Data{viewk}.tri(i,3),Data{viewk}.tri(i,2))=1;
    Data{viewk}.adjMatrix(Data{viewk}.tri(i,1),Data{viewk}.tri(i,3))=1;
    Data{viewk}.adjMatrix(Data{viewk}.tri(i,3),Data{viewk}.tri(i,1))=1;
end
end
% % sparsify the affinity matrix,和IJCV的策略一样，取三个图delaunay的交集,feng 的fgm也用了delaunay
if bedgesparse==0
    maskDelaunay = Data{1}.adjMatrix+Data{2}.adjMatrix;
else
   maskDelaunay = Data{1}.adjMatrix.*Data{2}.adjMatrix; 
end
% maskDelaunay = maskDelaunayAB+maskDelaunayAC+maskDelaunayBC;%这里又用了下并集
% maskDelaunay = Data{1}.adjMatrix.*Data{2}.adjMatrix.*Data{3}.adjMatrix;
maskDelaunay= logical(maskDelaunay);
for viewk=1:2
Data{viewk}.edge(~maskDelaunay) = NaN;
end
end
%更新adjMatrix，取交集
Data{1}.adjMatrix = maskDelaunay;
Data{2}.adjMatrix = maskDelaunay;
% Data{3}.adjMatrix = maskDelaunay;   
% Data{1}.adjMatrix = logical(Data{1}.adjMatrix);
% Data{2}.adjMatrix = logical(Data{2}.adjMatrix);
% Data{3}.adjMatrix = logical(Data{3}.adjMatrix);   
Data{1}.nE = sum(sum(maskDelaunay));
Data{2}.nE = Data{1}.nE;
% Data{3}.nE = Data{1}.nE;
% Data{1}.nE = sum(sum(Data{1}.adjMatrix));
% Data{2}.nE = sum(sum(Data{2}.adjMatrix));
% Data{3}.nE = sum(sum(Data{3}.adjMatrix));
end