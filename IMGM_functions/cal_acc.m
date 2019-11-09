function acc = cal_acc(P,nOutlier,GT)
%P: Nx*Ny, GT: Nx*Ny
% nOutlier = 0;
[Nx Ny]=size(P);
if nargin<3
    GT = zeros(Nx,Ny);
    GT(:,1:Nx) = diag(ones(1,Nx));
end
if Ny==Nx%������outlier
    P = P(1:end-nOutlier,1:end-nOutlier);%ֻ��Nx �ڵ㣬�������?
    Nx = size(P,1);
    GT = GT(1:end-nOutlier,1:end-nOutlier);
end
acc = (Nx-sum(sum(abs(P-GT)))/2)/Nx;
% acc = sum(sum(abs(P-GT),2)==0)/Nx;%��ÿ�����?,ֻ��ÿ������Ԫ�ض���ͬ������Ϊƥ��׼ȷ����ֹ����������£�P(1:nInlier)��û��1�����?
% acc = sum(GT(:).*P(:))/Nx;%��ÿ�����?,ֻ��ÿ������Ԫ�ض���ͬ������Ϊƥ��׼ȷ����ֹ����������£�P(1:nInlier)��û��1�����?