function dataset = readCMUData(DataPath,category,nodeCnt,bLoadImg,bPermute,bSCF)
if bSCF
    load house.mat
    len = length(data);
else
    ptFiles = dir(fullfile([DataPath,'/point/'],[category,'*']));
    len = length(ptFiles);
end    
for viewk=1:len
        if bSCF
            dataset(viewk).point = cell2mat(data(viewk));
            dataset(viewk).scf = cell2mat(scf(viewk));
            dataset(viewk).scf = dataset(viewk).scf(1:nodeCnt,:);
        else
            dataset(viewk).point = load([DataPath,'/point/',category,num2str(viewk)]);
        end
        dataset(viewk).point = dataset(viewk).point(1:nodeCnt,:);
        
        if viewk>1 && bPermute>0
            switch bPermute
                case 2
                    dataset(viewk).gt = randperm(length(dataset(viewk).point));
                case 1
                    dataset(viewk).gt = nodeCnt:-1:1;
            end
            dataset(viewk).point = dataset(viewk).point(dataset(viewk).gt,:);
            if bSCF
                dataset(viewk).scf = dataset(viewk).scf(dataset(viewk).gt,:);
            end
        else
            dataset(viewk).gt = 1:length(dataset(viewk).point);
        end
        I = eye(nodeCnt);
        %下面求Pi1
        dataset(viewk).gtPerm = I(:,dataset(viewk).gt)';%%%%%%%%注意取转置才是Pi1

        if bLoadImg
            dataset(viewk).img = imread([DataPath,'/img/',category,num2str(viewk),'.png']);
        end

    end
end