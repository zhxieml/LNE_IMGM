nTest = 1;
baseGraphList = [1:15];

acc = sum(accResult,3)/nTest;
scr = sum(scrResult,3)/nTest;
con = sum(conResult,3)/nTest;
tim = sum(timResult,3)/nTest;
metricName = {'Accuracy','Affinity Score','Consistency','Time(s)'};
metricData = {acc,scr,con,tim};
algName = {'raw CAO','incremental CAO','IMGM-D','IMGM-R','incremental ConMatch','mSync'};
lineStyleList = {'-.','-','-.','-','--','-'};
markerList = {'o','s','d','*','p','+'};
colorList = {'b','y','k','g','r','m'};

for i = 1:4
    figure; hold on; title(metricName{i});
    for j = 1:length(algName)
        p{j} = plot(baseGraphList,metricData{i}(j,:),'LineWidth',2,'LineStyle',lineStyleList{j},'Marker',markerList{j},'Color',colorList{j});
    end
    legend([p{1},p{2},p{3},p{4},p{5},p{6}],algName{1},algName{2},algName{3},algName{4},algName{5},algName{6});
    ylabel(['\fontname{times new roman}',metricName{i}],'FontSize',15);
    xlabel('\fontname{times new roman}Arriving graph','FontSize',15);
    set(gca,'xtick',baseGraphList);
    a = 1;
    hold off;
end