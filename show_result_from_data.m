load('exp_online_willow_outlier_batch_1.mat');
graphRange = target.config.graphRange;
baseGraphCnt = target.config.baseGraphCnt;
ave.accuracy = accAveFull;
ave.score = scrAveFull;
ave.consistency = conPairAveFull;
ave.time = timAveFull;
ave.matchingNumber = countPairAveFull;
legendOff = 0;
fields = fieldnames(ave);
figure(1);
for ifield = 1:length(fields)
    xtag='Arriving graph';ytag=[fields{ifield}];
    plotResult_new(legendOff,graphRange-baseGraphCnt+1, getfield(ave,fields{ifield}), algSet, xtag, ytag);
end