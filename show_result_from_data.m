graphMinCnt = 20;
graphMaxCnt = 30;
graphStep = 1;
baseGraphCnt = graphMinCnt;
graphRange = baseGraphCnt:graphStep:graphMaxCnt-graphStep;
load('willow_outlier.mat');
ave.accuracy = accAveFull;
ave.score = scrAveFull;
ave.consistency = conPairAveFull;
ave.time = timAveFull;
ave.matchingNumber = countPairAveFull;
set_
legendOff = 0;
fields = fieldnames(ave);
figure(1);
for ifield = 1:length(fields)
    xtag='Arriving graph';ytag=[fields{ifield}, '(',target.config.category,')'];
    plotResult_new(legendOff,graphRange-baseGraphCnt+1, getfield(ave,fields{ifield}), algSet, xtag, ytag);
end