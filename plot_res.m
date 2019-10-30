%% Initialization
clear *;clear -global *;close all;clc; clear all;

resDir = './res';
class = 'Duck';

inCnt = 10;
outCnt = 0;
test_i = 0;

loadPath = fullfile(resDir, class, sprintf('in%02d_out%02d_%02d.mat', inCnt, outCnt, test_i));
if isfile(loadPath)
    load(loadPath, 'accResult', 'timeResult', 'scrResult', 'conResult', 'matchResult');
end

incrementCnt = size(accResult, 2);

setPlotColor;

algNameSepSpace = '                    ';
algSet.algNameSet = {'cao','cao_raw','IMGM','IMGM_new'};
algSet.algEnable = [1,1,1,1];
algSet.algColor = {caoClr,cao_rawClr,IMGMClr,IMGM_newClr};
algSet.algLineStyle = {'--',':','-','-.'};
algSet.algMarker = {'.','.','.','.'};%
[~,caoIdx] = ismember('cao',algSet.algNameSet);
[~,rawCaoIdx] = ismember('cao_raw',algSet.algNameSet);
[~,IMGMIdx] = ismember('IMGM',algSet.algNameSet);
[~,IMGMNewIdx] = ismember('IMGM_new',algSet.algNameSet);
algSet.algNameSetDisplay{caoIdx} = 'cao';
algSet.algNameSetDisplay{rawCaoIdx} = 'raw_cao';
algSet.algNameSetDisplay{IMGMIdx} = 'IMGM';
algSet.algNameSetDisplay{IMGMNewIdx} = 'IMGM_new';

%% Plot
plotResult_new(0, 1:incrementCnt, timeResult', algSet, 'increment #', 'time');
plotResult_new(0, 1:incrementCnt, accResult', algSet, 'increment #', 'acc');
plotResult_new(0, 1:incrementCnt, scrResult', algSet, 'increment #', 'src');
plotResult_new(0, 1:incrementCnt, conResult', algSet, 'increment #', 'con');

plot(1:incrementCnt, matchResult(:, 1), '-o');


pause;