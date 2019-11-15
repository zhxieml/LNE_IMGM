
clear;
init_path;
setPlotColor;
algNameSet = { 'Synthetic', 'Car', 'Duck', 'Winebottle'};%4

algColor   = { 'r',         'y',    'b',    'g'};

algPath = [...
    "C:\Users\zchen\Desktop\exp1.4\exp_partition_synthetic_deform.mat", ...
    "C:\Users\zchen\Desktop\exp1.4\exp_partition_car_willow_outlier.mat",...
    "C:\Users\zchen\Desktop\exp1.4\exp_partition_Duck_willow_outlier.mat",...
    "C:\Users\zchen\Desktop\exp1.4\exp_partition_Winebottle_willow_outlier"...
];

algLineStyle = {'-','-','-','-','-', '-','-','-','-','-', '-','-','-','-','-', '-','-','-','-','-', '-','-','-','-','-', '-'};
algMarker =    {'.','.','.','.','.', '.','.','.','.','.', '.','.','.','.','.', '.','.','.','.','.', '.','.','.','.','.', '.'};

for alg = 1:4
    load(algPath(alg));
    [~, idx] = ismember(algSet.algNameSet, 'tbimgm_cao_c');
    idx = logical(idx);
    newAlgSet.algNameSetDisplay{alg} = algNameSet{alg};
    newAlgSet.algColor{alg} = algColor{alg};
    newAlgSet.algLineStyle{alg} = algLineStyle{alg};
    newAlgSet.algMarker{alg} = algMarker{alg};
    newAlgSet.algEnable(alg) = true;
    ave.accuracy(:, alg)       = accAveFull(:, idx);
    ave.score(:, alg)          = scrAveFull(:, idx);
    ave.consistency(:, alg)    = conPairAveFull(:, idx);
    ave.time(:, alg)           = timAveFull(:, idx);
    ave.matchingNumber(:, alg) = countPairAveFull(:, idx);
end

partitionRange = target.config.partitionRange;
baseGraphCnt = target.config.baseGraphCnt;
legendOff = 0;
fields = fieldnames(ave);
xTick = partitionRange;
% xNum = length(xTick);
% stepSize = xNum / 5;
% randomChoose = 1:5:xNum;
% xTick = xTick(randomChoose);

for ifield = 1:length(fields)
    xtag='M';
    name=[fields{ifield}];
    switch name
        case 'accuracy'
            ytag = 'Accuracy';
        case 'score'
            ytag = 'Affinity Score';
        case 'consistency'
            ytag = 'Consistency';
        case 'time'
            ytag = 'Time(s)';
        case 'matchingNumber'
            ytag = 'Pair';
        otherwise
            error('undefined field\n');
    end
    % data= getfield(ave,fields{ifield});
    % data = data(2:end, :);
    % xTick = xTick(1:end-1);
    plotResult(legendOff, xTick, getfield(ave,fields{ifield}), newAlgSet, xtag, ytag);
end
