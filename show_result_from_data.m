
clear;
init_path;
savepath = 'C:\Users\zchen\Desktop\offline\duck\res_Duck.mat';
load(savepath);
setPlotColor;
algNameSet = { 'cao_pc_inc', 'cao_pc_raw', 'cao_c_inc','cao_c_raw','imgm_d','imgm_r','tbimgm_cao_c',...8
                'tbimgm_cao_pc','tbimgm_cao_cst', 'tbimgm_qm', 'tbimgm_matchALS', 'anc_imgm_a',...5
                'anc_imgm', 'anc_imgm_a4', 'anc_imgm_a8', 'anc_imgm_a16', 'anc_imgm_a32', 'anc_imgm_d4', 'anc_imgm_d8', 'anc_imgm_d16', 'anc_imgm_d32'...9
                'imgm_raw','imgm_cao','anc_imgm_raw','anc_imgm_cao'};%4

algReNameSet = {'cao_pc-inc', 'cao-pc-raw', 'cao-c-inc','cao-c-raw','imgm-d','imgm-r','lne-imgm',...8
                'lne-imgm-cao-pc','lne-imgm-cao-cst', 'lne-imgm-qm', 'lne-imgm-matchALS', 'lne-imgm-a',...5
                'lne-imgm', 'lne-imgm-a4', 'lne-imgm-a8', 'lne-imgm-a16', 'lne-imgm-a32', 'lne-imgm-d4', 'lne-imgm-d8', 'lne-imgm-d16', 'lne-imgm-d32'...9
                'imgm-raw','imgm-cao','lne-imgm-raw','lne-imgm-cao'};%4

algColor   = { cao_pcClr, cao_pc_rawClr, cao_cClr, cao_c_rawClr,imgm_dClr, imgm_rClr, lne_imgmClr, ...
               tbimgm_cao_pcClr, tbimgm_cao_cstClr, tbimgm_qmClr, tbimgm_matchALSClr, tbimgm_cao_c_adaClr...
               anc_imgmClr, anc_imgmaClr4, anc_imgmaClr8, anc_imgmaClr16, anc_imgmaClr32, anc_imgmdClr4, anc_imgmdClr8, anc_imgmdClr16, anc_imgmdClr32...
               imgm_dClr, imgm_caoClr, lne_imgmClr, anc_imgm_caoClr};

algLineStyle = {'-','-','-','-','-', '-','-','-','-','-', '-','-','-','-','-', '-','-','-','-','-', '-','-','-','-','-', '-'};
algMarker =    {'.','.','.','.','.', '.','.','.','.','.', '.','.','.','.','.', '.','.','.','.','.', '.','.','.','.','.', '.'};

for alg = find(algSet.algEnable)
    [Lia, idx] = ismember(algSet.algNameSet{alg}, algNameSet);
    if Lia
        algSet.algNameSetDisplay{alg} = algReNameSet{idx};
        algSet.algColor{alg} = algColor{idx};
        algSet.algLineStyle{alg} = algLineStyle{idx};
        algSet.algMarker{alg} = algMarker{idx};
    else
        error("undefined algorithms");
    end
end

%graphRange = target.config.graphRange;
%baseGraphCnt = target.config.baseGraphCnt;
graphRange = 16:4:32;
baseGraphCnt = 1;
ave.accuracy = accAveFull;
ave.score = scrAveFull;
ave.consistency = conPairAveFull;
ave.time = timAveFull;
ave.matchingNumber = countPairAveFull;
legendOff = 0;
fields = fieldnames(ave);

xTick = graphRange-baseGraphCnt+1;
% xNum = length(xTick);
% stepSize = xNum / 5;
% randomChoose = 1:5:xNum;
% xTick = xTick(randomChoose);

for ifield = 1:length(fields)
    xtag='Arriving graph';
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
    %data= getfield(ave,fields{ifield});
    %data = data(2:end, :);
    %xTick = xTick(1:end-1);
    plotResult(legendOff, xTick, getfield(ave,fields{ifield}), algSet, xtag, ytag);
end
