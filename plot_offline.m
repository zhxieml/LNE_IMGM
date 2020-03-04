%% Initialization
clear *;clear -global *;close all;clc;
global target
init_path;
setPlotColor;
algpar = setPairwiseSolver();
setObsoleteVariables;
accAveFull = zeros(30, 8);

% target.config.graphMinCnt=19; 
% target.config.graphMaxCnt=50; 
target.config.testCnt = 20;% v
target.config.maxNumSearch = 20;
target.config.batchSize = 4;

target.config.graphMinCnt=16; 
target.config.graphMaxCnt=32; 
graphMaxCnt = target.config.graphMaxCnt;
graphMinCnt = target.config.graphMinCnt;
batchSize = target.config.batchSize;

graphCntRange = graphMinCnt:batchSize:graphMaxCnt;
% graphStep = 1;
target.config.database = "willow"; % "willow", "synthetic", "CMU-sequence"
target.config.category = 'outlier';

% set algorithms
algNameSepSpace = '                    ';
algSet.algNameSet = {'cao_c_raw','imgm_d','imgm_r','anc_imgm','anc_imgm_a16','anc_imgm_d16', 'anc_imgm_a32','anc_imgm_d32'};
algSet.algEnable =  [     1,        1,       1,         1,          0,              0,             1,              1];
algSet.algColor = { cao_c_rawClr,imgm_dClr,imgm_rClr,anc_imgmClr,anc_imgmaClr16,anc_imgmdClr16, anc_imgmaClr16,anc_imgmdClr16};
algSet.algLineStyle = {'-',        '--',      '--',      '-.',        '-.',            '-.',          '-.'            '-.'};
algSet.algMarker = {   '.',        '.',       '.',       '.',        '.',            '.',              '.',            '.'};

% rename the algorithm names to be more friendly
for i=1:length(algSet.algNameSet)
    algSet.algNameSetDisplay{i} = strrep(algSet.algNameSet{i},'_','-');
    if algSet.algNameSetDisplay{i}(end)=='-'
        algSet.algNameSetDisplay{i}(end) = '*';
    end
end

legendOff = 0;
ave.accuracy = accAveFull;
fields = fieldnames(ave);
for ifield = 1:length(fields)
    xtag='number of graph';ytag=[fields{ifield}, '(',target.config.category,')'];
    plotResult_new(legendOff,graphCntRange, getfield(ave,fields{ifield}), algSet, xtag, ytag);
end