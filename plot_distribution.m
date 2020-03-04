

distriPath = "distribution_Duck_from20to50_30.mat";
distriPicFile = distriPath(1:end-4);
load(distriPath, 'distribution');
% distribution = distribution(1:20, :);

% for i = 1:50
%     disp(mean(distribution(i, 30, 10), 2));
% end
% for i = 1:30
%     disp(distribution(49, i, 10));
% end
distriAve = zeros(30, 10);
distriAve(:, :) = mean(distribution(1:5, :, :), 1);

algSet.algNameSet = {'[0.0, 0.1]', '(0.1, 0.2]', '(0.2, 0.3]', '(0.3, 0.4]', '(0.4, 0.5]', '(0.5, 0.6]', '(0.6, 0.7]', '(0.7, 0.8]', '(0.8, 0.9]', '(0.9, 1.0]'};
algSet.algEnable =  [ 1,            1,             1,            1,             1,              1,              1,             1,              1,                1];
algSet.algColor = { 'black', [0.04,0.04,0.72], [0.00,0.27,1.00], [0.00,0.55,1.00], [0.00,0.84,1.00], [1.00,0.84,0.00], [0.85,0.55,0.00], [0.92,0.42,0.23], [0.83,0.27,0.27], [1, 0, 0]};
algSet.algLineStyle = {'-',       '-',           '-',            '-',           '-',         '-.',             '-.',         '-.',           '-.',            '--'};
algSet.algMarker =    {'.',       '.',            '.',             '.',            '.',          '.',              '.',          '.',            '.',           '.'};

% rename the algorithm names to be more friendly
for i=1:length(algSet.algNameSet)
    if strcmp(target.config.testType,'massOutlier')
        algSet.algNameSetDisplay{cao_Idx} = 'cao^{cst}';
        algSet.algNameSetDisplay{i} = strrep(algSet.algNameSet{i},'_s','^{sim}');
        algSet.algNameSetDisplay{i} = strrep(algSet.algNameSetDisplay{i},'o_','o-');
        algSet.algNameSetDisplay{i} = strrep(algSet.algNameSetDisplay{i},'c_','c^{cst}');
    else
        algSet.algNameSetDisplay{i} = strrep(algSet.algNameSet{i},'_','-');
        if algSet.algNameSetDisplay{i}(end)=='-'
            algSet.algNameSetDisplay{i}(end) = '*';
        end
    end
end

legendOff = 0;
ave.Frequency = distriAve;
fields = fieldnames(ave);
for ifield = 1:length(fields)
    xtag='Arriving graph';ytag=fields{ifield};
    plotResult_new(legendOff,(20:39)-baseGraphCnt+1, getfield(ave,fields{ifield}), algSet, xtag, ytag);
end