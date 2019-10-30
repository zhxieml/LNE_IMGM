%% dataset papameters 
conf.dataset = 'WILLOWOBJ';
conf.dataDir = 'C:\Users\ChenZixuan\Documents\data\WILLOW-ObjectClass-dataset';
conf.imgDir=[conf.dataDir '/WILLOW-ObjectClass'];
conf.annoDir=[conf.dataDir '/WILLOW-ObjectClass'];
conf.featDir = [conf.dataDir '/feature'];
conf.affinityDir=[conf.dataDir '/affinity'];
conf.gtDir = [conf.dataDir '/ground_truth'];
conf.resDir = './res';
conf.tmpDir = './tmp';
conf.numInlier = 10;
conf.numOutlier = 0;
conf.distRatioTrue = 0.15; % threshold ratio to decide the true match
                           % ratio of the object size 
%conf.overwrite = false ;
conf.prefix = 'DEMO' ;
conf.learningType = 'HARG';           % HARG, DW, SW, non, SW_marius
conf.referenceType = 'null';      % sample, source, null

%classes = dir(conf.imgDir);
%classes = classes([classes.isdir]) ;
%classes = {classes(3:end).name} ;

%conf.class = 'Face';
conf.class = 'Duck';
%conf.class = 'Car';
%conf.class = 'Motorbike';
%conf.class = 'Winebottle';

conf.nExp = 1; % # of random trials for the class

%randn('state',conf.randSeed) ;
%rand('state',conf.randSeed) ;
%vl_twister('state',conf.randSeed) ;

conf.bOW_dataset = false;  % overwrite dataset initialization
conf.bOW_model = true;     % overwrite model learning
conf.bOW_result = true;    % overwrite test results

conf.verbose = false;
