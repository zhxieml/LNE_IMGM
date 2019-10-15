clear all;
init_path;
set_conf_img;

if 1
[ classes dataset ] = init_dataset_img(conf);
    save(getDatasetPath(conf), 'classes', 'dataset');
else
    load(getDatasetPath(conf), 'classes', 'dataset');
end

nameExp = 'demo';
target_class = conf.class;
nExp = conf.nExp;

conf.bOW_dataset = false;
conf.bOW_model = true;
conf.bOW_result = true;
conf.verbose = false;

conf.referenceType = 'sample';  
conf.learningType = 'non';        
confs{1} = conf;
conf.referenceType = 'null';    
conf.learningType = 'HARG';
confs{2} = conf;
nConf = 2;

table_exp_TRAIN = zeros(nConf,2,nExp);
table_exp_TEST = zeros(nConf,2,nExp);

for iConf = 1:nConf
    if ~isempty(confs{iConf})
        [ accuracy error_EP score_GM ] = run_HARG_img('conf',confs{iConf},'overwrite',false);
        table_exp_TRAIN(iConf,1,:) = accuracy(1,:);
        table_exp_TRAIN(iConf,2,:) = error_EP(1,:);
        table_exp_TEST(iConf,1,:) = accuracy(2,:);
        table_exp_TEST(iConf,2,:) = error_EP(2,:);
    end
end

fprintf('\n\n *** %10s SUMMARY: average performance of %d experiments\n', target_class, nExp);
for iConf = 1:nConf
    if ~isempty(confs{iConf})
        fprintf('%10s +%8s on TRAIN - accuracy: %.4f    error: %.4f\n', ...
            confs{iConf}.learningType,confs{iConf}.referenceType, mean(table_exp_TRAIN(iConf,1,:)), mean(table_exp_TRAIN(iConf,2,:)) );
        fprintf('%10s +%8s on TEST  - accuracy: %.4f    error: %.4f\n', ...
            confs{iConf}.learningType, confs{iConf}.referenceType, mean(table_exp_TEST(iConf,1,:)), mean(table_exp_TEST(iConf,2,:)) );
    end
end
fprintf('\n\n');

save(getExpPath(conf, nameExp), 'classes', 'dataset', 'confs', 'table_exp_TRAIN', 'table_exp_TEST' );

