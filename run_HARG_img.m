function [ accuracy error_EP score_GM ] = run_HARG_img( varargin )
% RUN_HARG_IMG
% by Minsu CHO

%randn('state',0) ;
%rand('state',0) ;
init_path; 
set_conf_img;
bOverwrite = false;  % rerun and overwrite data

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'conf'
      conf = arg;    
    case 'overwrite'
      bOverwrite = arg; 
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

% --------------------------------------------------------------------
%                                                          Setup data
% --------------------------------------------------------------------

%initialize dataset and save it for reuse
if ~exist(getDatasetPath(conf)) || bOverwrite || conf.bOW_dataset
    [ classes dataset ] = init_dataset_img(conf);
    save(getDatasetPath(conf), 'classes', 'dataset');
else
    load(getDatasetPath(conf));
end

% ------------------------------------------------------------------
%                                                Learning SVM struct
% ------------------------------------------------------------------
if ~exist(getModelPath(conf)) || bOverwrite || conf.bOW_model
    load(getDatasetPath(conf), 'classes', 'dataset');
    for iExp=1:min(conf.nExp,numel(dataset))
        model = init_model_img( dataset{iExp}, conf );
        model = init_svm_img( model ); 
        save(getModelPath(conf, iExp),'model');
        model = learn_demo_img(dataset{iExp}, getModelPath(conf, iExp), conf);
        save(getModelPath(conf, iExp),'model');
        if conf.verbose
            show_modelA_img( iExp, 'conf', conf);
            drawnow;
        end
    end
end

% ------------------------------------------------------------------
%                                                              Plots
% ------------------------------------------------------------------

if ~exist(getResultPath(conf)) || bOverwrite || conf.bOW_result
      for iExp=1:min(conf.nExp,numel(dataset))
        fprintf('\n\n--- learn ''%s'' with %s on the TRAIN set\n', conf.learningType, conf.referenceType);
        
        if conf.verbose % for comparison
            tconf = conf;
            tconf.learningType = 'non';
            tconf.referenceType = 'sample';
            fprintf('sample w/o learning\n');
            test_demo_img(iExp,'datatype',0,'conf',tconf);
        end
        
         % actual performance
        fprintf('after %s learning\n', conf.learningType);
        [ accuracy(1,iExp) error_EP(1,iExp) score_GM(1,iExp) ] = ...
            test_demo_img(iExp,'datatype',0,'conf',conf);
        
        fprintf('\n--- learn ''%s'' with %s on the TEST set\n', conf.learningType, conf.referenceType);
        
        if conf.verbose % for comparison
            tconf = conf;
            tconf.learningType = 'non';
            tconf.referenceType = 'sample';
            fprintf('sample w/o learning\n');
            test_demo_img(iExp,'datatype',1,'conf',tconf);
        end
        
        % actual performance
        fprintf('after %s learning\n',conf.learningType);
        [ accuracy(2,iExp) error_EP(2,iExp) score_GM(2,iExp) ] = ...
            test_demo_img(iExp,'datatype',1,'conf',conf);
        
        %save(getResultPath(conf,classes{ci}),'res');
      end
      
      fprintf('\n\n === average performance of %d experiments\n', size(accuracy,2));
      fprintf('TRAIN - accuracy: %.4f    error: %.4f\n', mean(accuracy(1,:)), mean(error_EP(1,:)) );
      fprintf('TEST  - accuracy: %.4f    error: %.4f\n', mean(accuracy(2,:)), mean(error_EP(2,:)) );
      
end

end

