
function [] = adi_nested_cv(X_train, trialinfo)

trialinfo = rmfield(trialinfo, 'X_test');

for kk=1:length(trialinfo.X_train)
    trialinfo.subject{kk} = trialinfo.X_train{kk}{4,2};
    trialinfo.subject_num(kk) = str2double(trialinfo.subject{kk}(end-1:end));
    if 1 == strcmp(trialinfo.X_train{kk}{2,2}, 'like') || 1 == strcmp(trialinfo.X_train{kk}{2,2}, 'Neu_Like')
       trialinfo.condition(kk) = 1;
    else 
       trialinfo.condition(kk) = 2;
    end
end
num_holdout = 1; % es soll nur 1 Proband pro like/dislike im Testset sein
CV_nested = adi_crossval(trialinfo.subject_num(find(trialinfo.condition==1)), trialinfo.subject_num(find(trialinfo.condition==2)), 'holdout', num_holdout);

X_train_like = X_train(find(trialinfo.condition==1),:,:);
X_train_dislike = X_train(find(trialinfo.condition==2),:,:);
clear X_train
trialinfo.like.subject_num = trialinfo.subject_num(find(trialinfo.condition==1));
trialinfo.dislike.subject_num = trialinfo.subject_num(find(trialinfo.condition==2));
trialinfo.like.info = trialinfo.X_train(find(trialinfo.condition==1));
trialinfo.dislike.info = trialinfo.X_train(find(trialinfo.condition==2));



for cc_nested = 1:3%length(CV_nested.pairs)
    X_test_nested  = cat(1, X_train_like(find(ismember(trialinfo.like.subject_num, CV_nested.like(CV_nested.pairs(cc_nested,1)).test(1))),:,:), ...
              X_train_dislike(find(ismember(trialinfo.dislike.subject_num, CV_nested.dislike(CV_nested.pairs(cc_nested,2)).test(1))),:,:));

    trialinfo.nested(cc_nested).X_test  = cat(2, trialinfo.like.info(find(ismember(trialinfo.like.subject_num, CV_nested.like(CV_nested.pairs(cc_nested,1)).test(1)))), ...
              trialinfo.dislike.info(find(ismember(trialinfo.dislike.subject_num, CV_nested.dislike(CV_nested.pairs(cc_nested,2)).test(1)))));              

    clabel_test_nested = [ones(1, size(X_train_like(find(ismember(trialinfo.like.subject_num, CV_nested.like(CV_nested.pairs(cc_nested,1)).test(1)))),2)) ... ...
              2*ones(1,size(X_train_dislike(find(ismember(trialinfo.dislike.subject_num, CV_nested.dislike(CV_nested.pairs(cc_nested,2)).test(1)))),2))];

    X_train_like_nested = X_train_like;
    X_train_like_nested(find(ismember(trialinfo.like.subject_num, CV_nested.like(CV_nested.pairs(cc_nested,1)).test(1))),:,:) = []; 
    trialinfo.nested(cc_nested).X_train_like = trialinfo.like.info;
    trialinfo.nested(cc_nested).X_train_like(find(ismember(trialinfo.like.subject_num, CV_nested.like(CV_nested.pairs(cc_nested,1)).test(1)))) = []; 

    X_train_dislike_nested = X_train_dislike;
    X_train_dislike_nested(find(ismember(trialinfo.dislike.subject_num, CV_nested.dislike(CV_nested.pairs(cc_nested,2)).test(1))),:,:) = []; 
    trialinfo.nested(cc_nested).X_train_dislike = trialinfo.dislike.info;
    trialinfo.nested(cc_nested).X_train_dislike(find(ismember(trialinfo.dislike.subject_num, CV_nested.dislike(CV_nested.pairs(cc_nested,2)).test(1)))) = []; 

    X_train_nested = cat(1, X_train_like_nested, X_train_dislike_nested) ;
    clabel_train_nested = [ones(1, size(X_train_like_nested,1)) 2*ones(1, size(X_train_dislike_nested,1))];
    
  
    %% PCA % Frage an Stefan: PCA für jede Epoche separat?  oder feature selection by adding a sparsity promoting penalty term like l1-cost?
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4040248/ => kann man
% searchlight für feature reduction verwenden anstelle von PCA =>
% accuracies in p-WErte umwandeln, und nur signifikante p-werte verwenden
% => ist aber nicht so sensitiv wie andere tests...
% => aber dann nur für relevante Zeitbereiche?
% feature weight analysis for svm => aber nur für interessierende
% Zeitpunkte?d
    
    
for pp = 1:size(X_train_nested,1)
    train_nested.trial{pp} = squeeze(X_train_nested(pp,:,:));
    train_nested.time{pp} = trialinfo.time;
end
    train_nested.label = trialinfo.label;
    train_nested.fsample = trialinfo.fsample;
    train_nested.grad = trialinfo.grad;
    train_nested.grad.chantype(249:end)=[];
    train_nested.grad.chantype = train_nested.grad.chantype';
    train_nested.grad.chanunit(249:end)=[];
    
       
for pp = 1:size(X_test_nested,1)
    test_nested.trial{pp} = squeeze(X_test_nested(pp,:,:));
    test_nested.time{pp} = trialinfo.time;
end
    test_nested.label = trialinfo.label;
    test_nested.fsample = trialinfo.fsample;
    test_nested.grad = trialinfo.grad;
    test_nested.grad.chantype(249:end)=[];
    test_nested.grad.chantype = train_nested.grad.chantype;
    test_nested.grad.chanunit(249:end)=[];
    
    
    %% PCA
    
    cfg = [];
    cfg.method = 'pca';
    cfg.demean = 'yes';
    [comp_train, E_train, d_train] = ft_componentanalysis_kh(cfg, train_nested);
    
    
    %%  trials auch demeanen:
    
    
    for kk = 1:length(train_nested.trial)
        for pp = 1:size(train_nested.trial{kk},1)
            train_nested.trial{kk}(pp,:) = train_nested.trial{kk}(pp,:) - mean(train_nested.trial{kk}(pp,:));
        end
    end
    
    for kk = 1:length(train_nested.trial)
        trial{kk} = comp_train.unmixing*train_nested.trial{kk};
    end
    
    train_nested.trial = trial;
    clearvars trial X_train_nested X_train_like_nested X_train_dislike_nested comp_train
    
    %% test set: falsch!
     for kk = 1:length(test_nested.trial)
        for pp = 1:size(test_nested.trial{kk},1)
            test_nested.trial{kk}(pp,:) = test_nested.trial{kk}(pp,:) - mean(test_nested.trial{kk}(pp,:));
        end
    end
    
    for kk = 1:length(test_nested.trial)
        trial{kk} = comp_test.unmixing*test_nested.trial{kk};
    end
    
    test_nested.trial = trial;

    
    
    
    
    %% %% Get default hyperparameters for the logreg and lda classifier
   
    param_logreg = mv_get_classifier_param('logreg');
    
    param_lda = mv_get_classifier_param('lda');
    
    param_svm = mv_get_classifier_param('svm');
    
    data_train_nested = zeros(length(train_nested.trial), size(train_nested.trial{1,1},1), size(train_nested.trial{1,1},2));
    for kk = 1:length(train_nested.trial)
        data_train_nested(kk,:,:) = train_nested.trial{kk};
    end
    
    data_test_nested = zeros(length(test_nested.trial), size(test_nested.trial{1,1},1), size(test_nested.trial{1,1},2));
    for kk = 1:length(test_nested.trial)
        data_test_nested(kk,:,:) = test_nested.trial{kk};
    end
    
    %% 
    
    
    for tt = 1:length(trialinfo.time)
        cf_logreg{cc_nested, tt} = train_logreg(param_logreg, data_train_nested(:,:,tt), clabel_train_nested);
        [predlabel, dval, prob] = test_logreg(cf_logreg{cc_nested, tt}, data_test_nested(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean

        % Calculate AUC and Accuracy
        perf.logreg.auc(cc_nested, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test_nested);
        perf.logreg.accuracy(cc_nested, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test_nested);

    end
    
%%  lda for single time points  
    for tt = 1:length(trialinfo.time)
        cf_lda{cc_nested, tt} = train_lda(param_lda, data_train_nested(:,:,tt), clabel_train_nested);
        [predlabel, dval] = test_lda(cf_lda{cc_nested, tt}, data_test_nested(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.lda.auc(cc_nested, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test_nested);
        perf.lda.accuracy(cc_nested, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test_nested);
        
    end    
    
%%  svm for single time points  
    for tt = 1:length(trialinfo.time)
        cf_svm{cc_nested, tt} = train_svm(param_svm, data_train_nested(:,:,tt), clabel_train_nested);
        [predlabel, dval] = test_svm(cf_svm{cc_nested,tt}, data_test_nested(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.svm.auc(cc_nested, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test_nested);
        perf.svm.accuracy(cc_nested, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test_nested);
        
    end        
    
end


    %% 









end