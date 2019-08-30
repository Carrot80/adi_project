load 'E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data\gbf_ball\zscore\MVPA\holdout_subj\feature_reduction\tests\eigenvectors_E_train.mat';
load 'E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data\gbf_ball\zscore\MVPA\holdout_subj\feature_reduction\tests\sorted_eigenvalues_d_train.mat';

param_logreg = mv_get_classifier_param('logreg');  
param_lda = mv_get_classifier_param('lda'); 
param_svm = mv_get_classifier_param('svm');

explained_variance = [0.95 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1]; 
sum_d_train = sum(d_train(:,2));
ratio_d_train = d_train(:,2)/sum_d_train;
ratio_d_train_sum = sum(ratio_d_train);

for kk=1:length(explained_variance)
    
    sum_explained = 0;
    idx = 0;
    while sum_explained < explained_variance(kk)
        idx = idx + 1;
        sum_explained = sum_explained + ratio_d_train(idx);
    end
    kept_comp_train.(['exp_var_' num2str(explained_variance(kk)*100)]) = d_train(1:idx,1);

end


%% test

sum_d_test = sum(d_test(:,2));
ratio_d_test = d_test(:,2)/sum_d_test;
ratio_d_test_sum = sum(ratio_d_test);

for kk_explained_var=1:length(explained_variance)
    
    sum_explained = 0;
    idx = 0;
    while sum_explained < explained_variance(kk_explained_var)
        idx = idx + 1;
        sum_explained = sum_explained + ratio_d_test(idx);
    end
    kept_comp_test.(['exp_var_' num2str(explained_variance(kk_explained_var)*100)]) = d_test(1:idx,1);

end
%% demeaning
    % train  
    for ff = 1:length(train_nested.trial)
        for pp = 1:size(train_nested.trial{ff},1)
            train_nested_demean.trial{ff}(pp,:) = train_nested.trial{ff}(pp,:) - mean(train_nested.trial{ff}(pp,:));
        end
    end

    for ff = 1:length(test_nested.trial)
        for pp = 1:size(test_nested.trial{ff},1)
            test_nested_demean.trial{ff}(pp,:) = test_nested.trial{ff}(pp,:) - mean(test_nested.trial{ff}(pp,:));
        end
    end
    
%%
% return the desired number of principal components

for kk_explained_var = 1:length(explained_variance)
    
    unmixing_train = E_train(:,d_train(1:length(kept_comp_train.(['exp_var_' num2str(explained_variance(kk_explained_var)*100)]))))';
    unmixing_test = E_test(:,d_test(1:length(kept_comp_train.(['exp_var_' num2str(explained_variance(kk_explained_var)*100)]))))';

    for ff = 1:length(train_nested_demean.trial)
        trials_train{ff} = unmixing_train*train_nested_demean.trial{ff};
    end
    
%     train_nested_demean.trial = trial;
%     clear trial
    
    
        %% test set
    
    for ff = 1:length(test_nested_demean.trial)
        trials_test{ff} = unmixing_test*test_nested_demean.trial{ff};
    end
    
%     test_nested_demean.trial = trial;
%     clear trial
%%
  
    data_train_nested = zeros(length(trials_train), size(trials_train{1,1},1), size(trials_train{1,1},2));
    for ff = 1:length(trials_train)
        data_train_nested(ff,:,:) = trials_train{ff};
    end
    
    data_test_nested = zeros(length(trials_test), size(trials_test{1,1},1), size(trials_test{1,1},2));
    for ff = 1:length(trials_test)
        data_test_nested(ff,:,:) = trials_test{ff};
    end

    clear trials_test trials_train
    
     for tt = 1:length(trialinfo.time)
        cf_lda{kk_explained_var, tt} = train_lda(param_lda, data_train_nested(:,:,tt), clabel_train_nested);
        [predlabel, dval] = test_lda(cf_lda{kk_explained_var, tt}, data_test_nested(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.lda.auc(kk_explained_var, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test_nested);
        perf.lda.accuracy(kk_explained_var, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test_nested);
        
    end    
    
    
    
end


figure; hold on;
axis tight
plot(trialinfo.time, 0.5*ones(1, length(trialinfo.time)), '--')
plot(trialinfo.time, perf.lda.auc(1,:), 'r')
plot(trialinfo.time, perf.lda.auc(2,:), 'k')
plot(trialinfo.time, perf.lda.auc(3,:), 'g')
plot(trialinfo.time, perf.lda.auc(4,:), 'y')
plot(trialinfo.time, perf.lda.auc(5,:), 'b')

plot(trialinfo.time, perf.lda.auc)
plot(trialinfo.time, mean(perf_acc_svm))
axis tight
ylim([0.3 0.8])
legend({' ', 'lda', 'logreg', 'svm'})
title('mean accuracy')
savefig([path2save '\mean_accuracy.fig'])