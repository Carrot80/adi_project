function [] = adi_mvpa_sensorspace_holdout(balldesign, fn_input, input_path, fn2save, path2save, config)


%% hold out: leave two subjects out

if ~exist(path2save, 'dir')
    mkdir(path2save)
end

if ~exist([path2save fn2save], 'file')
    load ([input_path fn_input])
    if 0 == isfield(sensordata_all_subj_zscore, 'trial_info')
        sensordata_all_subj_zscore = setfield(sensordata_all_subj_zscore,{1}, 'trial_no', []);
        sensordata_all_subj_zscore = setfield(sensordata_all_subj_zscore,{1}, 'trial_info', []);
        for k = 1:length(sensordata_all_subj_zscore)
            for p = 1:length(sensordata_all_subj_zscore(k).trial)
                sensordata_all_subj_zscore(k).trial_no{1, p} = num2str(p);
                sensordata_all_subj_zscore(k).trial_info{1, p} = {'trial_no' sensordata_all_subj_zscore(k).trial_no{1, p}; ...
                'response_label' sensordata_all_subj_zscore(k).response_label{1, p}; 'balldesign_short' sensordata_all_subj_zscore(k).balldesign_short{1, p}; ...
                'subject' sensordata_all_subj_zscore(k).subject{1, p}; 'run' sensordata_all_subj_zscore(k).run{1, p} };
            end
        end
        save([input_path fn_input], 'sensordata_all_subj_zscore')
    end
    
   session = sensordata_all_subj_zscore(1);
   for k = 2:length(sensordata_all_subj_zscore)
       session.trial = cat(2,session.trial, sensordata_all_subj_zscore(k).trial);
       session.time = cat(2,session.time, sensordata_all_subj_zscore(k).time);
       session.response_label = cat(2,session.response_label, sensordata_all_subj_zscore(k).response_label);
       session.balldesign_short = cat(2,session.balldesign_short, sensordata_all_subj_zscore(k).balldesign_short);
       session.subject = cat(2,session.subject, sensordata_all_subj_zscore(k).subject);
       session.trial_info = cat(2,session.trial_info, sensordata_all_subj_zscore(k).trial_info);
   end


    for k=1:length(session.response_label)
        switch session.response_label{k}
            case 'like' % 'Volley'
                session.labels(k) = 1;
            case 'Neu_Like' % 'Volley'
                session.labels(k) = 1;
            case 'dislike' % 'Space'
                session.labels(k) = 2;
            case 'Neu_Dislike' % 'Space'
                session.labels(k) = 2;
        end
    end
    clear sensordata_all_subj_zscore
    
    mvpa_cv_feature_red(session, path2save, fn2save, config) 
else
    mk_figures(fn2save, path2save) 
end

end


function mvpa_cv_feature_red(session, path2save, fn2save, config)

if isempty(session)
    mk_figures(path2save)
    
else

[like, dislike] = split_data(session, 'like_dislike');
clear session

data_like = zeros(length(like.trial), size(like.trial{1,1},1), size(like.trial{1,1},2));
for k=1:length(like.trial)
    data_like(k,:,:) = like.trial{k};
end
% like = rmfield(like, 'trial'); 

data_dislike = zeros(length(dislike.trial), size(dislike.trial{1,1},1), size(dislike.trial{1,1},2));
for k=1:length(dislike.trial)
    data_dislike(k,:,:) = dislike.trial{k};
end
% dislike = rmfield(dislike, 'trial'); 
rng('default')

%% Crossvalidation folds
num_holdout = 2;
CV = adi_crossval(like.subject_num, dislike.subject_num, 'holdout', num_holdout);

%% Get default hyperparameters for the logreg and lda classifier
   
    param_logreg = mv_get_classifier_param('logreg');
    
    param_lda = mv_get_classifier_param('lda');
    
    param_svm = mv_get_classifier_param('svm');

    
%% built train and test data based on CV folds
cf_logreg = cell(length(CV.pairs),length(like.time{1}));
cf_lda = cell(length(CV.pairs),length(like.time{1}));
cf_svm = cell(length(CV.pairs),length(like.time{1}));

for cc = 1:length(CV.pairs)
    
    if 2 == length(CV.like(CV.pairs(cc,1)).test)  && 1 == length(CV.dislike(CV.pairs(cc,2)).test)
        X_test  = cat(1, data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:,:), ...
                  data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2))),:,:), ...
                  data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:,:));   
              
        X_test_trialinfo  = cat(1, like.trial_info(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))), ...
                  like.trial_info(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2))),:,:), ...
                  dislike.trial_info(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:,:));             
              
        
        clabel_test = [ones(1, size(data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))),2)) ...
                  ones(1, size(data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))),2)) ...
                  2*ones(1,size(data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))),2))];
        X_train_like = data_like;
        X_train_like([find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))) find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))],:,:) = []; 
        X_train_like_trialinfo = like.trialinfo;
        X_train_like_trialinfo([find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))) find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))]) = []; 
   
        X_train_dislike = data_dislike;
        X_train_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:,:) = []; 
        X_train_dislike_trialinfo = dislike.trialinfo;
        X_train_dislike_trialinfo(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))) = []; 
        
    elseif 2 == length(CV.dislike(CV.pairs(cc,2)).test) && 1 == length(CV.like(CV.pairs(cc,1)).test)
        X_test  = cat(1, data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:,:), ...
                  data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:,:), ...
                  data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2))),:,:));
              
        X_test_trialinfo  = cat(1, like.trial_info(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))), ...
                  dislike.trial_info(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))), ...
                  dislike.trial_info(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))));                 
              
        
        clabel_test = [ones(1, size(data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))),2)) ...
                      2*ones(1,size(data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))),2))...
                      2*ones(1,size(data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))),2))];
                  
        X_train_like = data_like;
        X_train_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:,:) = []; 
        X_train_like_trialinfo = like.trial_info;
        X_train_like_trialinfo(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))) = []; 

        X_train_dislike = data_dislike;
        X_train_dislike([find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))) find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))],:,:) = []; 
        X_train_dislike_trialinfo = dislike.trial_info;
        X_train_dislike_trialinfo([find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))) find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))]) = []; 

        
    elseif 1 == length(CV.dislike(CV.pairs(cc,2)).test) && 1 == length(CV.like(CV.pairs(cc,1)).test)
        
        X_test  = cat(1, data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:,:), ...
                  data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:,:));
        
        X_test_trialinfo  = cat(1, like.trial_info(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))), ...
                  dislike.trial_info(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))));              
              
        clabel_test = [ones(1, size(data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))),2)) ... ...
                  2*ones(1,size(data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))),2))];
 
        X_train_like = data_like;
        X_train_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:,:) = []; 
        X_train_like_trialinfo = like.trial_info;
        X_train_like_trialinfo(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))) = []; 
        
        X_train_dislike = data_dislike;
        X_train_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:,:) = []; 
        X_train_dislike_trialinfo = dislike.trial_info;
        X_train_dislike_trialinfo(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))) = []; 
        
    else
        X_test  = cat(1, data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:,:), ...
                  data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2))),:,:), ...
                  data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:,:), ...
                  data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2))),:,:));
        
        X_test_trialinfo = cat(2, like.trial_info(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))), ...
                  like.trial_info(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))), ...
                  dislike.trial_info(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))), ...
                  dislike.trial_info(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))));
          
              
        clabel_test = [ones(1, size(data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))),2)) ...
                  ones(1, size(data_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))),2)) ...
                  2*ones(1,size(data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))),2))...
                  2*ones(1,size(data_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))),2))];
        
        X_train_like = data_like;
        X_train_like([find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))) find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))],:,:) = []; 
        X_train_like_trialinfo = like.trial_info;
        X_train_like_trialinfo([find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))) find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))]) = []; 

        X_train_dislike = data_dislike;
        X_train_dislike([find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))) find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))],:,:) = []; 
        X_train_dislike_trialinfo = dislike.trial_info;
        X_train_dislike_trialinfo([find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))) find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))]) = []; 
      
    end
    
    trialinfo(cc).X_train = cat(2, X_train_like_trialinfo, X_train_dislike_trialinfo) ;
    trialinfo(cc).X_test = X_test_trialinfo;
    X_train = cat(1, X_train_like, X_train_dislike) ;
    clabel_train = [ones(1, size(X_train_like,1)) 2*ones(1, size(X_train_dislike,1))];
%     clear X_train_like X_train_dislike
%     clear X_train_like_trialinfo X_train_dislike_trialinfo X_test_trialinfo
 
%%  balance data  (evtl. nach PCA durchführen)
if strcmp(config.balance, 'undersample')
    [X_train, clabel_train, trialinfo(cc).X_train, labelidx_train] = kh_mv_balance_classes(X_train,clabel_train,'undersample',[], trialinfo(cc).X_train);
    [X_test, clabel_test, trialinfo(cc).X_test, labelidx_test] = kh_mv_balance_classes(X_test,clabel_test,'undersample',[], trialinfo(cc).X_test);
end

%% PCA in nested crossvalidation:

trialinfo.time = like.time{1,1};
trialinfo.fsample = like.fsample;
trialinfo.grad = like.grad;
trialinfo.label = like.label;
adi_nested_cv(X_train, trialinfo);

    
    %% pca with matlab

%     Ntrials  = length(data_train.trial);
%     Nsamples = zeros(1,Ntrials);
%     Nchans = size(data_train.trial{1,1},1);
%     % determine the size of each trial, they can be variable length
% 
%     for trial=1:Ntrials
%       Nsamples(trial) = size(data_train.trial{trial},2);
%     end
%     
%         % concatenate data:
%       dat = zeros(Nchans, sum(Nsamples));
%       for trial=1:Ntrials
%         ft_info('.');
%         begsample = sum(Nsamples(1:(trial-1))) + 1;
%         endsample = sum(Nsamples(1:trial));
%         dat(:,begsample:endsample) = data_train.trial{trial};
%       end
%       
%       [coeff,score,latent,tsquared,explained] = pca(dat);
%       
%     % ind the number of components required to explain at least 95% variability by using a while loop  
%     
%     sum_explained = 0;
%     idx = 0;
%     while sum_explained < 95
%         idx = idx + 1;
%         sum_explained = sum_explained + explained(idx);
%     end
%     idx
%     
%     data_train_95
    
     
        
%%  logistic regression for single time points  
    for tt = 1:length(like.time{1})
        cf_logreg{cc, tt} = train_logreg(param_logreg, X_train(:,:,tt), clabel_train);
        [predlabel, dval, prob] = test_logreg(cf_logreg{cc, tt}, X_test(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.logreg.auc(cc, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test);
        perf.logreg.accuracy(cc, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test);
        
    end
    
%%  lda for single time points  
    for tt = 1:length(like.time{1})
        cf_lda{cc, tt} = train_lda(param_lda, X_train(:,:,tt), clabel_train);
        [predlabel, dval] = test_lda(cf_lda{cc, tt}, X_test(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.lda.auc(cc, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test);
        perf.lda.accuracy(cc, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test);
        
    end    
    
%%  svm for single time points  
    for tt = 1:length(like.time{1})
        cf_svm{cc, tt} = train_svm(param_svm, X_train(:,:,tt), clabel_train);
        [predlabel, dval] = test_svm(cf_svm{cc,tt}, X_test(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.svm.auc(cc, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test);
        perf.svm.accuracy(cc, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test);
        
    end        
    
    clear X_train clabel_train Xtest clabel_test
    
end
end

result_feature_red = [];
result_feature_red.logreg.cf = cf_logreg;
result_feature_red.logreg.perf = perf.logreg;
result_feature_red.lda.cf = cf_lda;
result_feature_red.lda.perf = perf.lda;
result_feature_red.svm.cf = cf_svm;
result_feature_red.svm.perf = perf.svm;
result_feature_red.CV = CV;
result_feature_red.balance = 'undersample';
result_feature_red.time = like.time{1};
result_feature_red.trialinfo = trialinfo;

save([path2save 'result.mat'], 'result')

%% figures 

if 1 == iscell(result_feature_red.lda.perf.accuracy)
    for k=1:length(result_feature_red.CV.pairs)
        perf_acc_lda(k,:) = [result_feature_red.lda.perf.accuracy{k,:}];
        perf_auc_lda(k,:) = [result_feature_red.lda.perf.auc{k,:}];
        perf_acc_logreg(k,:) = [result_feature_red.logreg.perf.accuracy{k,:}];
        perf_auc_logreg(k,:) = [result_feature_red.logreg.perf.auc{k,:}];
        perf_acc_svm(k,:) = [result_feature_red.svm.perf.accuracy{k,:}];
        perf_auc_svm(k,:) = [result_feature_red.svm.perf.auc{k,:}];
    end
else
    for k=1:length(result_feature_red.CV.pairs)
        perf_acc_lda(k,:) = result_feature_red.lda.perf.accuracy(k,:);
        perf_auc_lda(k,:) = result_feature_red.lda.perf.auc(k,:);
        perf_acc_logreg(k,:) = result_feature_red.logreg.perf.accuracy(k,:);
        perf_auc_logreg(k,:) = result_feature_red.logreg.perf.auc(k,:);
        perf_acc_svm(k,:) = result_feature_red.svm.perf.accuracy(k,:);
        perf_auc_svm(k,:) = result_feature_red.svm.perf.auc(k,:);
    end
end

std_perf_acc_lda = std(perf_acc_lda);
std_perf_auc_lda = std(perf_auc_lda);
std_perf_acc_logreg = std(perf_acc_logreg);
std_perf_auc_logreg = std(perf_auc_logreg);
std_perf_acc_svm = std(perf_acc_svm);
std_perf_auc_svm = std(perf_auc_svm);

%%
figure; hold on;
axis tight
plot(result.time, 0.5*ones(1, length(result_feature_red.time)), '--')
plot(result.time, mean(perf_acc_lda))
plot(result.time, mean(perf_acc_logreg))
plot(result.time, mean(perf_acc_svm))
axis tight
ylim([0.3 0.8])
legend({' ', 'lda', 'logreg', 'svm'})
title('mean accuracy')
savefig([path2save '\mean_accuracy.fig'])
hold off
close 

figure; hold on;
axis tight
plot(result.time, 0.5*ones(1, length(result.time)), '--')
plot(result.time, mean(perf_auc_lda))
hold on
plot(result.time, mean(perf_auc_logreg))
hold on
plot(result.time, mean(perf_auc_svm))
legend({'', 'lda', 'logreg', 'svm'})
title('mean AUC')
ylim([0.3 0.8])
savefig([path2save '\mean_auc.fig'])
close 
% 
% figure;
% plot(result.time, mean(result.lda.perf.auc))
% hold on
% plot(result.time, mean(result.logreg.perf.auc))
% hold on
% plot(result.time, mean(result.svm.perf.auc))
% legend({'lda', 'logreg', 'svm'})
% title('mean area under curve')
% savefig([path2save 'mean_auc.fig'])
% close

% fprintf('Classification accuracy: %0.2f\n', stat.accuracy)

end



function [] = mk_figures(fn2save, path2save)

load ([path2save 'result.mat'])

if 1 == iscell(result.lda.perf.accuracy)
    for k=1:length(result.CV.pairs)
        perf_acc_lda(k,:) = [result.lda.perf.accuracy{k,:}];
        perf_auc_lda(k,:) = [result.lda.perf.auc{k,:}];
        perf_acc_logreg(k,:) = [result.logreg.perf.accuracy{k,:}];
        perf_auc_logreg(k,:) = [result.logreg.perf.auc{k,:}];
        perf_acc_svm(k,:) = [result.svm.perf.accuracy{k,:}];
        perf_auc_svm(k,:) = [result.svm.perf.auc{k,:}];
    end
else
    for k=1:length(result.CV.pairs)
        perf_acc_lda(k,:) = result.lda.perf.accuracy(k,:);
        perf_auc_lda(k,:) = result.lda.perf.auc(k,:);
        perf_acc_logreg(k,:) = result.logreg.perf.accuracy(k,:);
        perf_auc_logreg(k,:) = result.logreg.perf.auc(k,:);
        perf_acc_svm(k,:) = result.svm.perf.accuracy(k,:);
        perf_auc_svm(k,:) = result.svm.perf.auc(k,:);
    end
end

std_perf_acc_lda = std(perf_acc_lda);
std_perf_auc_lda = std(perf_auc_lda);
std_perf_acc_logreg = std(perf_acc_logreg);
std_perf_auc_logreg = std(perf_auc_logreg);
std_perf_acc_svm = std(perf_acc_svm);
std_perf_auc_svm = std(perf_auc_svm);

% standard error:
std_error_perf_acc_lda = std(perf_acc_lda)./sqrt(size(perf_acc_lda,1));
std_error_perf_auc_lda = std(perf_auc_lda)./sqrt(size(perf_auc_lda,1));
std_error_perf_acc_logreg = std(perf_acc_logreg)./sqrt(size(perf_acc_logreg,1));
std_error_perf_auc_logreg = std(perf_auc_logreg)./sqrt(size(perf_auc_logreg,1));
std_error_perf_acc_svm = std(perf_acc_svm)./sqrt(size(perf_acc_svm,1));
std_error_perf_auc_svm = std(perf_auc_svm)./sqrt(size(perf_auc_svm,1));


%% plot AUC of LDA => min, max, std
transparency = 0.5;
close all
figure; hold on;
filled_min_max=[max(perf_auc_lda),fliplr(min(perf_auc_lda))];
xpoints=[result.time,fliplr(result.time)];
fillhandle_min_max=fill(xpoints,filled_min_max,'green');%plot the data
set(fillhandle_min_max,'EdgeColor',[0.1 0.1 0.1],'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
filled_std=[(mean(perf_auc_lda)+std(perf_auc_lda)),fliplr(mean(perf_auc_lda)-std(perf_auc_lda))];
fillhandle_std=fill(xpoints,filled_std,'blue');%plot the data
set(fillhandle_std,'EdgeColor',[0.1 0.1 0.1],'FaceAlpha',0.2,'EdgeAlpha',0.2);%set edge color
plot(result.time, mean(perf_auc_lda), 'Color', 'black', 'LineWidth',1.5);
plot(result.time, 0.5*ones(1, length(result.time)),  'Color','blue', 'LineStyle', '--');
ylim([0 1])
line([0 0],ylim, 'Color','blue', 'LineStyle','--')
axis tight
% plot(result.time, max(perf_auc_lda),'Color', 'green');
% plot(result.time, min(perf_auc_lda),'Color', 'green');
set(gcf,'Position',[722 376 1060 420]);
% legend_fig = legend({'lda', 'std'}, 'Location','northeast', 'box', 'off');
title('LDA - mean AUC, std, min, max')
savefig([path2save '\LDA_mean_AUC.fig'])
hold off
close 


%% plot AUC of LDA => Convidence Intervals

for kk = 1:length(perf_auc_lda)
    ci_lda_auc(:,kk) = ci(perf_auc_lda(:,kk));
end

transparency = 0.5;
close all
figure; hold on;
filled_ci=[ci_lda_auc(2,:),fliplr(ci_lda_auc(1,:))];
xpoints=[result.time,fliplr(result.time)];
fillhandle_ci=fill(xpoints,filled_ci,'blue');%plot the data
set(fillhandle_ci,'EdgeColor',[0.1 0.1 0.1],'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
plot(result.time, mean(perf_auc_lda), 'Color', 'black', 'LineWidth',1.5);
plot(result.time, 0.5*ones(1, length(result.time)),  'Color','blue', 'LineStyle', '--');
ylim([0 1])
line([0 0],ylim, 'Color','blue', 'LineStyle','--')
axis tight
% plot(result.time, max(perf_auc_lda),'Color', 'green');
% plot(result.time, min(perf_auc_lda),'Color', 'green');
set(gcf,'Position',[722 376 1060 420]);
% legend_fig = legend({'lda', 'std'}, 'Location','northeast', 'box', 'off');
title('LDA - mean AUC + CI')
savefig([path2save '\LDA_mean_AUC_plusCI.fig'])
hold off
close 

figure;
plot(result.time, perf_auc_lda)
title('LDA - AUC per CV fold')
savefig([path2save '\LDA_AUC_per_fold.fig'])


end



