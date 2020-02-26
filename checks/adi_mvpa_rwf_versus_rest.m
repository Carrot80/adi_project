function [] = adi_leave_out_exemplar_singlsubj(subjectdir, balldesign, timerange)

%%
% created: 06.02.2020
path2save = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\group_rwf_vs_non_rwf\';
if 2 == exist([path2save 'session_all_subj_1.mat'], 'file')
    mvpa_leave_out_subject(path2save, balldesign, timerange)
    
else
    
    for ii = 1:length(subjectdir)

        if 2 == exist([subjectdir(ii).folder filesep subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\session.mat'], 'file')

            load ([subjectdir(ii).folder filesep subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\session.mat'])

        else

            dir_data = dir([subjectdir(ii).folder filesep subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\02_interpolated\' '*.mat']);
            ind_ = [];
            for kk = 1:length(dir_data)
                 ind_(kk) = contains(dir_data(kk).name, 'Dontcare');
            end

            indx_dontcare = find(ind_);
            if ~isempty(indx_dontcare)
                dir_data(indx_dontcare) = [];
            end
            for kk = 1:length(dir_data)

                load ([dir_data(kk).folder filesep dir_data(kk).name], 'cleanMEG_interp')
                for pp = 1:length(cleanMEG_interp.trial)
                    cleanMEG_interp.trialinfo.run{pp} = dir_data(kk).name(end-4);
                end
                if 1 == isfield(cleanMEG_interp, 'additional_cleaning')
                    cleanMEG_interp = rmfield(cleanMEG_interp,  'additional_cleaning');
                end
                [data_bpfreq_res_sel] = adi_bpfilter(cleanMEG_interp, 'bp1_45Hz');

                data(kk) = data_bpfreq_res_sel;
                clear cleanMEG_interp data_bpfreq_res_sel
            end

            session = data(1);
            session.balldesign = data(1).trialinfo.balldesign_short;
            session = rmfield(session, 'trialinfo');
            session = rmfield(session, 'ChannelFlag_Bst');
            session = rmfield(session, 'cfg');
            session = rmfield(session, 'grad');

            for kk = 2:length(data)
                session.trial = cat(2, session.trial, data(kk).trial);
                session.time = cat(2, session.time, data(kk).time);
                session.balldesign = cat(2, session.balldesign, data(kk).trialinfo.balldesign_short);
            end
            clear data



          %% z-transform trials  

            [session] = adi_ztrans_sensorspace(session);

        end

        for kk = 1:length(session.balldesign) 
            balldesign_array(kk)= session.balldesign{kk};
        end

        for kk = 1:numel(session.balldesign)
            session.subject{1, kk} = subjectdir(ii).name;
        end

        for kk = 1:numel(session.balldesign)
            switch balldesign_array{kk}
                case 'rwf'
                    session.condition_rwf{1, kk} = 'rwf';
                    session.label_rwf{1, kk} = 1;
                otherwise
                    session.condition_rwf{1, kk} = 'non_rwf';
                    session.label_rwf{1, kk} = 0;
            end
        end

        switch ii
            case 1
                session_all_subj = session;
            otherwise
                session_all_subj.trial = cat(2, session_all_subj.trial, session.trial);
                session_all_subj.time = cat(2, session_all_subj.time, session.time);
                session_all_subj.condition_rwf = cat(2, session_all_subj.condition_rwf, session.condition_rwf);
                session_all_subj.label_rwf = cat(2, session_all_subj.label_rwf, session.label_rwf);
                session_all_subj.subject = cat(2, session_all_subj.subject, session.subject);
        end

        clear session
    end

    session_all_subj_1 = session_all_subj;
    session_all_subj_1.trial(2001:end) = [];
    session_all_subj_1.time(2001:end) = [];
    session_all_subj_1.subject(2001:end) = [];
    session_all_subj_1.condition_rwf(2001:end) = [];
    session_all_subj_1.label_rwf(2001:end) = [];
    save('session_all_subj_1.mat', 'session_all_subj_1')
    session_all_subj_2 = session_all_subj;
    session_all_subj_2.trial(1:2000) = [];
    session_all_subj_2.time(1:2000) = [];
    session_all_subj_2.subject(1:2000) = [];
    session_all_subj_2.condition_rwf(1:2000) = [];
    session_all_subj_2.label_rwf(1:2000) = [];
    save('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\group_rwf_vs_non_rwf\session_all_subj_2.mat', 'session_all_subj_2')
   
end

tic
load([path2save 'session_all_subj_1.mat'], 'session_all_subj_1')
load([path2save 'session_all_subj_2.mat'], 'session_all_subj_2')
toc
session = session_all_subj_1; 
session.trial = cat(2, session.trial, session_all_subj_2.trial);
session.condition_rwf = cat(2, session.condition_rwf, session_all_subj_2.condition_rwf);
session.subject = cat(2, session.subject, session_all_subj_2.subject);
session.label_rwf = cell2mat(cat(2, session.label_rwf, session_all_subj_2.label_rwf));
if any(session.label_rwf < 1) 
   session.label_rwf(find(session.label_rwf<1)) = 2;
end
clear session_all_subj_1 session_all_subj_2

data_trials = kh_trial2dat(session.trial);
save('data_trials.mat', 'data_trials', '-v7.3')
MatObj = matfile('data_trials.mat', 'Writable',true);



%% MVPA
%    [perf, perf_perm] = mvpa_leave_out_balldesign(session, balldesign, timerange);
%     perf.config.methods = 'perf_without_pca_bootstrap_4pseudotrials_per_condition_in_trainfold';
%     perf.config.comment = 'da testset trials mit sowohl like als auch dislike-Antworten enthalten kann, wurden nur trials aus dem trainfold zur Bildung von Pseudotrials gemittelt. Um unbalanced data zu balancieren, wurde oversampling der kleineren Bedingung (Like/dislike) durchgef�hrt';
%     perf.config.feature_selection = 'no feature selection';
%     perf.ratings = like_dislike_ratings;
%     
%     perm_gbf = perf_perm.gbf;
%     perm_gbs = perf_perm.gbs;
%     perm_gbv = perf_perm.gbv;
%     perm_ggf = perf_perm.ggf;
%     perm_ggs = perf_perm.ggs;
%     perm_ggv = perf_perm.ggv;
%     perm_rwf = perf_perm.rwf;
%     perm_rws = perf_perm.rws;
%     perm_rwv = perf_perm.rwv;
%     
%     if ~exist(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing'], 'dir')
%         mkdir(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing'])
%     end
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_real_data.mat'], 'perf')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_gbf.mat'], 'perm_gbf')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_gbs.mat'], 'perm_gbs')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_gbv.mat'], 'perm_gbv')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_ggf.mat'], 'perm_ggf')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_ggs.mat'], 'perm_ggs')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_ggv.mat'], 'perm_ggv')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_rwf.mat'], 'perm_rwf')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_rws.mat'], 'perm_rws')
%     save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_rwv.mat'], 'perm_rwv')
%         
%     clear perm_gbf perm_gbs perm_gbv perm_ggf perm_ggs perm_ggv perm_rwf perm_rws perm_rwv 
%  
%     
%     %%
%     clear session perf perf_perm
      
   
   end





function [] = mvpa_leave_out_subject(path2save, balldesign, timerange)


%% Get default hyperparameters for the logreg and lda classifier
load([path2save 'session.mat'], 'session')

%% built train and test data based on CV folds
%% Crossvalidation folds
CV = adi_crossval_leaveExemplarOut(session.subject, session.condition_rwf);

%%
MatObj = matfile([path2save 'data_trials.mat'], 'Writable',true);
time = session.time{1};

%% data_bootstrap_trainfold verwenden f�r die Klassifikation, so dass nur labels getauscht werden und die Daten nicht ver�ndert werden
perf = struct('mean_f1score', [], 'mean_acc', [], 'CI_acc', [], 'SEM', []);
perf_perm = struct('mean_f1score', [], 'mean_acc', [], 'CI_acc', [], 'SEM', []);
parfor ii = 1:numel(time)   
   
    [perf(ii), perf_perm(ii)] = compute_folds_and_mvpa(CV, MatObj, timerange, path2save, ii)
    
end
perf.CV = CV;
 


end
function  [perf_ii, perf_perm_ii] = compute_folds_and_mvpa(CV, MatObj, timerange, path2save, ii)
load([path2save 'session.mat'], 'session')
param_lda = mv_get_hyperparameter('lda');
time = session.time{1};
indx_time_1 = nearest(time, timerange(1)); 
indx_time_2 = nearest(time, timerange(2)); 

subject = unique(session.subject);
data_trials = MatObj.data_trials(:,:,ii); 

for kk = 1:length(CV) 

    % labels des testset werden beibehalten
    test_fold = data_trials(CV(kk).testset,:);
    clabel_test_fold = session.label_rwf(CV(kk).testset);

    train_fold = data_trials;
    train_fold(CV(kk).testset,:) = [];  
    clabel_train_fold = session.label_rwf;
    clabel_train_fold(CV(kk).testset) = [];

    trials_trainfold_balldesign1 = train_fold(clabel_train_fold == 1,:);
    trials_trainfold_balldesign2_9 = train_fold(clabel_train_fold == 2,:);
    max_numel_trials = max([size(trials_trainfold_balldesign1,1), size(trials_trainfold_balldesign2_9,1)]);
    train_fold = [];
    
    cfg = [];
    cfg.mode = 'bootstrap';
    cfg.averages = 4;
    cfg.repetitions = max_numel_trials;

    bootstrap_trainfold_balldesign1 = fte_subaverage(cfg, trials_trainfold_balldesign1);
    bootstrap_trainfold_balldesign2_9 = fte_subaverage(cfg, trials_trainfold_balldesign2_9);
    trials_trainfold_balldesign1 = [];
    trials_trainfold_balldesign2_9 = [];

    data_bootstrap_trainfold = cat(1, bootstrap_trainfold_balldesign1, bootstrap_trainfold_balldesign2_9);
    clabel_train_fold = cat(2, ones(1, size(bootstrap_trainfold_balldesign1,1)), 2*ones(1, size(bootstrap_trainfold_balldesign2_9,1)));
    bootstrap_trainfold_balldesign1 = [];
    bootstrap_trainfold_balldesign2_9 = [];

%%  lda for single time points  
    disp(['starting lda on real data for testdata ' subject{kk}]) 

    cf_lda = train_lda(param_lda, data_bootstrap_trainfold(:,:), clabel_train_fold);

    [predlabel, dval] = test_lda(cf_lda, test_fold(:,:)); 
    % To calculate classification accuracy, compare the predicted labels to
    % the true labels and take the mean

    % Calculate AUC and Accuracy
    accuracy = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    confusion_percentage = mv_calculate_performance('confusion', 'clabel', predlabel', clabel_test_fold);
    confusion = confusionmat(clabel_test_fold, predlabel');
    f1score = mv_calculate_performance('f1', 'clabel', predlabel', clabel_test_fold);

    lda{kk}.confusion = confusion;
    lda{kk}.accuracy = accuracy;
    lda{kk}.confusion_percentage = confusion_percentage;
    lda{kk}.f1score = f1score;
%     perf.(subject{kk}).number_of_trials.testset = numel(clabel_test_fold);  
%     perf.(subject{kk}).number_of_trials.trainingsset = numel(clabel_train_fold);

%%  
    size_random_dataset = 1000;
    clabel_train_fold_rand = zeros(size_random_dataset, length(clabel_train_fold));

    for pp = 1:size_random_dataset
        clabel_train_fold_rand(pp,:) = clabel_train_fold(randperm(length(clabel_train_fold)));
    end
    if ii >= indx_time_1 && ii <= indx_time_2
        disp(['creating null distribution for subject ' subject{kk}])

        accuracy_perm_ii{kk} = nan(size_random_dataset,1);
        f1score_perm_ii{kk} = nan(size_random_dataset,1);

        for pp = 1:size_random_dataset
            cf_perm_lda = train_lda(param_lda, data_bootstrap_trainfold, clabel_train_fold_rand(pp,:));
            [predlabel_perm, dval_perm] = test_lda(cf_perm_lda, test_fold); 

            % Calculate AUC and Accuracy and other performance measures:
            accuracy_perm_ii{kk}(pp) = mv_calculate_performance('acc', 'dval', dval_perm, clabel_test_fold);
%             perf_perm_ii.(subject{kk}).confusion_percentage{pp} = mv_calculate_performance('confusion', 'clabel', predlabel_perm', clabel_test_fold);
%             confusion_perm_ii.(subject{kk}){pp} = confusionmat(clabel_test_fold, predlabel_perm');
            f1score_perm_ii{kk}(pp) = mv_calculate_performance('f1', 'clabel', predlabel_perm', clabel_test_fold);  
%             perf_perm_ii.(subject{kk}).number_of_trials.testset = numel(clabel_test_fold);  
%             perf_perm_ii.(subject{kk}).number_of_trials.trainingsset = numel(clabel_train_fold);
        end
        mean_accuracy_perm{kk} = nanmean(accuracy_perm_ii{kk});
        mean_f1_perm{kk} = nanmean(f1score_perm_ii{kk});
        
    else
        perf_perm_ii = struct('mean_f1score', [], 'mean_acc', [], 'CI_acc', [], 'SEM', []);
    end
    test_fold = [];
    clabel_train_fold = [];
    clabel_test_fold = [];
end

for kk = 1:length(CV)
    accuracy(kk) = lda{kk}.accuracy;
    f1score(kk) = lda{kk}.f1score;
end

perf_ii.mean_f1score = nanmean(f1score);
perf_ii.SEM = std(accuracy)/sqrt(length(accuracy));  
perf_ii.mean_acc = nanmean(accuracy);
perf_ii.CI_acc = ci(accuracy);

if exist('mean_accuracy_perm', 'var')
    for kk = 1:length(CV)
        perm_accuracy(kk) = mean_accuracy_perm{kk};
        perm_f1score(kk) = mean_f1_perm{kk};
    end
    perf_perm_ii.mean_acc = nanmean(perm_accuracy);
    perf_perm_ii.mean_f1score = nanmean(perm_f1score);
    perf_perm_ii.CI_acc = ci(perm_accuracy);
    perf_perm_ii.SEM = std(perm_accuracy)/sqrt(length(perm_accuracy));  
end
end


function [CV] = adi_crossval_leaveExemplarOut(subject_array, condition)

[no_trls_subject, subject] = histcounts(categorical(subject_array), categorical(unique(subject_array)));

for kk = 1:length(subject)
    index_design.(subject{kk}) = find(strcmp(subject_array, subject(kk)));
end

for kk = 1:length(subject)
    CV(kk).testset = index_design.(subject{kk});
    CV(kk).trainingsset = 1:length(subject_array);
    CV(kk).trainingsset(CV(kk).testset) = [];
    CV(kk).design = subject{kk};
end

for kk = 1:length(subject)
    CV(kk).labels_trainingsset = condition;
    CV(kk).labels_trainingsset(CV(kk).testset) = [];
    CV(kk).labels_testset = condition(index_design.(subject{kk}));
    CV(kk).condition_trainingsset = subject_array;
    CV(kk).condition_trainingsset(CV(kk).testset) = [];
    CV(kk).trialnumer_rwfs_trainingsset = numel(find(strcmp(CV(kk).labels_trainingsset, 'rwf')));
    CV(kk).trialnumer_non_rwfs_trainingsset = numel(find(strcmp(CV(kk).labels_trainingsset, 'non_rwf')));
    CV(kk).ratio_trialnumer_likes_dislikes_trainingsset = CV(kk).trialnumer_rwfs_trainingsset./(CV(kk).trialnumer_rwfs_trainingsset+CV(kk).trialnumer_non_rwfs_trainingsset);
end

end


function [data_bpfreq_res_sel] = adi_bpfilter(filename, bpname)


switch bpname
    case 'bp1_95Hz'
       data_bpfreq_res_sel =  filename;
        return
    case 'bp1_45Hz'
        bpfreq = [1 45];
    case'1_5_45Hz'
        bpfreq = [1.5 45];
    case 'bp2_45Hz'
        bpfreq = [2 45];
    case 'bp3_45Hz'
        bpfreq = [3 45];
    case 'delta'
        bpfreq = 4;
    case 'theta'
        bpfreq = [4 8];
    case 'alpha'
        bpfreq = [8 13];
    case 'beta'
        bpfreq = [13 25];
    case 'low_gamma'
        bpfreq = [25 45];
    case 'high_gamma'
        bpfreq = [55 90];
    case 'bp10-45Hz'
        bpfreq = [10 45];
end

cfg = [];
cfg.keeptrials = 'yes';
cfg.vartrllength = 2;

cfg.trials  = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') %|| 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

try
    [data_bpfreq] = ft_preprocessing(cfg, filename); 
    [warnMsg, warnID] = lastwarn;
    if ~isempty(warnMsg)
       warnMsg
    end
catch
    
    warnMsg
end

cfg =[];
cfg.resamplefs = 256;
% cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);

for k = 1:length(data_bpfreq_res_sel.trial)
    data_bpfreq_res_sel.grad.label(249:end) = [];
    data_bpfreq_res_sel.grad.chanori(249:end, :) = [];
    data_bpfreq_res_sel.grad.chanpos(249:end, :) = [];
    data_bpfreq_res_sel.grad.tra(249:end, :) = [];
    data_bpfreq_res_sel.label(249:end) = [];
end

fn_data_bpfreq_sel_res = fieldnames(data_bpfreq_res_sel);
fn_filename = fieldnames(filename);

diff_fieldnames = setdiff(fn_filename, fn_data_bpfreq_sel_res);

for k=1:length(diff_fieldnames)
    data_bpfreq_res_sel.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
end
% fn_filename{end+1}='cfg';
data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, fn_filename);
clearvars filename data_bpfreq data_bpfreq_res 

end



