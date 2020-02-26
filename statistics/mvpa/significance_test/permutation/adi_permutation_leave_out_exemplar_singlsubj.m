function [] = adi_leave_out_exemplar_singlsubj(subjectdir, balldesign, timerange)

%%
% created: 21.01.2020


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
        session.response_label = data(1).trialinfo.response_label;
        session.run = data(1).trialinfo.run;
        session.balldesign = data(1).trialinfo.balldesign_short;
        session = rmfield(session, 'trialinfo');

        for kk = 2:length(data)
            session.trial = cat(2, session.trial, data(kk).trial);
            session.time = cat(2, session.time, data(kk).time);
            session.cfg = cat(2, session.cfg, data(kk).cfg);
            session.response_label = cat(2, session.response_label, data(kk).trialinfo.response_label);
            session.balldesign = cat(2, session.balldesign, data(kk).trialinfo.balldesign_short);
            session.run = cat(2, session.run, data(kk).trialinfo.run);
        end
        clear data
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

      %% z-transform trials  

        [session] = adi_ztrans_sensorspace(session);
    
    end
    like_dislike_ratings = [];
    for kk = 1:length(session.balldesign) 
        balldesign_array(kk)= session.balldesign{kk};
    end
    
    [count_balldesign, balls] = histcounts(categorical(balldesign_array), categorical(unique(balldesign_array)));
    
    for kk = 1:length(balls)
        like_dislike_ratings.(balls{kk}) = session.labels(find(strcmp(balldesign_array, balls{kk})));
    end
    clear balldesign_array
    
    %% MVPA  
   [perf, perf_perm] = mvpa_leave_out_balldesign(session, balldesign, timerange);
    perf.config.methods = 'perf_without_pca_bootstrap_4pseudotrials_per_condition_in_trainfold';
    perf.config.comment = 'da testset trials mit sowohl like als auch dislike-Antworten enthalten kann, wurden nur trials aus dem trainfold zur Bildung von Pseudotrials gemittelt. Um unbalanced data zu balancieren, wurde oversampling der kleineren Bedingung (Like/dislike) durchgeführt';
    perf.config.feature_selection = 'no feature selection';
    perf.ratings = like_dislike_ratings;
    
    perm_gbf = perf_perm.gbf;
    perm_gbs = perf_perm.gbs;
    perm_gbv = perf_perm.gbv;
    perm_ggf = perf_perm.ggf;
    perm_ggs = perf_perm.ggs;
    perm_ggv = perf_perm.ggv;
    perm_rwf = perf_perm.rwf;
    perm_rws = perf_perm.rws;
    perm_rwv = perf_perm.rwv;
    
    if ~exist(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing'], 'dir')
        mkdir(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing'])
    end
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_real_data.mat'], 'perf')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_gbf.mat'], 'perm_gbf')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_gbs.mat'], 'perm_gbs')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_gbv.mat'], 'perm_gbv')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_ggf.mat'], 'perm_ggf')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_ggs.mat'], 'perm_ggs')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_ggv.mat'], 'perm_ggv')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_rwf.mat'], 'perm_rwf')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_rws.mat'], 'perm_rws')
    save(['E:\adidas\fieldtrip_Auswertung\single_subjects\' subjectdir(ii).name filesep 'MEG\mvpa\significance_testing\performance_permutation_null_distribution_rwv.mat'], 'perm_rwv')
        
    clear perm_gbf perm_gbs perm_gbv perm_ggf perm_ggs perm_ggv perm_rwf perm_rws perm_rwv 
 
    
    %%
    clear session perf perf_perm
      
   
   end




end

function [perf, perf_perm] = mvpa_leave_out_balldesign(session, balldesign, peaktime)

%% Get default hyperparameters for the logreg and lda classifier
param_lda = mv_get_hyperparameter('lda');

    %% built train and test data based on CV folds
cf_lda = cell(length(balldesign),length(session.time{1}));


%% Crossvalidation folds

CV = adi_crossval_leaveExemplarOut(session.balldesign, session.labels);
data_trials = kh_trial2dat(session.trial);

%%

indx_time_1 = nearest(session.time{1,1}, peaktime(1)); 
indx_time_2 = nearest(session.time{1,1}, peaktime(2));  

n_CV = length(CV);

%% data_bootstrap_trainfold verwenden für die Klassifikation, so dass nur labels getauscht werden und die Daten nicht verändert werden
   
for kk = 1:n_CV 
      
    % labels des testset werden beibehalten
    test_fold = data_trials(CV(kk).testset,:,:);
    clabel_test_fold = session.labels(CV(kk).testset);
   
    train_fold = data_trials;
    train_fold(CV(kk).testset,:,:) = [];  
    clabel_train_fold = session.labels;
    clabel_train_fold(CV(kk).testset) = [];
    
    %% MVPA 
       
        %% PCA (z-transformation noch einbauen)
%         [pca_tt.coeff,pca_tt.score,pca_tt.latent,pca_tt.tsquared,pca_tt.explained] = pca(train_fold(:,:,tt), 'rows', 'all');
        
%         sum_explained = 0;
%         idx = 0;
%         while sum_explained < 95 
%             idx = idx + 1;
%             sum_explained = sum_explained + pca_tt.explained(idx);
%         end

%         train_pca = train_fold(:,:,tt)* pca_tt.coeff(:,1:idx);
%         test_pca = test_fold(:,:,tt)* pca_tt.coeff(:,1:idx);
        
%%  bootstrap trainfold    
    [trials_trainfold] = kh_trial2dat(train_fold);
 
    trials_trainfold_like = [];
    trials_trainfold_like.trial = trials_trainfold(clabel_train_fold == 1);
    trials_trainfold_like.time = session.time(1:length(trials_trainfold_like.trial));
    trials_trainfold_like.fsample = session.fsample;
    trials_trainfold_like.dimord = 'rpt_chan_time';
    trials_trainfold_like.istimelock = 1;
    trials_trainfold_like.label = session.label;
    
    trials_trainfold_dislike = [];
    trials_trainfold_dislike.trial = trials_trainfold(clabel_train_fold == 2);
    trials_trainfold_dislike.time = session.time(1:length(trials_trainfold_dislike.trial));
    trials_trainfold_dislike.fsample = session.fsample;
    trials_trainfold_dislike.dimord = 'rpt_chan_time';
    trials_trainfold_dislike.istimelock = 1;
    trials_trainfold_dislike.label = session.label;
    
    max_numel_trials = max([numel(trials_trainfold_like.trial), numel(trials_trainfold_dislike.trial)]);
    
    cfg = [];
    cfg.mode = 'bootstrap';
    cfg.averages = 4;
%     cfg.repetitions = max_numel_trials*2;
    cfg.repetitions = max_numel_trials;
    bootstrap_trainfold_like = fte_subaverage(cfg, trials_trainfold_like);
    bootstrap_trainfold_dislike = fte_subaverage(cfg, trials_trainfold_dislike);
    
    data_bootstrap_trainfold_like = kh_trial2dat(bootstrap_trainfold_like.trial);
    data_bootstrap_trainfold_dislike = kh_trial2dat(bootstrap_trainfold_dislike.trial);
    
    data_bootstrap_trainfold = cat(1, data_bootstrap_trainfold_like, data_bootstrap_trainfold_dislike);
    clabel_train_fold = cat(2, ones(1, size(data_bootstrap_trainfold_like,1)), 2*ones(1, size(data_bootstrap_trainfold_dislike,1)));
    
        %%  lda for single time points  
    disp(['starting lda on real data for balldesign ' balldesign{kk}]) 
    
    accuracy = nan(1, size(session.time{1},2));
%     auc = nan(1, size(session.time{1},2));
    confusion_percentage = cell(1, size(session.time{1},2));
    confusion = cell(1, size(session.time{1},2));
    f1score = cell(1, size(session.time{1},2));
%     sensitivity = cell(1, size(session.time{1},2));
%     precision_PPV = cell(1, size(session.time{1},2));
%     specificity_TNR = cell(1, size(session.time{1},2));
%     balanced_accuracy = cell(1, size(session.time{1},2));
%     dprime = cell(1, size(session.time{1},2));
%     tval = cell(1, size(session.time{1},2));
%     FDR = cell(1, size(session.time{1},2));
%     predicted_labels = cell(1, size(session.time{1},2));
%     testlabels = cell(1, size(session.time{1},2));
%     dvals = cell(1, size(session.time{1},2));
    
    tic
    parfor tt = 1:numel(session.time{1,1})
        cf_lda = train_lda(param_lda, data_bootstrap_trainfold(:,:,tt), clabel_train_fold);
        
        [predlabel, dval] = test_lda(cf_lda, test_fold(:,:,tt)); 
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
%         auc(1,tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
        accuracy(1, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
        confusion_percentage{1, tt} = mv_calculate_performance('confusion', 'clabel', predlabel', clabel_test_fold);
        confusion{1, tt} = confusionmat(clabel_test_fold, predlabel');
        f1score{1, tt} = mv_calculate_performance('f1', 'clabel', predlabel', clabel_test_fold);
%         sensitivity{1, tt} = mv_calculate_performance('recall', 'clabel', predlabel', clabel_test_fold);
%         precision_PPV{1, tt} = mv_calculate_performance('precision', 'clabel', predlabel', clabel_test_fold); 
%         try
%             specificity_TNR{1, tt} = confusion{1, tt}(2,2)./(confusion{1, tt}(1,2)+confusion{1, tt}(2,2)); 
%         catch
%             specificity_TNR{1, tt} = NaN;     
%         end
%         balanced_accuracy{1, tt} = (sensitivity{1, tt}+specificity_TNR{1, tt})./2;
%         tval{1, tt} = mv_calculate_performance('tval', 'dval', dval, clabel_test_fold);
%         dprime{1, tt} = norminv(sensitivity{1, tt})-norminv(1-specificity_TNR{1, tt});
%         FDR{1, tt} = 1-precision_PPV{1, tt};
%         predicted_labels{tt} = predlabel;
%         testlabels{tt} = clabel_test_fold;
%         dvals{1,tt} = dval; 
    end    
    toc
%     perf.(balldesign{kk}).lda.auc = auc;
    perf.(balldesign{kk}).lda.confusion = confusion;
    perf.(balldesign{kk}).lda.accuracy = accuracy;
    perf.(balldesign{kk}).lda.confusion_percentage = confusion_percentage;
    perf.(balldesign{kk}).lda.f1score = f1score;
%     perf.(balldesign{kk}).lda.sensitivity = sensitivity;
%     perf.(balldesign{kk}).lda.precision_PPV = precision_PPV;
%     perf.(balldesign{kk}).lda.specificity_TNR = specificity_TNR;
%     perf.(balldesign{kk}).lda.balanced_accuracy = balanced_accuracy;
%     perf.(balldesign{kk}).lda.tval = tval;
%     perf.(balldesign{kk}).lda.dprime = dprime;
%     perf.(balldesign{kk}).lda.FDR = FDR;
%     perf.(balldesign{kk}).lda.predicted_labels = predicted_labels;
%     perf.(balldesign{kk}).lda.testlabels = testlabels;
%     perf.(balldesign{kk}).lda.dvals = dvals;
    perf.(balldesign{kk}).number_of_trials.testset = numel(clabel_test_fold);  
    perf.(balldesign{kk}).number_of_trials.trainingsset = numel(clabel_train_fold);
    
    
   %%  
    size_random_dataset = 1000;
    clabel_train_fold_rand = zeros(size_random_dataset, length(clabel_train_fold));
    
    for pp = 1:size_random_dataset
        clabel_train_fold_rand(pp,:) = clabel_train_fold(randperm(length(clabel_train_fold)));
    end

    disp(['creating null distribution for balldesign ' balldesign{kk}])
    
    accuracy = nan(size_random_dataset, size(session.time{1},2));
%     auc = nan(size_random_dataset, size(session.time{1},2));
    confusion_percentage = cell(size_random_dataset, size(session.time{1},2));
    confusion = cell(size_random_dataset, size(session.time{1},2));
    f1score = cell(size_random_dataset, size(session.time{1},2));
%     sensitivity = cell(size_random_dataset, size(session.time{1},2));
%     precision_PPV = cell(size_random_dataset, size(session.time{1},2));
%     specificity_TNR = cell(size_random_dataset, size(session.time{1},2));
%     balanced_accuracy = cell(size_random_dataset, size(session.time{1},2));
%     dprime = cell(size_random_dataset, size(session.time{1},2));
%     tval = cell(size_random_dataset, size(session.time{1},2));
%     FDR = cell(size_random_dataset, size(session.time{1},2));

    tic
    parfor pp = 1:size(clabel_train_fold_rand,1)
        for tt = indx_time_1:indx_time_2
            cf_perm_lda = train_lda(param_lda, data_bootstrap_trainfold(:,:,tt), clabel_train_fold_rand(pp,:));
            [predlabel_perm, dval_perm] = test_lda(cf_perm_lda, test_fold(:,:,tt)); 

            % Calculate AUC and Accuracy and other performance measures:
%             auc(pp, tt) = mv_calculate_performance('auc', 'dval', dval_perm, clabel_test_fold);
            accuracy(pp, tt) = mv_calculate_performance('acc', 'dval', dval_perm, clabel_test_fold);
            confusion_percentage{pp, tt} = mv_calculate_performance('confusion', 'clabel', predlabel_perm', clabel_test_fold);
            confusion{pp, tt} = confusionmat(clabel_test_fold, predlabel_perm');
            f1score{pp, tt} = mv_calculate_performance('f1', 'clabel', predlabel_perm', clabel_test_fold);
%             sensitivity{pp, tt} = mv_calculate_performance('recall', 'clabel', predlabel_perm', clabel_test_fold);
%             precision_PPV{pp, tt} = mv_calculate_performance('precision', 'clabel', predlabel_perm', clabel_test_fold); 
%             try
%                 specificity_TNR{pp, tt} = confusion{pp, tt}(2,2)./(confusion{pp, tt}(1,2)+confusion{pp, tt}(2,2)); 
%             catch
%                 specificity_TNR{pp, tt} = NaN;
%             end
%             balanced_accuracy{pp, tt} = (sensitivity{pp, tt}+specificity_TNR{pp, tt})./2;
%             tval{pp, tt} = mv_calculate_performance('tval', 'dval', dval_perm, clabel_test_fold);
%             dprime{pp, tt} = norminv(sensitivity{pp, tt})-norminv(1-specificity_TNR{pp, tt});
%             FDR{pp, tt} = 1-precision_PPV{pp, tt};
%             predicted_label_perm{pp, tt} = predlabel_perm;
%             dvals_perm{pp, tt} = dval_perm;
                    
        end     
    end
    toc
%     perf_perm.(balldesign{kk}).auc = auc;
    perf_perm.(balldesign{kk}).accuracy = accuracy;
    perf_perm.(balldesign{kk}).confusion_percentage = confusion_percentage;
    perf_perm.(balldesign{kk}).confusion = confusion;
    perf_perm.(balldesign{kk}).f1score = f1score;
%     perf_perm.(balldesign{kk}).sensitivity = sensitivity;
%     perf_perm.(balldesign{kk}).precision_PPV = precision_PPV;
%     perf_perm.(balldesign{kk}).specificity_TNR = specificity_TNR;
%     perf_perm.(balldesign{kk}).balanced_accuracy = balanced_accuracy;
%     perf_perm.(balldesign{kk}).tval = tval;
%     perf_perm.(balldesign{kk}).dprime = dprime;
%     perf_perm.(balldesign{kk}).FDR = FDR;
%     perf_perm.(balldesign{kk}).timerange = indx_time_1:indx_time_2;
%     perf_perm.(balldesign{kk}).predicted_label_perm = predicted_label_perm;
%     perf_perm.(balldesign{kk}).dval_perm = dvals_perm;
%     
    perf_perm.(balldesign{kk}).number_of_trials.testset = numel(clabel_test_fold);  
    perf_perm.(balldesign{kk}).number_of_trials.trainingsset = numel(clabel_train_fold);
  
    clear test_fold train_fold clabel_train_fold clabel_test_fold
    
end

perf.CV = CV;
 
end

function [CV] = adi_crossval_leaveExemplarOut(balldesign_short, labels)

for kk=1:length(balldesign_short) 
    balldesign_array(kk) = balldesign_short{kk};
end

[count_balldesign, balls] = histcounts(categorical(balldesign_array), categorical(unique(balldesign_array)));

for kk = 1:length(balls)
    index_design.(balls{kk}) = find(strcmp(balldesign_array, balls(kk)));
end


for kk = 1:length(balls)
    CV(kk).testset = index_design.(balls{kk});
    CV(kk).trainingsset = 1:length(balldesign_short);
    CV(kk).trainingsset(CV(kk).testset) = [];
    CV(kk).design = balls{kk};
end

for kk = 1:length(balls)
    CV(kk).labels_trainingsset = labels;
    CV(kk).labels_trainingsset(CV(kk).testset) = [];
    CV(kk).labels_testset = labels(index_design.(balls{kk}));
    CV(kk).balldesign_trainingsset = balldesign_short;
    CV(kk).balldesign_trainingsset(CV(kk).testset) = [];
    CV(kk).trialnumer_likes_trainingsset = numel(find(CV(kk).labels_trainingsset==1));
    CV(kk).trialnumer_dislikes_trainingsset = numel(find(CV(kk).labels_trainingsset==2));
    CV(kk).ratio_trialnumer_likes_dislikes_trainingsset = CV(kk).trialnumer_likes_trainingsset/CV(kk).trialnumer_dislikes_trainingsset ;
    
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



