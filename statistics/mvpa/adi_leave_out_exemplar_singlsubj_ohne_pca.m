function [] = adi_leave_out_exemplar_singlsubj(mainpath, subjectdir, balldesign)


for ii = 1:6% 10:length(subjectdir)
    
    if 2 == exist([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\session.mat'], 'file')
        
        load ([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\session.mat'])
        
    else
    
        dir_data = dir([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\02_interpolated\' '*.mat']);

        for kk = 1:length(dir_data)
            load ([dir_data(kk).folder filesep dir_data(kk).name], 'cleanMEG_interp')
            for pp = 1:length(cleanMEG_interp.trial)
                cleanMEG_interp.trialinfo.run{pp} = dir_data(kk).name(end-4);
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
  
    %% MVPA  
   [perf] = mvpa_leave_out_balldesign(session, balldesign);
%     perf.config.methods = 'bootstrap 2*trialnumber';
%     perf.config.methods = 'bootstrap_pseudotrials';
    perf.config.methods = 'perf_without_pca_bootstrap_4pseudotrials_trainfold_only';
    perf.config.comment = ' da alle trials mit den uneindeutigen Antworten zur Klassifikation genutzt worden, wurden nur trials aus dem trainfold zur Bildung von Pseudotrials gemittelt ';
    figure; hold on;
    plot(session.time{1,1}, mean(perf.lda.accuracy), 'k')
    plot(session.time{1,1}, mean(perf.svm.accuracy), 'b')
    plot(session.time{1,1}, mean(perf.logreg.accuracy), 'r')
    plot(session.time{1,1}, 0.5*ones(1, length(session.time{1,1})), '--')
    title('Accuracy without PCA')
    legend('lda', 'svm', 'logreg', ' ')
   
    
   %% Konfidenzintervall
   for kk = 1:length(perf.lda.accuracy)
       CI(kk).lda = ci(perf.lda.accuracy(:,kk));
       CI(kk).svm = ci(perf.svm.accuracy(:,kk));
       CI(kk).logreg = ci(perf.logreg.accuracy(:,kk));
   end
    
    if ~exist([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\'], 'dir')
        mkdir([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\'])
    end
    savefig([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\perf_without_pca_bootstrap_pseudotrials_trainfold_only.fig'])
   
    save([mainpath subjectdir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\perf_without_pca_bootstrap_pseudotrials_trainfold_only.mat'], 'perf')
    clear session
end




end

function [perf] = mvpa_leave_out_balldesign(session, balldesign)

%% Get default hyperparameters for the logreg and lda classifier
param_logreg = mv_get_hyperparameter('logreg');

param_lda = mv_get_hyperparameter('lda');

param_svm = mv_get_hyperparameter('svm');

    %% built train and test data based on CV folds
cf_logreg = cell(length(balldesign),length(session.time{1}));
cf_lda = cell(length(balldesign),length(session.time{1}));
cf_svm = cell(length(balldesign),length(session.time{1}));


%% Crossvalidation folds

CV = adi_crossval_leaveExemplarOut(session.balldesign);
data_trials = kh_trial2dat(session.trial);
%%

for kk = 1:length(balldesign)
    
    test_fold = data_trials(CV(kk).testset,:,:);
    clabel_test_fold = session.labels(CV(kk).testset);
    % hier überprüfen, ob testset eindeutig einer condition zugeordnet ist:
    % evtl. morgen weitermachen
%     if ismember
%         size(unique(clabel_test_fold),2) > 1
%         C = categorical(clabel_test_fold,[1 2 0],{'like','dislike','dontcare'});
%         if str
%         [group_count, b] = histcounts(C);
%         exist()
%     end
%     
    
    
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
    
    cfg = [];
    cfg.mode = 'bootstrap';
    cfg.averages = 4;
    cfg.repetitions = numel(trials_trainfold_like.trial);
    bootstrap_trainfold_like = fte_subaverage(cfg, trials_trainfold_like);

    cfg.repetitions = numel(trials_trainfold_dislike.trial);
    bootstrap_trainfold_dislike = fte_subaverage(cfg, trials_trainfold_dislike);
    
    
    data_bootstrap_trainfold_like = kh_trial2dat(bootstrap_trainfold_like.trial);
    data_bootstrap_trainfold_dislike = kh_trial2dat(bootstrap_trainfold_dislike.trial);
    
    data_bootstrap_trainfold = cat(1, data_bootstrap_trainfold_like, data_bootstrap_trainfold_dislike);
    clabel_train_fold = cat(2, ones(1, size(data_bootstrap_trainfold_like,1)), 2*ones(1, size(data_bootstrap_trainfold_dislike,1)));
    
%     clabel_train_fold(find(clabel_train_fold == 1))

    %% bootstrap testfold: geht nur bei eindeutigen Antworten, deshalb erst einmal weglassen:
%     [trls_testfold] = kh_trial2dat(test_fold);
%  
%     trials_testfold = [];
%     trials_testfold.trial = trls_testfold;
%     trials_testfold.time = session.time(1:length(trials_testfold.trial));
%     trials_testfold.fsample = session.fsample;
%     trials_testfold.dimord = 'rpt_chan_time';
%     trials_testfold.istimelock = 1;
%     trials_testfold.label = session.label;
%     
%     cfg = [];
%     cfg.mode = 'bootstrap';
%     cfg.averages = 4;
%     cfg.repetitions = numel(trials_testfold.trial);
%     bootstrap_testfold = fte_subaverage(cfg, trials_testfold);
%     data_bootstrap_testfold = kh_trial2dat(bootstrap_testfold.trial);
    
    
    %%  lda for single time points  
    for tt = 1:size(session.time{1,1},2)
        cf_lda{kk, tt} = train_lda(param_lda, data_bootstrap_trainfold(:,:,tt), clabel_train_fold);
        [predlabel, dval] = test_lda(cf_lda{kk, tt}, test_fold(:,:,tt)); 
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.lda.auc(kk, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
        perf.lda.accuracy(kk, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    end     
    
        %%  svm for single time points  
    for tt = 1:size(session.time{1,1},2)   
        cf_svm{kk, tt} = train_svm(param_svm, data_bootstrap_trainfold(:,:,tt), clabel_train_fold);
        [predlabel, dval] = test_svm(cf_svm{kk,tt}, test_fold(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.svm.auc(kk, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
        perf.svm.accuracy(kk, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    end   
    
       %% logreg for single time points  
    for tt = 1:size(session.time{1,1},2)   
        cf_logreg{kk, tt} = train_logreg(param_logreg, data_bootstrap_trainfold(:,:,tt), clabel_train_fold);
        [predlabel, dval, prob] = test_logreg(cf_logreg{kk, tt}, test_fold(:,:,tt));
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean

        % Calculate AUC and Accuracy
        perf.logreg.auc(kk, tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
        perf.logreg.accuracy(kk, tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    end
        
        
        
        
    
    clear test_fold train_fold clabel_train_fold clabel_test_fold
    
end

end

function [CV] = adi_crossval_leaveExemplarOut(balldesign_short)

for kk=1:length(balldesign_short) 
    balldesign_array(kk)=balldesign_short{kk};
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



