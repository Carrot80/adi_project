function [] = adi_leave_out_exemplar_singlsubj(subjectdir, filename, balldesign, comp, neighbours)

%% searchlight für einzelne runs einiger Probanden, die wegen Kopfpositionsänderung in der Pause zwischen runs (z.B. Absacken oder nach Vorne kippen) ausgeschlossen wurden
for ii = 2:length(subjectdir)
    

    dir_data = dir([subjectdir{ii} '*.mat']);
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
    for kk = 1:length(session.response_label)
        switch session.response_label{kk}
            case 'like' % 'Volley'
                session.labels(kk) = 1;
            case 'Neu_Like' % 'Volley'
                session.labels(kk) = 1;
            case 'dislike' % 'Space'
                session.labels(kk) = 2;
            case 'Neu_Dislike' % 'Space'
                session.labels(kk) = 2;
        end
    end
   
        
  %% z-transform trials  

    [session] = adi_ztrans_sensorspace(session);
    

  
    %% MVPA  
   [perf] = mvpa_leave_out_balldesign(session, balldesign, comp, neighbours);
   
    if ~exist([subjectdir{ii} 'searchlight\'], 'dir')
        mkdir([subjectdir{ii} '\searchlight\'])
    end
%     perf.config.methods = 'bootstrap 2*trialnumber';
%     perf.config.methods = 'bootstrap_pseudotrials';
    perf.number_of_trials = numel(session.trial);
    perf.config.methods = 'perf_without_pca_bootstrap_4pseudotrials_per_condition_trainfold_only_x2';
    perf.config.comment = ' da alle trials mit den uneindeutigen Antworten zur Klassifikation genutzt worden, wurden nur trials aus dem trainfold zur Bildung von Pseudotrials gemittelt; es wurde trialanzahl für like und dislike angeglichen ';
    perf.component = [num2str(comp(1)) '_' num2str(comp(2)) '_ms'];
    perf.features = session.label;
    
    
   %% mean accuracy
   for kk = 1:length(perf.lda.accuracy)
       perf.lda.mean_accuracy(kk) = mean(perf.lda.accuracy(:,kk));
       perf.svm.mean_accuracy(kk,1) = mean(perf.svm.accuracy(:,kk));
       perf.logreg.mean_accuracy(kk,1) = mean(perf.logreg.accuracy(:,kk));
   end
   
   %% figure result lda:
   
   
   stat = [];
   stat.time = comp(1);
   stat.mvpa.perf = perf.lda.mean_accuracy;
   stat.accuracy = perf.lda.mean_accuracy;
   stat.label = session.label;
   stat.dimord = 'time';

   cfg              = [];
   cfg.parameter    = 'accuracy';
   cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
   cfg.colorbar     = 'yes';
   cfg.title = 'Searchlight LDA'; 
   figure; ft_topoplotER(cfg, stat); 
   savefig([subjectdir{ii} num2str(comp(1)) '_' num2str(comp(2)) '_ms_lda.fig'])
   close 
   
   %% figure result svm:
   
   stat = [];
   stat.time = comp(1);
   stat.mvpa.perf = perf.svm.mean_accuracy';
   stat.accuracy = perf.svm.mean_accuracy';
   stat.label = session.label;
   stat.dimord = 'time';
   cfg              = [];
   cfg.parameter    = 'accuracy';
   cfg.layout       = '4D248_helmet.mat';  
   cfg.title = 'Searchlight SVM';
%    cfg.xlim         = [0, 0];
   cfg.colorbar     = 'yes';
   try
        ft_topoplotER(cfg, stat);
   catch
       
   end
   savefig([subjectdir{ii}  num2str(comp(1)) '_' num2str(comp(2)) '_ms_svm.fig'])

   
    %% figure result logreg:
   
   stat = [];
   stat.time = comp(1);
   stat.mvpa.perf = perf.logreg.mean_accuracy';
   stat.accuracy = perf.logreg.mean_accuracy';
   stat.label = session.label;
   stat.dimord = 'time';

   cfg              = [];
   cfg.parameter    = 'accuracy';
   cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
   cfg.colorbar     = 'yes';
   cfg.title = 'Searchlight logreg';
   try
        ft_topoplotER(cfg, stat);
   catch
       
   end
   savefig([subjectdir{ii} num2str(comp(1)) '_' num2str(comp(2)) '_ms_logreg.fig'])

   
 
  
    %% Konfidenzintervall
%    for kk = 1:length(perf.lda.accuracy)
%        perf.lda.ci_acc(kk) = ci(perf.lda.accuracy(:,kk));
%        perf.svm.ci_acc(kk) = ci(perf.svm.accuracy(:,kk));
%        perf.logreg.ci_acc(kk) = ci(perf.logreg.accuracy(:,kk));
%    end
    
    save([subjectdir{ii} num2str(comp(1)) '_' num2str(comp(2)) '_ms.mat'], 'perf')
    clear session
end




end

function [perf] = mvpa_leave_out_balldesign(session, balldesign, comp, neighbours)

%% Get default hyperparameters for the logreg and lda classifier
param_logreg = mv_get_hyperparameter('logreg');

param_lda = mv_get_hyperparameter('lda');

param_svm = mv_get_hyperparameter('svm');

    %% built train and test data based on CV folds
cf_logreg = cell(length(balldesign), length(session.label));
cf_lda = cell(length(balldesign), length(session.label));
cf_svm = cell(length(balldesign), length(session.label));


%% Crossvalidation folds

CV = adi_crossval_leaveExemplarOut(session.balldesign);
data_trials = kh_trial2dat(session.trial);
%%

for kk = 1:length(CV)
    
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
    
    max_numel_trials = max([numel(trials_trainfold_like.trial), numel(trials_trainfold_dislike.trial)]);
    
    cfg = [];
    cfg.mode = 'bootstrap';
    cfg.averages = 4;
    cfg.repetitions = max_numel_trials*2;
    bootstrap_trainfold_like = fte_subaverage(cfg, trials_trainfold_like);
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
%%

% We want to classify on time window of components 
time_idx = find(session.time{1} >= comp(1)  &  session.time{1} <= comp(2));
% nachbarschaft 5 cm laut Stefan bzw. noch besser template von Fieldtrip verwenden:

X_train = mean(data_bootstrap_trainfold(:,:,time_idx),3);
X_test = mean(test_fold(:,:,time_idx),3);

 %% LDA - Loop across features
for ff = 1:length(session.label)
    
    current_sensor = session.label(ff);
    current_neighb = neighbours(str2num(current_sensor{1}(2:end))).neighblabel;
    for pp = 1:length(current_neighb)
        neighb_index(pp) = find(strcmp(session.label, current_neighb(pp)));
    end
    nb = [ff neighb_index];

    fprintf('Classifying using feature %d with neighbours %s\n', ff, mat2str(setdiff(nb,ff)))
 
    Xfeat_train = reshape(X_train(:,nb,:), size(X_train,1), []);
    Xfeat_test = reshape(X_test(:,nb,:), size(X_test,1), []);
    
    %% LDA: 
    % Store the current state of the random number generator
    rng_state = rng;
    cf_lda{kk, ff} = train_lda(param_lda, Xfeat_train, clabel_train_fold);
    [predlabel, dval] = test_lda(cf_lda{kk, ff}, Xfeat_test); 
    % To calculate classification accuracy, compare the predicted labels to
    % the true labels and take the mean
    % Calculate AUC and Accuracy
    perf.lda.auc(kk, ff) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
    perf.lda.accuracy(kk, ff) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    
    %% SVM:
    % Store the current state of the random number generator
    rng_state = rng;
    cf_svm{kk, ff} = train_svm(param_svm, Xfeat_train, clabel_train_fold);
    [predlabel, dval] = test_svm(cf_svm{kk,ff}, Xfeat_test);
    % To calculate classification accuracy, compare the predicted labels to
    % the true labels and take the mean

    % Calculate AUC and Accuracy
    perf.svm.auc(kk, ff) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
    perf.svm.accuracy(kk, ff) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    
    %% logreg
    % Store the current state of the random number generator
        rng_state = rng;
     cf_logreg{kk, ff} = train_logreg(param_logreg, Xfeat_train, clabel_train_fold);
    [predlabel, dval, prob] = test_logreg(cf_logreg{kk, ff}, Xfeat_test);
    % To calculate classification accuracy, compare the predicted labels to
    % the true labels and take the mean

    % Calculate AUC and Accuracy
    perf.logreg.auc(kk, ff) = mv_calculate_performance('auc', 'dval', dval, clabel_test_fold);
    perf.logreg.accuracy(kk, ff) = mv_calculate_performance('acc', 'dval', dval, clabel_test_fold);
    
end

    clear Xfeat_test Xfeat_train clabel_train_fold clabel_test_fold neighb_index nb
    
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
    cfg.demean = 'yes';
    cfg.baselinewindow  = [-0.5 -0.030];
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

if ~isfield(data_bpfreq_res_sel, 'additional_cleaning')
        data_bpfreq_res_sel = setfield(data_bpfreq_res_sel, 'additional_cleaning', 'no');
end

right_fieldorder = {'label'; 'time';'trial'; 'ChannelFlag_Bst'; 'additional_cleaning'; 'cfg'; 'dimord'; 'fsample'; 'grad'; 'trialinfo'};
data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, right_fieldorder);

clearvars filename data_bpfreq data_bpfreq_res 

end



