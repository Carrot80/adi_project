function stats(like_training, dislike_training, like_cv, dislike_cv, time, path)

% evtl. lambda auf null setzen, so dass wir bei Searchlight nur Daten
% beschreiben?
% evtl. cfg.normalize = 'no'?
% bei sortierung der labels und trials leicht andere Accuracy-Werte; unklar
% weshalb, da Objecte verwendet werden
% was sagt mir die area under curve (sie ist gegenüber Klassenimbalancen weniger fehleranfällig)? Wo sehe ich, welcher lambda-wert
% verwendet wurde? 
% relevante Zeitbereiche: aus rms erstellen oder aus lda?
% 20 % der Daten aufheben, um Modell nochmal an neuen Daten zu testen?

% normalisierung ausschalten

%% speichere später result mit größeninfo


%%
% zuerst ohne regularisierung:  

if isempty(time)
    fprintf('no time input, computing classification across time...')

end
rng default

data_cv = [like_cv.trial dislike_cv.trial];

cfg = [] ;

for p = 1:length(data_cv)
    temp = data_cv{1,p};
    cfg.cv_data(p,:,:) = temp;
    clear temp
end

cfg.cv_labels_cond = [like_cv.labels_cond dislike_cv.labels_cond];

cfg.method          = 'mvpa_kh';
cfg.mvpa.classifier = 'logreg_kh'; %  multi-class Linear Discriminant Analysis (LDA)
cfg.mvpa.metric     = {'accuracy'; 'auc'};
cfg.mvpa.cv = 'testset'; %
% cfg.mvpa.param.bias = 1;
cfg.mvpa.param.lambda = 'auto';
cfg.mvpa.param.reg = 'l2';%'l2';
cfg.mvpa.param.plot = 0;
% cfg.mvpa.k          = 5;
% cfg.mvpa.normalise  = 'none'; % accuracy mit normlisierung etwas besser
% cfg.mvpa.balance = 'undersample'; % undersampling?
% cfg.latency         = [-0.5 1];
if ~isempty(time)
   cfg.latency     = time; 
   cfg.avgovertime = 'yes';
else
    cfg.latency     = [-0.5 1]; 
    cfg.avgovertime     = 'no';
end
cfg.design          = [ones(1, size(like_training.trial,2)) 2*ones(1,size(dislike_training.trial,2))]';
stat = ft_timelockstatistics_kh(cfg, like_training, dislike_training);

fprintf('Classification accuracy: %0.2f\n', stat.mvpa.perf.accuracy(:).acc)
number_of_subjects = size(unique(like_training.subject),2);
number_of_trials = size(like_training.trial,2);


figure
plot(like_training.time{1}, cat(1, stat.mvpa.perf.accuracy(:).acc))
fig_title = [num2str(number_of_subjects) ' subjects ' num2str(number_of_trials) ' trials'];
title(fig_title);
savefig([path fig_title '.fig'])
save ([path 'result_' num2str(number_of_subjects) ' subjects' num2str(number_of_trials) 'trials'], 'stat')
close all

end