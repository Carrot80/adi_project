function stats(like, dislike, time)

if isempty(time)
    fprintf('no time input, computing classification across time...')

end
rng default

cfg_logreg = [] ;
cfg_logreg.method          = 'mvpa';
cfg_logreg.mvpa.classifier = 'logreg'; %  'logreg' oder'lda' => multi-class Linear Discriminant Analysis (LDA)
cfg_logreg.mvpa.metric     = {'accuracy'; 'auc'};
cfg_logreg.mvpa.param.lambda =  [0.0001 0.0006 0.0036 0.026 0.13 0.1 0.2 0.3 0.4 0.5 0.6 0.77 1 4.6416 27.82559 166.81 1000];
cfg_logreg.mvpa.param.reg = 'l2';%
cfg_logreg.mvpa.param.plot = 0;
cfg_logreg.mvpa.param.repeat = 5;
cfg_logreg.mvpa.repeat = 5;
cfg_logreg.mvpa.normalise = 'demean';
if ~isempty(time)
   cfg_logreg.latency     = time; 
   cfg_logreg.avgovertime = 'yes';
else
    cfg_logreg.latency     = [-0.5 1]; 
    cfg_logreg.avgovertime     = 'no';
end
cfg_logreg.design          = [ones(1, size(like.trial,2)) 2*ones(1,size(dislike.trial,2))]';
tic
    stat = ft_timelockstatistics(cfg_logreg, like, dislike);
toc

fprintf('Classification accuracy: %0.2f\n', stat.accuracy)

figure
plot(stat.time, stat.accuracy)

mv_plot_result(stat.mvpa, like.time{1})


cfg_logreg.timextime   = 'yes';
stat_logreg_timeXtime = ft_timelockstatistics(cfg_logreg, like, dislike);
figure
mv_plot_result(stat_logreg_timeXtime.mvpa)



%% lda
cfg_LDA = [] ;
cfg_LDA.method          = 'mvpa';
cfg_LDA.mvpa.classifier = 'lda'; %  'logreg' oder'lda' => multi-class Linear Discriminant Analysis (LDA)
cfg_LDA.mvpa.metric     = {'accuracy'; 'auc'};
% cfg_LDA.mvpa.param.lambda =  [0.0001 0.0006 0.0036 0.026 0.13 0.1 0.2 0.3 0.4 0.5 0.6 0.77 1 4.6416 27.82559 166.81 1000];
% cfg_LDA.mvpa.param.reg = 'shrink';%
cfg_LDA.mvpa.param.k = 5;
% cfg_LDA.mvpa.param.plot = 0;
% cfg_LDA.mvpa.param.repeat = 5;
% cfg_LDA.mvpa.repeat = 5;
cfg_LDA.cv              = 'kfold';
cfg_LDA.k               = 5;
cfg_LDA.repeat          = 5;
cfg_LDA.mvpa.normalise = 'demean'; % 'none' or 'demean' or 'zscore'
if ~isempty(time)
   cfg_LDA.latency     = time; 
   cfg_LDA.avgovertime = 'yes';
else
    cfg_LDA.latency     = [-0.5 1]; 
    cfg_LDA.avgovertime     = 'no';
end
cfg_LDA.design          = [ones(1, size(like.trial,2)) 2*ones(1,size(dislike.trial,2))]';
tic
    stat_LDA = ft_timelockstatistics(cfg_LDA, like, dislike);
toc

save('stat_LDA', 'stat_LDA', '-v7.3')

fprintf('Classification accuracy: %0.2f\n', stat_LDA.accuracy)

figure
plot(stat_LDA.time, stat_LDA.accuracy)

mv_plot_result(stat_LDA.mvpa, like.time{1})



%% Plot classification accuracy across time
close all
mv_plot_result({stat_LDA.mvpa, stat_logreg.mvpa}, like.time{1}) % second argument is optional

%% time generalization:


cfg_LDA.timextime   = 'yes';
stat_LDA_timeXtime = ft_timelockstatistics(cfg_LDA, like, dislike);
figure
mv_plot_result(stat_LDA_timeXtime.mvpa)


figure
cfg_plot= [];
cfg_plot.x   = dat.time;
cfg_plot.y   = cfg_plot.x;
mv_plot_2D(cfg_plot, acc);
colormap jet
title('Accuracy')

figure
mv_plot_2D(cfg_plot, auc);
colormap jet
title('AUC')


end









% %% searchlight analyse: where
% 
% % comp1: Zeitintervall extrahiert aus rms
% 
% rng default
% 
% cfg = [];
% cfg.method      = 'triangulation';
% cfg.layout      = '4D248.lay';    
% cfg.channel     = dislike_training.label;
% neighbours = ft_prepare_neighbours(cfg);
% 
% cfg = [] ;  
% cfg.method      = 'mvpa';
% cfg.searchlight = 'yes';
% cfg.mvpa.neighbours  = neighbours;
% cfg.mvpa.classifier  = 'logreg'; % logreg, lda, svm, ensemble, kernel_fda (for nonlinear data), libsvm, liblinear 
% cfg.mvpa.param.reg = 'l2';
% cfg.mvpa.param.lambda = 0; %'auto'; % [0 0.01 0.03 0.1 0.3 1];
% cfg.mvpa.param.plot = 1;
% cfg.mvpa.correct_bias = 0;
% cfg.mvpa.bias = 0;
% cfg.mvpa.param.bias = 0;
% cfg.mvpa.metric      = [];
% cfg.mvpa.cv = 'none';
% cfg.design          = [ones(1, size(like_training.trial,2)) 2*ones(1,size(dislike_training.trial,2))]';
% cfg.latency     = [0.05 0.16];
% cfg.avgovertime = 'yes';
% stat_comp1 = ft_timelockstatistics(cfg, like_training, dislike_training);
% 
% cfg_fig              = [];
% cfg_fig.parameter    = 'auc'; %oder accuracy
% cfg_fig.layout       = '4D248.lay';            
% cfg_fig.xlim         = [0, 0];
% cfg_fig.colorbar     = 'yes';
% ft_topoplotER(cfg_fig, stat_comp1);
% title('comp1')
% 
% % comp2:
% cfg.latency     = [0.15 0.5];
% stat_comp2 = ft_timelockstatistics(cfg, session_sorted);
% ft_topoplotER(cfg_fig, stat_comp2);
% 
% % comp3:
% cfg.latency     = [0.5 0.8];
% stat_comp3 = ft_timelockstatistics(cfg, session_sorted);
% ft_topoplotER(cfg_fig, stat_comp3);
% 
% 
% 
% 
% 
% %
% 
% 
% 
% 
% 
% %%
% clusterbased_perm(sensordata_all_subj, time)
% % statistics:
% cfg=[];
% for p=1:length(sensordata_all_subj)
%     avg=ft_timelockanalysis(cfg, sensordata_all_subj(p));
%     sensordata_all_subj(p).avg = avg.avg;
%     clear avg
% end
% 
% for p=1:length(sensordata_all_subj)
%     n=1;
%     switch sensordata_all_subj(p).response_label{1}
%         case 'Neu_Like'
%             sensordata_all_subj(p).condition = 'like';
%         case 'like'
%             sensordata_all_subj(p).condition = 'like';
%          case 'Neu_Dislike'
%             sensordata_all_subj(p).condition = 'dislike';
%         case 'dislike'
%             sensordata_all_subj(p).condition = 'dislike';           
%     end
% %     like_single_subj = 
%     n=n+1;
% end
% 
% 
% 
% % like_avg = cell('trial', [], 'time', [], 'response_label',  [], 'balldesign_short',  [], 'label',  [], 'fsample',  [], 'grad',  [], 'cfg',  [], 'condition',  [], 'avg',  [], 'dimord', []);
% % dislike_avg = cell('trial', [], 'time', [], 'response_label',  [], 'balldesign_short',  [], 'label',  [], 'fsample',  [], 'grad',  [], 'cfg',  [], 'condition',  [], 'avg',  [], 'dimord', []);
% 
%   n=1;
%   o=1;
% for p=1:length(sensordata_all_subj)
%   
%     if 1==strcmp(sensordata_all_subj(p).condition, 'like')
%         sensordata_all_subj(p).dimord = 'chan_time';
%         like_avg{n} = sensordata_all_subj(p);
%         like_avg{n}.time = like_avg{n}.time{1,1};             
%         like_avg{n} = removefields(like_avg{n}, 'trial');    
%         n=n+1;
%     elseif 1==strcmp(sensordata_all_subj(p).condition, 'dislike')
%         sensordata_all_subj(p).dimord = 'chan_time';
%         dislike_avg{o} = sensordata_all_subj(p);
%         dislike_avg{o}.time = dislike_avg{o}.time{1,1};             
%         dislike_avg{o} = removefields(dislike_avg{o}, 'trial');   
%         o=o+1;
%     end
% end
% 
% 
% cfg = [];
% cfg.channel   = 'all';
% cfg.latency   = 'all';
% cfg.parameter = 'avg';
% cfg.keepindividual='yes';
% grandavg_like=ft_timelockgrandaverage(cfg, like_avg{:});
% grandavg_dislike=ft_timelockgrandaverage(cfg, dislike_avg{:});
% 
% 
% 
% 
% cfg = [];
% cfg.showlabels  = 'yes';
% % cfg.channel = sign_channels2;
% cfg.layout    	= 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\4D248.lay';
% figure; ft_multiplotER(cfg,grandavg_like, grandavg_dislike)
% 
% 
% 
% 
% % cfg = [];
% % cfg.channel = 'MLT12';
% % figure; ft_singleplotER(cfg,GA_FC, GA_FIC)
% 
% 
% 
% %dependent samples ttest
% FCminFIC = values_FC - values_FIC;
% [h,p,ci,stats] = ttest(FCminFIC, 0, 0.05) % H0: mean = 0, alpha 0.05
% 
% 
% % define the parameters for the statistical comparison
% cfg = [];
% % cfg.channel     = sign_channels; %'MEG';
% cfg.latency     = [0.7 0.8]; % comp 2: 0.57 0.7s; comp1: 0.0 0.220
% cfg.avgovertime = 'yes';
% cfg.parameter   = 'avg';
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05;
% cfg.correctm    = 'no';
% cfg.correcttail = 'prob';
% cfg.numrandomization = 1000;
% 
% 
% Nsub = 11;
% cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
% 
% stat = ft_timelockstatistics(cfg,like_avg{:},dislike_avg{:});   % don't forget the {:}!
% 
% ind_prob = find(stat.prob<0.05);
% sign_channels = stat.label(ind_prob);
% 
% load ('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20181213\template\layout\4D248_helmet.mat');
% 
% figure
% scatter(lay.pos(1:248,1), lay.pos(1:248,2))
% hold on
% ind=zeros(1, length(sign_channels));
% for i =1:length(sign_channels)
%     
%     temp = find(strcmp(lay.label, sign_channels{i}));
%     ind(i)=temp;
%     clear temp;
%     
% end
% 
% scatter(lay.pos(ind,1), lay.pos(ind,2), 'filled')
% title('sign. channels like vs. dislike 0.57-0.8s p<0.05')
% 
% squeeze(rms(grandavg_dislike.individual));
% figure
% boxplot(squeeze(rms(grandavg_like.individual(:,275:309))))
% figure
% boxplot(squeeze(rms(grandavg_dislike.individual(:,275:309))))
% 
% x1=squeeze(rms(grandavg_like.individual(:,275:309)));
% x2=squeeze(rms(grandavg_dislike.individual(:,275:309)));
% figure
% boxplot([x1',x2'],'Notch','on','Labels',{'like','dislike'},'Whisker',1)
% title('rms sensor data 0.57 bis 0.7s')
% 
% 
% boxplot([squeeze(squeeze(rms(grandavg_like.individual(:,275:309)))), squeeze(squeeze(rms(grandavg_dislike.individual(:,275:309))))],'Notch','on','Labels',{'like','dislike'})
% boxplot([squeeze(rms(grandavg_like.individual(:,275:309))), squeeze(rms(grandavg_dislike.individual(:,275:309)))],'Notch','on','Labels',{'like','dislike'})
% title('rms sensor data 0.57 bis 0.7s')
% [p,h,stats] = ranksum(squeeze(rms(grandavg_like.individual(:,275:309))),squeeze(rms(grandavg_dislike.individual(:,275:309))))
% [p,h,stats] = ranksum(squeeze(rms(grandavg_like.individual(:,275:309))),squeeze(rms(grandavg_dislike.individual(:,130:185))))
% 
% 
% 
% %%
% 
% time = [0.3 0.7];
% % Scaling of the vertical axis for the plots below
% ymax = 1.9e-13;
% isub = 12;
% 
% figure;
% for nsubj = 1:length(dislike_avg)
% subplot(3,4,nsubj)
%   % use the rectangle to indicate the time range used later
% %     rectangle('Position',[time(1) 0 (time(2)-time(1)) ymax],'FaceColor',[0.7 0.7 0.7]);    
%      hold on;
%     % plot the lines in front of the rectangle
%     plot(dislike_avg{nsubj}.time, mean(dislike_avg{nsubj}.avg(:,:)));
%     plot(like_avg{nsubj}.time,mean(like_avg{nsubj}.avg(:,:)),'r');
%     axis tight
% %     title(['subject ',num2str(nsubj)])
% %     ylim([0 1.9e-13])
% %     xlim([-1 2])
% end
% 
% axis off
% 
% %%
% 
% figure
% plot(grandavg_dislike.time. mean(grandavg_dislike.avg))
% 
% 
% end
% 
% function clusterbased_perm(session, time)
% 
% % das ganze mit einzelnen trials machen: 2 Bedingungen, like mit trials und
% % dislike mit trials: Dienstag
% ind_like = find(session.labels==1);
% like = [];
% like.trial=session.trial(ind_like);
% like.time=session.time(ind_like);
% like.grad=session.grad;
% like.label = session.label;
% like.dimord = 'chan_time';
% 
% ind_dislike = find(session.labels==2);
% dislike = [];
% dislike.trial=session.trial(ind_dislike);
% dislike.time=session.time(ind_dislike);
% dislike.grad=session.grad;
% dislike.label = session.label;
% dislike.dimord = 'chan_time';
% 
% cfg = [];
% cfg.correctm = 'cluster';
% cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
% cfg.statistic = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to
% 
% cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that
%                                % will be used for thresholding
% cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
%                                % permutation distribution.
% cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
%                                % required for a selected sample to be included
%                                % in the clustering algorithm (default=0).
% % cfg.neighbours = neighbours;   % see below
% cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
% cfg.clustertail = 0;
% cfg.alpha = 0.05;               % alpha level of the permutation test
% cfg.numrandomization = 100;      % number of draws from the permutation distribution
% 
% 
% design = zeros(1,size(like.trial,2) + size(dislike.trial,2));
% design(1,1:size(like.trial,2)) = 1;
% design(1,size(like.trial,2)+1:(size(like.trial,2)+size(dislike.trial,2))) = 2;
% 
% 
% cfg.design = design;             % design matrix
% cfg.ivar  = 1;                   % number or list with indices indicating the independent variable(s)
% 
% cfg_neighb        = [];
% cfg_neighb.method = 'distance';         
% neighbours        = ft_prepare_neighbours(cfg_neighb, session);
% 
% 
% cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
%                                % which other sensors it can form clusters
% cfg.channel       = {'MEG'};     % cell-array with selected channel labels
% cfg.latency       = time;       % time interval over which the experimental
%                                % conditions must be compared (in seconds)
% 
% [stat] = ft_timelockstatistics(cfg, like, dislike);
% 
% % Then take the difference of the averages using ft_math
% avg_dislike = grandavg_dislike;
% avg_dislike = rmfield(avg_dislike, 'individual');
% avg_dislike.label = grandavg_dislike.label;
% avg_like = grandavg_like;
% avg_like = rmfield(avg_like, 'individual');
% 
% avg_dislike.avg = squeeze(mean(grandavg_dislike.individual));
% avg_like.avg = squeeze(mean(grandavg_like.individual));
% avg_dislike.dimord = 'chan_time';
% avg_like.dimord = 'chan_time';
% 
% time1=nearest(avg_like.time, time(1));
% time2=nearest(avg_like.time, time(2));
% 
% avg_like_selected = avg_like;
% avg_like_selected.avg=avg_like.avg(:,time1:time2);
% avg_like_selected.time=avg_like.time(1,time1:time2);
% 
% avg_dislike_selected = avg_dislike;
% avg_dislike_selected.avg=avg_dislike.avg(:,time1:time2);
% avg_dislike_selected.time=avg_dislike.time(1,time1:time2);
% 
% cfg  = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% raweffectLike_vs_Dislike = ft_math(cfg,avg_like_selected,avg_dislike_selected);
% 
% 
% % Then, find which clusters are significant, outputting their indices as held in stat.posclusters
% % In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
% if ~isfield(stat.cfg,'alpha') 
%     stat.cfg.alpha = 0.025;
% 
% end % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.
% 
% for i=1:size(stat.posclusters,2)
%     poscluster_pval(i)= stat.posclusters(i).prob;
% end
% 
% pos_signif_clust = find(poscluster_pval < stat.cfg.alpha);
% % (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% % make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
% pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
% 
% % and now for the negative clusters...
% neg_cluster_pvals = [stat.negclusters(:).prob];
% neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
% neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
% 
% timestep = 0.005;	% timestep between time windows for each subplot (in seconds) 
% sampling_rate = session.fsample;	% Data has a temporal resolution of 300 Hz
% sample_count = length(stat.time); % number of temporal samples in the statistics object
% % j = [time(1):timestep:time(2)]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% j = [0.6:timestep:0.68]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% 
% m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in MEEG samples
% 
% % First ensure the channels to have the same order in the average and in the statistical output.
% % This might not be the case, because ft_math might shuffle the order  
% [i1,i2] = match_str(raweffectLike_vs_Dislike.label, stat.label);
% 
% 
% figure
% for k = 1:length(m)
%    subplot(4,round(length(m)/5),k);
%    cfg = [];   
%    cfg.xlim=[j(k) j(k+1)];   % time interval of the subplot
% %    cfg.zlim = [-2.5e-13 2.5e-13];
%    % If a channel reaches this significance, then
%    % the element of pos_int with an index equal to that channel
%    % number will be set to 1 (otherwise 0).
% 
%    % Next, check which channels are significant over the
%    % entire time interval of interest.
%    pos_int = zeros(numel(raweffectLike_vs_Dislike.label),1);
%    neg_int = zeros(numel(raweffectLike_vs_Dislike.label),1);
%    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
%    neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
% 
%    cfg.highlight = 'on';
%    % Get the index of each significant channel
%    cfg.highlightchannel = find(pos_int | neg_int);
%    cfg.comment = 'xlim';   
%    cfg.commentpos = 'title';   
% %    cfg.layout = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\4D248.lay';
%    cfg.layout = layout;
%    cfg.interactive = 'no';
%    ft_topoplotER(cfg, raweffectLike_vs_Dislike);   
% end
% 
% end
% 
% function x()
% 
% cfg=[];
% cfg.latency = [0.57 0.7]; 
% for i=1:length(dislike_avg)
%     data_selected = ft_selectdata(cfg, dislike_avg{i});
%     dislike_avg_selected{i} = data_selected;
%     clear data_selected
% end
% 
% for i=1:length(like_avg)
%     data_selected = ft_selectdata(cfg, like_avg{i})
%     like_avg_selected{i} = data_selected;
%     clear data_selected
% end
% 
% cfg=[];
% cfg.latency = [0.57 0.7]; 
% session_selected = ft_selectdata(cfg, session);
% session_sel = session;
% for i=1:length(session.trial)
%     data_selected = session.trial{i}(:, nearest(session.time{1}, 0.57): nearest(session.time{1}, 0.7));
%     session_sel.trial{i} = data_selected;
%     clear data_selected
% end
% for i=1:length(session.time)
%     data_selected = session.time{i}(1, nearest(session.time{1}, 0.57): nearest(session.time{1}, 0.7));
%     session_sel.time{i} = data_selected;
%     clear data_selected
% end
% 
% 
% end