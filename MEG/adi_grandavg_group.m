function  [grandavg] = grandavg_sensorspace(subjectpath, folderpath)

    
for ii = 1:length(subjectpath)
        %% like
        load([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_like.mat']);
        avg_like.subject =  subjectpath(ii).name;
        grandavg_like(1,ii).avg = avg_like;
        clear avg_like
        
        %% dislike
        load([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_dislike.mat']);
        avg_dislike.subject =  subjectpath(ii).name;
        grandavg_dislike(1,ii).avg = avg_dislike;
        clear avg_dislike
        
        % dontcare
        if exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_dontcare.mat'])
            load([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_dontcare.mat']);
            avg_dontcare.subject =  subjectpath(ii).name;
            grandavg_dontcare(1,ii).avg = avg_dontcare;
            clear avg_dontcare
        end
         
end

%% 

%% grandavg: 
          
cfg = [];
cfg.method  = 'within';
grandavg_like_group = ft_timelockgrandaverage(cfg, grandavg_like.avg);

grandavg_dislike_group = ft_timelockgrandaverage(cfg, grandavg_dislike.avg);
grandavg_dontcare_group = ft_timelockgrandaverage(cfg, grandavg_dontcare.avg);     

cfg = [];
cfg.method  = 'within';
grandavg = ft_timelockgrandaverage(cfg, grandavg_like_group, grandavg_dislike_group) 


save (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_like.mat' ], 'grandavg_like_group');
save (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_dislike.mat' ], 'grandavg_dislike_group');
save (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_dontcare.mat' ], 'grandavg_dontcare_group');

%% Erstellung Abbildung für Publikation:

cfg = [];
cfg.method = 'power';
[gmf_like] = ft_globalmeanfield(cfg, grandavg_like_group);
[gmf_dislike] = ft_globalmeanfield(cfg, grandavg_dislike_group);
[gmf_dontcare] = ft_globalmeanfield(cfg, grandavg_dontcare_group);


cfg = [];
cfg.method = 'power';
[gmf_like] = rms(cfg, grandavg_like_group);
[gmf_dislike] = rms(cfg, grandavg_dislike_group);
[gmf_dontcare] = rms(cfg, grandavg_dontcare_group);

figure
plot(gmf_like.time, gmf_like.avg, 'r', 'LineWidth',1)
axis tight
hold on
plot(gmf_dislike.time, gmf_dislike.avg, 'b', 'LineWidth',1)
hold on
plot(gmf_dontcare.time, gmf_dontcare.avg, 'k:', 'LineWidth',1)
legend({'like'; 'dislike'; 'dontcare'})
title('global field power ')
set(gcf,'color','w');
box off
legend({'like';  'dislike'; 'dontcare'}, 'boxoff')
legend('boxoff') 

figure
plot(gmf.time, gmf.avg, 'k', 'LineWidth',1.2)
hold on
plot(gmf.time, gmf_like.avg, 'r--', 'LineWidth',1.2)
hold on
plot(gmf.time, gmf_dislike.avg, 'b:', 'LineWidth',1.2)
axis tight
set(gcf,'color','w');
box off
legend({'all trials like and dislike'; 'like';  'dislike'}, 'boxoff')
legend('boxoff') 
       
 
%% topoplot für Publikation (like, dislike zusammen)


% baseline:
cfg = [];
cfg.xlim  = [-0.5 -0.03];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg);colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 1:
cfg = [];
cfg.xlim  = [0.03 0.15];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 2:

cfg = [];
cfg.xlim  = [0.15 0.56];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 3:
cfg = [];
cfg.xlim  = [0.56 0.74];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

%%


cfg = [];
cfg.xlim  = [0.0 0.1];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.1 0.2];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.2 0.3];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.3 0.4];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.4 0.5];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.5 0.6];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

cfg = [];
cfg.xlim  = [0.6 0.7];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );


%% components condition like
close all
% baseline:
cfg = [];
cfg.xlim  = [-0.5 -0.03];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_like_group);colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 1:
cfg = [];
cfg.xlim  = [0.03 0.15];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_like_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 2:

cfg = [];
cfg.xlim  = [0.15 0.56];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_like_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 3:
cfg = [];
cfg.xlim  = [0.56 0.74];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_like_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

%% components condition dislike
close all
% baseline:
cfg = [];
cfg.xlim  = [-0.5 -0.03];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_dislike_group);colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 1:
cfg = [];
cfg.xlim  = [0.03 0.15];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_dislike_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 2:

cfg = [];
cfg.xlim  = [0.15 0.56];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_dislike_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

% comp 3:
cfg = [];
cfg.xlim  = [0.56 0.74];
cfg.zlim = [-0.4 0.4];
figure; ft_topoplotER(cfg, grandavg_dislike_group); colorbar; set(gcf,'color','w'); %set(gcf, 'Position', [885   539   363   288] );

%%

figure
plot(grandavg_like_group.time, grandavg_like_group.avg)
axis tight
title('grandavg like')
figure
plot(grandavg_dislike_group.time, grandavg_dislike_group.avg)
axis tight
title('grandavg dislike')
        
savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_like.fig' ]);
        
%%
% comp1:

cfg = [];
cfg.channel   = 'all';
cfg.latency   = [.03 .15];
cfg.parameter = 'avg';
GA_like_comp1         = ft_timelockgrandaverage(cfg,grandavg_like.avg);
GA_dislike_comp1       = ft_timelockgrandaverage(cfg,grandavg_dislike.avg);

cfg = [];
cfg.showlabels  = 'yes';
cfg.layout      = '4D248.lay';
figure; ft_multiplotER(cfg,GA_like_comp1, GA_dislike_comp1)

% define the parameters for the statistical comparison
cfg = [];
cfg.channel     = 'all';
cfg.latency     = [0.56 0.74];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'FDR';

Nsub = 30;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike.avg);   % don't forget the {:}!


%%
cfg = [];
cfg.channel     = 'all';
% cfg.latency     = [0.03 0.15];
% cfg.latency     = [0.15 0.56];
cfg.latency     = [0.56 0.74];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;

Nsub = 30;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike.avg); 

cfg = [];
cfg.style     = 'blank';
cfg.layout    = '4D248.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_like_comp1)
title('Nonparametric: significant without multiple comparison correction')
   

%% grandavg dislike und dontcare zusammengenommen:
cfg = [];
cfg.method  = 'within';
grandavg_dislike_dontcare_group = ft_timelockgrandaverage(cfg, grandavg_dislike_group, grandavg_dontcare_group)

cfg = [];
cfg.channel     = 'all';
% cfg.latency     = [0.56 0.74];
cfg.latency     = [0.03 0.16];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'FDR';

Nsub = 30;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike_dontcare_group);   % don't forget the {:}!



%%
        
        
chan = 52;
time = [0.3 0.7];

% find the time points for the effect of interest in the grand average data
timesel_FIC = find(GA_FIC.time >= time(1) & GA_FIC.time <= time(2));
timesel_FC = find(GA_FC.time >= time(1) & GA_FC.time <= time(2));

% select the individual subject data from the time points and calculate the mean
for isub = 1:10
    values_FIC(isub)  = mean(allsubjFIC{isub}.avg(chan,timesel_FIC));
    values_FC(isub)  = mean(allsubjFC{isub}.avg(chan,timesel_FC));
end

% plot to see the effect in each subject
M = [values_FC',values_FIC'];
figure; plot(M','o-'); xlim([0.5 2.5])
legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
        'subj7', 'subj8', 'subj9', 'subj10'}, 'location','EastOutside');        
        
        
        
        
        
        
        
        
%%

         
         %% grandavg like and dislike together:
         
        cfg = [];
        cfg.method  = 'within';
        grandavg_group = ft_timelockgrandaverage(cfg, grandavg_like_group, grandavg_dislike_group);
        figure
        plot(grandavg_group.time, grandavg_group.avg)
        axis tight
        title('grandavg like, dislike')
        cfg = [];
        cfg.method = 'power';
        [gmf] = ft_globalmeanfield(cfg, grandavg_group);
        figure
        plot(gmf.time, gmf.avg,  'Color',[0.17, 0.17, 0.17], 'LineWidth',1)
        hold on
        
        cfg = [];
        cfg.method  = 'within';
        grandavg_group_alltrls = ft_timelockgrandaverage(cfg, grandavg_like_group, grandavg_dislike_group,grandavg_dontcare_group);
        cfg = [];
        cfg.method = 'power';
        [gmf_alltrls] = ft_globalmeanfield(cfg, grandavg_group_alltrls);
        
        figure
        plot(gmf_alltrls.time, gmf_alltrls.avg,  'k', 'LineWidth',1)
        set(gcf,'color','w');
        box off

        
        
 %% cluster permutation test:
 
cfg = [];
cfg.method  = 'within';
cfg.keepindividual = 'yes';
grandavg_like_group = ft_timelockgrandaverage(cfg, grandavg_like.avg);


cfg = [];
cfg.method  = 'within';
cfg.keepindividual = 'yes';
grandavg_dislike_group = ft_timelockgrandaverage(cfg, grandavg_dislike.avg);

cfg_neighb        = [];
cfg_neighb.method = 'distance';
neighbours        = ft_prepare_neighbours(cfg_neighb, grandavg_like_group);

cfg = [];
avglike = ft_timelockanalysis(cfg, grandavg_like_group);
avgdislike = ft_timelockanalysis(cfg, grandavg_dislike_group);
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectlike_vs_dislike = ft_math(cfg, avglike, avgdislike);

%% %% cluster permutation test:
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'depsamplesT'; % use the independent samples T-statistic as a measure to
                               % evaluate the effect at the sample level
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that
                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
                               % permutation distribution.
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
                               % required for a selected sample to be included
                               % in the clustering algorithm (default=0).
%cfg.neighbours = neighbours;   % see below
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 5000;      % number of draws from the permutation distribution

subj = size(grandavg_like_group.individual,1);
design = zeros(2,2*subj);
for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;
     
cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
                                 % which other sensors it can form clusters
% cfg.channel       = {'MEG'};     % cell-array with selected channel labels
cfg.latency       = [0.03 0.15];      % time interval over which the experimental
cfg.avgovertime = 'yes';                                 % conditions must be compared (in seconds)    
[stat] = ft_timelockstatistics(cfg, grandavg_like_group, grandavg_dislike_group);

% stat.posclusters(1)
% stat.negclusters(1)


% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
% if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.025; end; % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.
stat.cfg.alpha = 0.025; % da drei tests gerechnet werden (1 pro componente)


% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);


% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);

% pos = stat.posclusterslabelmat == 1; % or == 2, or 3, etc.
% neg = stat.negclusterslabelmat == 1;

% timestep = 0.01; % timestep between time windows for each subplot (in seconds)
% sampling_rate = 1017.25; % Data has a temporal resolution of 300 Hz
% sample_count = length(stat.time);
% number of temporal samples in the statistics object
% j = [0.03:timestep:0.15]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% m = [1:timestep*sampling_rate:sample_count+timestep*sampling_rate]; % t

[i1,i2] = match_str(raweffectlike_vs_dislike.label, stat.label);

% % plot
% for k = 1:length(j)-1;
%    subplot(3,round(length(j)/3),k);
%    cfg = [];
%    cfg.xlim=[j(k) j(k+1)];
%    pos_int = zeros(numel(raweffectlike_vs_dislike.label),1);
%    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
%    cfg.highlight = 'on';
%    cfg.highlightchannel = find(pos_int);
%    cfg.comment = 'xlim';
%    cfg.commentpos = 'title';
%    cfg.layout = '4D248.lay';
%    ft_topoplotER(cfg, raweffectlike_vs_dislike);
% end



% plot avg:
 
cfg = [];
cfg.xlim=[0.03 0.15]; 
cfg.zlim=[-0.15 0.15];
pos_int = zeros(numel(raweffectlike_vs_dislike.label),1);
pos_int(i1) = all(pos(i2, 1), 2);  
% neg_int = zeros(numel(raweffectlike_vs_dislike.label),1);
% neg_int(i1) = all(neg(i2, 1), 2); 
cfg.highlight = 'on';
cfg.highlightchannel = find(pos_int);
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = '4D248.lay';
figure;ft_topoplotER(cfg, raweffectlike_vs_dislike); set(gcf,'color','w'); colorbar
clear pos pos_int

end


