function  [Amplitude] = grandavg_sensorspace(subjectpath, folderpath, comp1)

    
for ii = 1:length(subjectpath)
        %% like
        load([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_like.mat']);
        
        % z-score-normalization:
        [avg_like_zscore] = subf_zscore(avg_like);
        
        time_beg = nearest(avg_like_zscore.time, comp1(1));
        time_end = nearest(avg_like_zscore.time, comp1(2));
        Amplitude.like.amp(ii) =  mean(mean(avg_like_zscore.avg(:,time_beg:time_end)));
        Amplitude.like.rms(ii) = mean(rms(avg_like_zscore.avg(:,time_beg:time_end)));
       
        
    
        %% dislike
        load([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_dislike.mat']);
         % z-score-normalization:
        [avg_dislike_zscore] = subf_zscore(avg_dislike);
        
        time_beg = nearest(avg_dislike_zscore.time, comp1(1));
        time_end = nearest(avg_dislike_zscore.time, comp1(2));
        Amplitude.dislike.amp(ii) =  mean(mean(avg_dislike_zscore.avg(:,time_beg:time_end)));
        Amplitude.dislike.rms(ii) = mean(rms(avg_dislike_zscore.avg(:,time_beg:time_end)));
       
         
end


end

function  [output_data_zscore] = subf_zscore(input_data)

output_data_zscore = input_data;

output_data_zscore.avg = [];
output_data_zscore.avg = zeros(size(input_data.avg, 1), size(input_data.avg, 2));
for oo = 1:size(input_data.avg,1)   
    M = mean(input_data.avg(oo,1:find(abs(input_data.time) == min(abs(input_data.time)))));
    STD = std(input_data.avg(oo,1:find(abs(input_data.time) == min(abs(input_data.time)))));
    output_data_zscore.avg(oo,:) = (input_data.avg(oo,:)- M)/STD;
    clearvars M STD
end


clear input_data

end

% %%
% % comp1:
% 
% cfg = [];
% cfg.channel   = 'all';
% cfg.latency   = [.03 .15];
% cfg.parameter = 'avg';
% GA_like_comp1         = ft_timelockgrandaverage(cfg,grandavg_like.avg);
% GA_dislike_comp1       = ft_timelockgrandaverage(cfg,grandavg_dislike.avg);
% 
% cfg = [];
% cfg.showlabels  = 'yes';
% cfg.layout      = '4D248.lay';
% figure; ft_multiplotER(cfg,GA_like_comp1, GA_dislike_comp1)
% 
% % define the parameters for the statistical comparison
% cfg = [];
% cfg.channel     = 'all';
% cfg.latency     = [0.56 0.74];
% cfg.avgovertime = 'yes';
% cfg.parameter   = 'avg';
% cfg.method      = 'analytic';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05;
% cfg.correctm    = 'FDR';
% 
% Nsub = 30;
% cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
% 
% stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike.avg);   % don't forget the {:}!
% 
% 
% %%
% cfg = [];
% cfg.channel     = 'all';
% % cfg.latency     = [0.03 0.15];
% % cfg.latency     = [0.15 0.56];
% cfg.latency     = [0.56 0.74];
% cfg.avgovertime = 'yes';
% cfg.parameter   = 'avg';
% cfg.method      = 'montecarlo';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05;
% cfg.correctm    = 'no';
% cfg.correcttail = 'prob';
% cfg.numrandomization = 1000;
% 
% Nsub = 30;
% cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
% 
% stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike.avg); 
% 
% cfg = [];
% cfg.style     = 'blank';
% cfg.layout    = '4D248.lay';
% cfg.highlight = 'on';
% cfg.highlightchannel = find(stat.mask);
% cfg.comment   = 'no';
% figure; ft_topoplotER(cfg, GA_like_comp1)
% title('Nonparametric: significant without multiple comparison correction')
%    
%          %% grandavg like 
%           
% cfg = [];
% cfg.method  = 'within';
% grandavg_like_group = ft_timelockgrandaverage(cfg, grandavg_like.avg);
% grandavg_dislike_group = ft_timelockgrandaverage(cfg, grandavg_dislike.avg);
% grandavg_dontcare_group = ft_timelockgrandaverage(cfg, grandavg_dontcare.avg);     
% 
% save (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_like.mat' ], 'grandavg_like_group');
% save (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_dislike.mat' ], 'grandavg_dislike_group');
% save (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_dontcare.mat' ], 'grandavg_dontcare_group');
% 
% %% grandavg dislike und dontcare zusammengenommen:
% cfg = [];
% cfg.method  = 'within';
% grandavg_dislike_dontcare_group = ft_timelockgrandaverage(cfg, grandavg_dislike_group, grandavg_dontcare_group)
% 
% cfg = [];
% cfg.channel     = 'all';
% % cfg.latency     = [0.56 0.74];
% cfg.latency     = [0.03 0.16];
% cfg.avgovertime = 'yes';
% cfg.parameter   = 'avg';
% cfg.method      = 'analytic';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05;
% cfg.correctm    = 'FDR';
% 
% Nsub = 30;
% cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
% 
% stat = ft_timelockstatistics(cfg,grandavg_like.avg ,grandavg_dislike_dontcare_group);   % don't forget the {:}!
% 
% 
% 
% %%
%         
%         
% chan = 52;
% time = [0.3 0.7];
% 
% % find the time points for the effect of interest in the grand average data
% timesel_FIC = find(GA_FIC.time >= time(1) & GA_FIC.time <= time(2));
% timesel_FC = find(GA_FC.time >= time(1) & GA_FC.time <= time(2));
% 
% % select the individual subject data from the time points and calculate the mean
% for isub = 1:10
%     values_FIC(isub)  = mean(allsubjFIC{isub}.avg(chan,timesel_FIC));
%     values_FC(isub)  = mean(allsubjFC{isub}.avg(chan,timesel_FC));
% end
% 
% % plot to see the effect in each subject
% M = [values_FC',values_FIC'];
% figure; plot(M','o-'); xlim([0.5 2.5])
% legend({'subj1', 'subj2', 'subj3', 'subj4', 'subj5', 'subj6', ...
%         'subj7', 'subj8', 'subj9', 'subj10'}, 'location','EastOutside');        
%         
%         
%         
%         
%         
%         
%         
%         
% %%
% cfg = [];
% cfg.method = 'power';
% [gmf_like] = ft_globalmeanfield(cfg, grandavg_like_group);
% [gmf_dislike] = ft_globalmeanfield(cfg, grandavg_dislike_group);
% [gmf_dontcare] = ft_globalmeanfield(cfg, grandavg_dontcare_group);
% 
% figure
% plot(gmf_like.time, gmf_like.avg)
% axis tight
% hold on
% plot(gmf_dislike.time, gmf_dislike.avg)
% hold on
% plot(gmf_dontcare.time, gmf_dontcare.avg)
% legend({'like'; 'dislike'; 'dontcare'})
% title('global field power ')
%         
%          figure
%         plot(gmf.time, gmf.avg, 'k', 'LineWidth',1.2)
%         hold on
%         plot(gmf.time, gmf_like.avg, 'r--', 'LineWidth',1.2)
%         hold on
%         plot(gmf.time, gmf_dislike.avg, 'b:', 'LineWidth',1.2)
%         axis tight
%         set(gcf,'color','w');
%         box off
%         legend({'all trials like and dislike'; 'like';  'dislike'}, 'boxoff')
%         legend('boxoff') 
%        
%          
%          %% grandavg like and dislike together:
%          
%         cfg = [];
%         cfg.method  = 'within';
%         grandavg_group = ft_timelockgrandaverage(cfg, grandavg_like_group, grandavg_dislike_group);
%         figure
%         plot(grandavg_group.time, grandavg_group.avg)
%         axis tight
%         title('grandavg like, dislike')
%         cfg = [];
%         cfg.method = 'power';
%         [gmf] = ft_globalmeanfield(cfg, grandavg_group);
%         figure
%         plot(gmf.time, gmf.avg,  'Color',[0.17, 0.17, 0.17], 'LineWidth',1)
%         hold on
%         
%         cfg = [];
%         cfg.method  = 'within';
%         grandavg_group_alltrls = ft_timelockgrandaverage(cfg, grandavg_like_group, grandavg_dislike_group,grandavg_dontcare_group);
%         cfg = [];
%         cfg.method = 'power';
%         [gmf_alltrls] = ft_globalmeanfield(cfg, grandavg_group_alltrls);
%         
%         figure
%         plot(gmf_alltrls.time, gmf_alltrls.avg,  'k', 'LineWidth',1)
%         set(gcf,'color','w');
%         box off
% 
%         
%         
%         figure
% %         plot(gmf.time, gmf.avg,  'Color',[0.17, 0.17, 0.17], 'LineWidth',1)
% %         hold on
%         plot(gmf.time, gmf_like.avg, 'r', 'LineWidth',1)
%         hold on
%         plot(gmf.time, gmf_dislike.avg, 'b', 'LineWidth',1)
%         hold on
%         plot(gmf.time, gmf_dontcare.avg, 'k:', 'LineWidth',1)
%         axis tight
%         set(gcf,'color','w');
%         box off
%         legend({'like';  'dislike'; 'dontcare'}, 'boxoff')
%         legend('boxoff') 
%        
% % plot(x,y1,'g',x,y2,'b--o',x,y3,'c*')
%         
%         figure
%         plot(grandavg_like_group.time, grandavg_like_group.avg)
%         axis tight
%         title('grandavg like')
%         figure
%         plot(grandavg_dislike_group.time, grandavg_dislike_group.avg)
%         axis tight
%         title('grandavg dislike')
%         
%         savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_like.fig' ]);
%         
%     
%         
% 
% end
% 
% 
