function  [grandavg] = grandavg_sensorspace(subjectpath, folderpath)

    
for ii = 1:length(subjectpath)
%         %% like
%         load([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_like.mat']);
%         avg_like.subject =  subjectpath(ii).name;
%         grandavg_like(1,ii).avg = avg_like;
%         clear avg_like
        
        %% dislike
        load([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_dislike.mat']);
        avg_dislike.subject =  subjectpath(ii).name;
        grandavg_dislike(1,ii).avg = avg_dislike;
        clear avg_dislike
        
        %% dontcare
%         if exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_dontcare.mat'])
%             load([subjectpath(ii).folder filesep subjectpath(ii).name filesep  folderpath 'grandavg_dontcare.mat']);
%             avg_dontcare.subject =  subjectpath(ii).name;
%             grandavg_dontcare(1,ii).avg = avg_dontcare;
%             clear avg_dontcare
%         end
         
end
       
   
         %% grandavg like 
          
        cfg = [];
        cfg.method  = 'within';
        grandavg_like_group = ft_timelockgrandaverage(cfg, grandavg_like.avg);
        grandavg_dislike_group = ft_timelockgrandaverage(cfg, grandavg_dislike.avg);
        grandavg_dontcare_group = ft_timelockgrandaverage(cfg, grandavg_dontcare.avg);     
        
        if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\'], 'dir')
            mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\'])
        end
        save (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_like.mat' ], 'grandavg_like_group');
        save (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_dislike.mat' ], 'grandavg_dislike_group');
        %%
        cfg = [];
        cfg.method = 'power';
        [gmf_like] = ft_globalmeanfield(cfg, grandavg_like_group);
        [gmf_dislike] = ft_globalmeanfield(cfg, grandavg_dislike_group);
         
        figure
        plot(gmf_like.time, gmf_like.avg)
        axis tight
        hold on
        plot(gmf_dislike.time, gmf_dislike.avg)
        legend({'like'; 'dislike'})
        title('global field power ')
         
         %% grandavg like and dislike together:
         
        cfg = [];
        cfg.method  = 'within';
        grandavg_group = ft_timelockgrandaverage(cfg, grandavg_like_group, grandavg_dislike_group);
        figure
        plot(grandavg_group.time, grandavg_group.avg)
        axis tight
        title('grandavg like and dislike')
        cfg = [];
        cfg.method = 'power';
        [gmf] = ft_globalmeanfield(cfg, grandavg_group);
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
       
% plot(x,y1,'g',x,y2,'b--o',x,y3,'c*')
        
        figure
        plot(grandavg_like_group.time, grandavg_like_group.avg)
        axis tight
        title('grandavg like')
        figure
        plot(grandavg_dislike_group.time, grandavg_dislike_group.avg)
        axis tight
        title('grandavg dislike')
        
        savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_like.fig' ]);
        
    
        

end


