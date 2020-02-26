function  [grandavg, trl_count] = grandavg_sensorspace(subjectpath, grandavg, delete_run, path2save, path2inputfile)

if 1 == isempty(delete_run)
%     counter_like = 1;
%     counter_dislike = 1;
%     counter_dontcare = 1;
    
    for ii = [10] %:length(subjectpath)
        %% like
        filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2inputfile 'Neu_Like*.mat']);
         counter = 1;
        for kk = 1:length(filename)
            load ([filename(kk).folder filesep filename(kk).name])
            [trials] = kh_trial2dat(cleanMEG_interp.trial);
            if any(any(any(isnan(trials))))
                error(['NaNs in ' filename(kk).name ' ' subjectpath(ii).name]);
            end
            clear trials
            cfg = [];
            cfg.lpfilter = 'yes';
            cfg.lpfreq   = 45;     
            cfg.demean = 'yes';
            cfg.baselinewindow  = [-0.5 -0.030];
            data = ft_preprocessing(cfg, cleanMEG_interp);
            
            %% select data
            
            cfg = [];
            cfg.latency = [-0.5 1];
            data = ft_selectdata(cfg, data);
            
            %% z-transformation:
             [data] = subfun_ztransform(data);
             
            %% timelockanalysis:
            cfg = [];
            avg = ft_timelockanalysis(cfg, data);
            avg.subject = subjectpath(ii).name; 
            avg.run = filename(kk).name(end-4);   
            avg.balldesign = cleanMEG_interp.trialinfo.balldesign_short;
            grandavg.like.(subjectpath(ii).name)(counter).avg = avg;
            trl_count.like(ii,counter) = length(cleanMEG_interp.trial);
            clear avg cleanMEG_interp data
            counter = counter + 1;
        end     
        counter = counter - 1;
        
        % grandavg like per subject 
        
        
        cfg = [];
        cfg.method  = 'within';
        switch counter
            case 3
                avg_like = ft_timelockgrandaverage(cfg, grandavg.like.(subjectpath(ii).name)(1).avg, grandavg.like.(subjectpath(ii).name)(2).avg, grandavg.like.(subjectpath(ii).name)(3).avg);
            case 2
                avg_like = ft_timelockgrandaverage(cfg, grandavg.like.(subjectpath(ii).name)(1).avg, grandavg.like.(subjectpath(ii).name)(2).avg);
            case 1
                avg_like = grandavg.like.(subjectpath(ii).name)(counter).avg;
        end
        if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save], 'dir')
            mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save])
        end
        save ([subjectpath(ii).folder filesep subjectpath(ii).name filesep path2save 'grandavg_like.mat' ], 'avg_like');
        figure
        plot(avg_like.time, avg_like.avg)
        title('grandavg like')
        axis tight
        savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'grandavg_like.fig' ]);
        
        
        %% dislike:
         filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2inputfile 'Neu_Dislike*.mat']);
         counter = 1;
        for kk = 1:length(filename)
            load ([filename(kk).folder filesep filename(kk).name])
            [trials] = kh_trial2dat(cleanMEG_interp.trial);
            if any(any(any(isnan(trials))))
                error(['NaNs in ' filename(kk).name ' ' subjectpath(ii).name]);
            end
            clear trials
            cfg = [];
            cfg.lpfilter = 'yes';
            cfg.lpfreq   = 45;     
            cfg.demean = 'yes';
            cfg.baselinewindow  = [-0.5 -0.030];
            data = ft_preprocessing(cfg, cleanMEG_interp);
            
            %% select data
            cfg = [];
            cfg.latency = [-0.5 1];
            data = ft_selectdata(cfg, data);
            
            %% z-transformation:
             [data] = subfun_ztransform(data);
            
            %% timelockanalysis
            
            cfg = [];
            avg = ft_timelockanalysis(cfg, data);
            avg.subject = subjectpath(ii).name; 
            avg.run = filename(kk).name(end-4);   
            avg.balldesign = cleanMEG_interp.trialinfo.balldesign_short;
            grandavg.dislike.(subjectpath(ii).name)(counter).avg = avg; 
            trl_count.dislike(ii,counter) = length(cleanMEG_interp.trial);
            clear avg cleanMEG_interp data
            counter = counter + 1;
        end     
         % grandavg dislike per subject 
         counter = counter - 1;
        cfg = [];
        cfg.method  = 'within';
        switch counter
            case 3
                avg_dislike = ft_timelockgrandaverage(cfg, grandavg.dislike.(subjectpath(ii).name)(1).avg, grandavg.dislike.(subjectpath(ii).name)(2).avg, grandavg.dislike.(subjectpath(ii).name)(3).avg);
            case 2
                avg_dislike = ft_timelockgrandaverage(cfg, grandavg.dislike.(subjectpath(ii).name)(1).avg, grandavg.dislike.(subjectpath(ii).name)(2).avg);
            case 1
                avg_dislike = grandavg.dislike.(subjectpath(ii).name)(counter).avg;
        end
        if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save], 'dir')
            mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save])
        end
        save ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'grandavg_dislike.mat' ], 'avg_dislike');
        figure
        plot(avg_dislike.time, avg_dislike.avg)
        title('grandavg dislike')
        savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'grandavg_dislike.fig' ]);
       
         
        
          %% dontcare:
         filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2inputfile 'Neu_Dontcare*.mat']);
         counter = 1;
         if ~isempty(filename)
            for kk = 1:length(filename)
                load ([filename(kk).folder filesep filename(kk).name])
                [trials] = kh_trial2dat(cleanMEG_interp.trial);
                if any(any(any(isnan(trials))))
                    error(['NaNs in ' filename(kk).name ' ' subjectpath(ii).name]);
                end
                clear trials
                cfg = [];
                cfg.lpfilter = 'yes';
                cfg.lpfreq   = 45;     
                cfg.demean = 'yes';
                cfg.baselinewindow  = [-0.5 -0.030];
                data = ft_preprocessing(cfg, cleanMEG_interp);
                
                %% select data
                cfg = [];
                cfg.latency = [-0.5 1];
                data = ft_selectdata(cfg, data);

                %% z-transformation:
                 [data] = subfun_ztransform(data);
            
            %% timelockanalysis
                cfg = [];
                avg = ft_timelockanalysis(cfg, data);
                avg.subject = subjectpath(ii).name; 
                avg.run = filename(kk).name(end-4);   
                avg.balldesign = cleanMEG_interp.trialinfo.balldesign_short;
                grandavg.dontcare.(subjectpath(ii).name)(counter).avg = avg; 
                trl_count.dontcare(ii,counter) = length(cleanMEG_interp.trial);
                clear avg cleanMEG_interp data
                counter = counter + 1;
            end  
         
             counter = counter - 1;
               % grandavg dontcare per subject 
            cfg = [];
            cfg.method  = 'within';
            switch counter
                case 3
                    avg_dontcare = ft_timelockgrandaverage(cfg, grandavg.dontcare.(subjectpath(ii).name)(1).avg, grandavg.dontcare.(subjectpath(ii).name)(2).avg, grandavg.dontcare.(subjectpath(ii).name)(3).avg);
                case 2
                    avg_dontcare = ft_timelockgrandaverage(cfg, grandavg.dontcare.(subjectpath(ii).name)(1).avg, grandavg.dontcare.(subjectpath(ii).name)(2).avg);
                case 1
                    avg_dontcare = grandavg.dontcare.(subjectpath(ii).name)(counter).avg;
            end
            if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save], 'dir')
                mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'grandavg\'])
            end
            save ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'grandavg_dontcare.mat' ], 'avg_dontcare');
            figure
            plot(avg_dontcare.time, avg_dontcare.avg)
            title('grandavg dontcare')
            savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'grandavg_dontcare.fig' ]);
         end
        close all
        
       %% global field power: 
        cfg = [];
        cfg.method = 'power';
        [gmf_like] = ft_globalmeanfield(cfg, avg_like);
        [gmf_dislike] = ft_globalmeanfield(cfg, avg_dislike);
        if exist('avg_dontcare', 'var')
            [gmf_dontcare] = ft_globalmeanfield(cfg, avg_dontcare);
            figure('Renderer', 'painters', 'Position', [300 300 800 200])
            plot(gmf_like.time, gmf_like.avg, 'r')
            axis tight
            hold on
            plot(gmf_dislike.time, gmf_dislike.avg, 'b')
            hold on
            plot(gmf_dontcare.time, gmf_dontcare.avg, 'k:')
            legend({'like'; 'dislike'; 'dontcare'})
            legend('boxoff') 
            title('global field power ')
         
            savefig([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'global_field_power_avg.fig' ]);
            clear avg_dontcare
        else
            figure('Renderer', 'painters', 'Position', [300 300 800 200])
            plot(gmf_like.time, gmf_like.avg)
            axis tight
            hold on
            plot(gmf_dislike.time, gmf_dislike.avg)
            legend({'like'; 'dislike'})
            title('global field power ')
            legend('boxoff') 
            savefig([subjectpath(ii).folder filesep subjectpath(ii).name filesep  path2save 'global_field_power_avg.fig' ]);
            
        end
        close all
        clear avg_like avg_dislike
         
    end
        
else
    counter_like=1;
    counter_dislike=1;
    cfg = [];
    for i=1:length(subjectpath)
        filename = dir([subjectpath(i).folder filesep subjectpath(i).name filesep  path2inputfile 'Neu_Like*.mat']);

        for k=1:length(filename)
            load ([filename(k).folder filesep filename(k).name])
            % delete balldesigns which were not clearly rated:
            delete_trials = zeros(1, length(cleanMEG_interp.trial));
            for p=1:length(cleanMEG_interp.trial)
                trial=cleanMEG_interp.trialinfo.balldesign_short{1,p}{1,1};
                % check if trial should be deleted:
                if 1==delete_run.(subjectpath(i).name).(['run' filename(k).name(end-4)]).(trial)
                   delete_trials(p) = 1;
                end
            end

            cleanMEG_interp.trial(find(delete_trials))=[];
            cleanMEG_interp.time(find(delete_trials))=[];
            cleanMEG_interp.trialinfo.balldesign_short(find(delete_trials))=[];
            grandavg.like(counter_like).subject = subjectpath(i).name; 
            grandavg.like(counter_like).run = filename(k).name(end-4); 
            grandavg.like(counter_like).trials = cleanMEG_interp.trialinfo.balldesign_short; 

            if ~isempty(cleanMEG_interp.trial)
                cfg = [];
                cfg.lpfilter = 'yes';
                cfg.lpfreq   = 45;     
                cfg.demean = 'yes';
                cfg.baselinewindow  = [-0.5 0];
                data = ft_preprocessing(cfg, cleanMEG_interp);
                cfg = [];
                avg = ft_timelockanalysis(cfg, data);
                grandavg.like(counter_like).avg = avg.avg; 
                clear avg data
            else
                grandavg.like(counter_like).avg = [];
            end

            clear cleanMEG_interp delete_trials
            counter_like=counter_like+1;
        end

        filename = dir([subjectpath(i).folder filesep subjectpath(i).name filesep  'MEG_analysis\noisereduced\1_95Hz\02_interpolated\Neu_Dislike*.mat']);
        for k=1:length(filename)
            load ([filename(k).folder filesep filename(k).name])
            % delete balldesigns which were not clearly rated:
            delete_trials = zeros(1, length(cleanMEG_interp.trial));
            for p=1:length(cleanMEG_interp.trial)
                trial=cleanMEG_interp.trialinfo.balldesign_short{1,p}{1,1};
                % check if trial should be deleted:
                if 1==delete_run.(subjectpath(i).name).(['run' filename(k).name(end-4)]).(trial)
                   delete_trials(p) = 1;
                end
            end

            cleanMEG_interp.trial(find(delete_trials))=[];
            cleanMEG_interp.time(find(delete_trials))=[];
            grandavg.dislike(counter_dislike).subject = subjectpath(i).name; 
            grandavg.dislike(counter_dislike).run = filename(k).name(end-4);
            cleanMEG_interp.trialinfo.balldesign_short(find(delete_trials))=[];
            grandavg.dislike(counter_dislike).trials = cleanMEG_interp.trialinfo.balldesign_short; 
            if ~isempty(cleanMEG_interp.trial) 
                cfg = [];
                cfg.lpfilter = 'yes';
                cfg.lpfreq   = 45;     
                cfg.demean = 'yes';
                cfg.baselinewindow  = [-0.5 0];
                data = ft_preprocessing(cfg, cleanMEG_interp);
                cfg = [];
                avg = ft_timelockanalysis(cfg, data);
                grandavg.dislike(counter_dislike).avg = avg.avg;
                clear avg
            else
                grandavg.dislike(counter_dislike).avg = [];
            end

            counter_dislike=counter_dislike+1;
            clear cleanMEG_interp delete_trials
        end

    end


end





end


function [output_data_zscore] = subfun_ztransform(input_data)

baseline_zero_samples = nearest(input_data.time{1,1},-0.03);
baseline_dur_samples = nearest(input_data.time{1,1}, -0.5);

output_data_zscore = input_data;

for pp=1:length(input_data)
    output_data_zscore(pp).trial = [];
    output_data_zscore(pp).trial = cell(1, length(input_data(pp).trial));
    for kk = 1:length(input_data(pp).trial)
        output_data_zscore(pp).trial{kk} = zeros(size(input_data(pp).trial{kk},1),size(input_data(pp).trial{kk},2));
        for oo = 1:size(input_data(pp).trial{kk},1)   
            M = mean(input_data(pp).trial{kk}(oo, baseline_dur_samples:baseline_zero_samples));
            STD = std(input_data(pp).trial{kk}(oo,baseline_dur_samples:baseline_zero_samples));
            output_data_zscore(pp).trial{kk}(oo,:) = (input_data(pp).trial{kk}(oo,:)- M)/STD;
            clearvars M STD
        end
    end
end

clear input_data


end