function  [grandavg] = grandavg_sensorspace(subjectpath, grandavg, delete_run)

if 1 == isempty(delete_run)
%     counter_like = 1;
%     counter_dislike = 1;
%     counter_dontcare = 1;
    
    for ii = 1:length(subjectpath)
        %% like
        filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\02_interpolated\Neu_Like*.mat']);
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
            cfg.baselinewindow  = [-0.5 0];
            data = ft_preprocessing(cfg, cleanMEG_interp);
            cfg = [];
            cfg.latency = [-0.5 1];
            avg = ft_timelockanalysis(cfg, data);
            avg.subject = subjectpath(ii).name; 
            avg.run = filename(kk).name(end-4);   
            avg.balldesign = cleanMEG_interp.trialinfo.balldesign_short;
            grandavg.like.(subjectpath(ii).name)(counter).avg = avg; 
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
        if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\'], 'dir')
            mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\'])
        end
        save ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_like.mat' ], 'avg_like');
        figure
        plot(avg_like.time, avg_like.avg)
        title('grandavg like')
        savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_like.fig' ]);
        
        
        %% dislike:
         filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\02_interpolated\Neu_Dislike*.mat']);
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
            cfg.baselinewindow  = [-0.5 0];
            data = ft_preprocessing(cfg, cleanMEG_interp);
            cfg = [];
            cfg.latency = [-0.5 1];
            avg = ft_timelockanalysis(cfg, data);
            avg.subject = subjectpath(ii).name; 
            avg.run = filename(kk).name(end-4);   
            avg.balldesign = cleanMEG_interp.trialinfo.balldesign_short;
            grandavg.dislike.(subjectpath(ii).name)(counter).avg = avg; 
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
        if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\'], 'dir')
            mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\'])
        end
        save ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_dislike.mat' ], 'avg_dislike');
        figure
        plot(avg_dislike.time, avg_dislike.avg)
        title('grandavg dislike')
        savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_dislike.fig' ]);
       
         
        
          %% dontcare:
         filename = dir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\02_interpolated\Neu_Dontcare*.mat']);
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
                cfg.baselinewindow  = [-0.5 0];
                data = ft_preprocessing(cfg, cleanMEG_interp);
                cfg = [];
                cfg.latency = [-0.5 1];
                avg = ft_timelockanalysis(cfg, data);
                avg.subject = subjectpath(ii).name; 
                avg.run = filename(kk).name(end-4);   
                avg.balldesign = cleanMEG_interp.trialinfo.balldesign_short;
                grandavg.dontcare.(subjectpath(ii).name)(counter).avg = avg; 
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
            if ~exist([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\'], 'dir')
                mkdir([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\'])
            end
            save ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_dontcare.mat' ], 'avg_dontcare');
            figure
            plot(avg_dontcare.time, avg_dontcare.avg)
            title('grandavg dontcare')
            savefig ([subjectpath(ii).folder filesep subjectpath(ii).name filesep  'MEG_analysis\noisereduced\1_95Hz\grandavg\grandavg_dontcare.fig' ]);
         end
        
    end
    
else
    counter_like=1;
    counter_dislike=1;
    cfg = [];
    for i=1:length(subjectpath)
        filename = dir([subjectpath(i).folder filesep subjectpath(i).name filesep  'MEG_analysis\noisereduced\1_95Hz\02_interpolated\Neu_Like*.mat']);

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