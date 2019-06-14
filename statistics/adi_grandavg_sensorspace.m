function  [grandavg] = grandavg_sensorspace(subjectpath, grandavg, delete_run)
counter_like=1;
counter_dislike=1;
cfg=[];
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