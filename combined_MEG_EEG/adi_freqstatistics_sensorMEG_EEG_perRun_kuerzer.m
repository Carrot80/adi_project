function adi_freqstatistics_sensorMEG_EEG_perRun_kuerzer (EEGpath, MEGpath, outPath, freqbandname, condition1, condition2, latency)

    NameCond1vs2 = ([condition1, '_vs_', condition2]);
    if ~exist ([outPath, NameCond1vs2, '500_1_', freqbandname, '.mat'], 'file') 
        adi_freqstatistics_sensor(EEGpath, MEGpath, outPath, freqbandname, condition1, condition2, latency, num2str(1))
    end
    if ~exist ([outPath, NameCond1vs2, '500_2_', freqbandname, '.mat'], 'file') 
        adi_freqstatistics_sensor(EEGpath, MEGpath, outPath, freqbandname, condition1, condition2, latency, num2str(2))
    end
    if ~exist ([outPath, NameCond1vs2, '500_3_', freqbandname, '.mat'], 'file') 
        adi_freqstatistics_sensor(EEGpath, MEGpath, outPath, freqbandname, condition1, condition2, latency, num2str(3))
    end

end

function adi_freqstatistics_sensor(EEGpath, MEGpath, outPath, freqbandname, condition1, condition2, latency, run)

    time =  latency(1,:);
    NameCond1vs2 = ([condition1, '_vs_', condition2]);
    
    if exist([EEGpath, 'like500_', run, '.mat']);
        EEG_condition1 = load([EEGpath, 'like500_', run, '.mat']);
    else
        return
    end
    if exist ([MEGpath, 'like500_', run, '.mat']);
        MEG_condition1 = load([MEGpath, 'like500_', run, '.mat']);
    else 
        return
    end
    if exist ([EEGpath, 'dislike500_', run, '.mat']);
        EEG_condition2 = load([EEGpath, 'dislike500_', run, '.mat']);
    else
        return
    end
    if exist ([MEGpath, 'dislike500_', run, '.mat']);
        MEG_condition2 = load([MEGpath, 'dislike500_', run, '.mat']);
    else
        return
    end

    
    % für MEG sampleinfo aufbauen, sofern keine vorhanden ist: 
    % condition1: 
    if ~isfield(MEG_condition1.cleanMEG_interp, 'sampleinfo')
         MEG_condition1.cleanMEG_interp.sampleinfo = MEG_condition1.cleanMEG_interp.sampleinfo_orig;
        if length(MEG_condition1.cleanMEG_interp.trial) ~= length(MEG_condition1.cleanMEG_interp.sampleinfo_orig) 
            MEG_condition1.cleanMEG_interp.sampleinfo(MEG_condition1.cleanMEG_interp.rejectedTrials, :) = [];
        end
    end
    
     % um MEG und EEG zusammenzufügen, müssen gleiche Anzahl an Trials vorhanden
    % sein:

    bad_trials_EEG = setdiff(MEG_condition1.cleanMEG_interp.sampleinfo(:,1), EEG_condition1.cleanEEG_interp.sampleinfo(:,1));
    bad_trials_MEG = setdiff(EEG_condition1.cleanEEG_interp.sampleinfo(:,1), MEG_condition1.cleanMEG_interp.sampleinfo(:,1));

    if ~isempty(bad_trials_EEG)
        for m = 1:length(bad_trials_EEG)
            [ind(m)] = find (MEG_condition1.cleanMEG_interp.sampleinfo(:,1) == bad_trials_EEG(m))
        end
        for m = 1:length(bad_trials_EEG)
            MEG_condition1.cleanMEG_interp.trial{1, ind(m)} = [];
            MEG_condition1.cleanMEG_interp.time{1, ind(m)} = [];
        end
        MEG_condition1.cleanMEG_interp.trial(cellfun('isempty',MEG_condition1.cleanMEG_interp.trial)) = [];
        MEG_condition1.cleanMEG_interp.time(cellfun('isempty',MEG_condition1.cleanMEG_interp.time)) = [];
        MEG_condition1.cleanMEG_interp.sampleinfo(ind',:) = [];
        clear ind bad_trials_EEG
    end
    if ~isempty(bad_trials_MEG)
         for m = 1:length(bad_trials_MEG)
            [ind(m)] = find (EEG_condition1.cleanEEG_interp.sampleinfo(:,1) == bad_trials_MEG(m))
         end
         for m = 1:length(bad_trials_MEG)
            EEG_condition1.cleanEEG_interp.trial{1, ind(m)} = [];
            EEG_condition1.cleanEEG_interp.time{1, ind(m)} = [];
         end
         EEG_condition1.cleanEEG_interp.trial(cellfun('isempty',EEG_condition1.cleanEEG_interp.trial)) = [];
         EEG_condition1.cleanEEG_interp.time(cellfun('isempty',EEG_condition1.cleanEEG_interp.time)) = [];
         EEG_condition1.cleanEEG_interp.sampleinfo(ind',:) = [];
         clear ind bad_trials_MEG
    end
    
     % condition 2:
     if ~isfield(MEG_condition2.cleanMEG_interp, 'sampleinfo')
         MEG_condition2.cleanMEG_interp.sampleinfo = MEG_condition2.cleanMEG_interp.sampleinfo_orig;
        if length(MEG_condition2.cleanMEG_interp.trial) ~= length(MEG_condition2.cleanMEG_interp.sampleinfo_orig)
            MEG_condition2.cleanMEG_interp.sampleinfo(MEG_condition2.cleanMEG_interp.rejectedTrials, :) = [];
        end
     end
    
    bad_trials_EEG = setdiff(MEG_condition2.cleanMEG_interp.sampleinfo(:,1), EEG_condition2.cleanEEG_interp.sampleinfo(:,1));
    bad_trials_MEG = setdiff(EEG_condition2.cleanEEG_interp.sampleinfo(:,1), MEG_condition2.cleanMEG_interp.sampleinfo(:,1));

    if ~isempty(bad_trials_EEG)
        for m = 1:length(bad_trials_EEG)
            [ind(m)] = find (MEG_condition2.cleanMEG_interp.sampleinfo(:,1) == bad_trials_EEG(m))
        end
        for m = 1:length(bad_trials_EEG)
            MEG_condition2.cleanMEG_interp.trial{1, ind(m)} = [];
            MEG_condition2.cleanMEG_interp.time{1, ind(m)} = [];
        end
        MEG_condition2.cleanMEG_interp.trial(cellfun('isempty',MEG_condition2.cleanMEG_interp.trial)) = [];
        MEG_condition2.cleanMEG_interp.time(cellfun('isempty',MEG_condition2.cleanMEG_interp.time)) = [];
        MEG_condition2.cleanMEG_interp.sampleinfo(ind',:) = [];
        clear ind bad_trials_EEG
    end   
    if ~isempty(bad_trials_MEG)
         for m = 1:length(bad_trials_MEG)
            [ind(m)] = find (EEG_condition2.cleanEEG_interp.sampleinfo(:,1) == bad_trials_MEG(m))
         end
         for m = 1:length(bad_trials_MEG)
            EEG_condition2.cleanEEG_interp.trial{1, ind(m)} = [];
            EEG_condition2.cleanEEG_interp.time{1, ind(m)} = [];
         end
         EEG_condition2.cleanEEG_interp.trial(cellfun('isempty',EEG_condition2.cleanEEG_interp.trial)) = [];
         EEG_condition2.cleanEEG_interp.time(cellfun('isempty',EEG_condition2.cleanEEG_interp.time)) = [];
         EEG_condition2.cleanEEG_interp.sampleinfo(ind',:) = [];
         clear ind bad_trials_MEG
    end
    
     % Referenzsensoren rauswerfen:
    cfg = [];
    cfg.channel = 'meg';
    MEG_condition1 = ft_selectdata(cfg, MEG_condition1.cleanMEG_interp)
    MEG_condition1 = rmfield(MEG_condition1, 'grad');
    MEG_condition2 = ft_selectdata(cfg, MEG_condition2.cleanMEG_interp)
    MEG_condition2 = rmfield(MEG_condition2, 'grad');

%     for m = 1:length(MEG_condition1.cleanMEG_interp.trial)
%         MEG_condition1.cleanMEG_interp.trial{1,m}([249:end],:) = [];    
%     end
% 
%     for m = 1:length(MEG_condition2.cleanMEG_interp.trial)
%         MEG_condition2.cleanMEG_interp.trial{1,m}([249:end],:) = [];
%     end
% 
%     MEG_condition1.cleanMEG_interp.label(249:end) = [];
%     MEG_condition2.cleanMEG_interp.label(249:end) = [];    
    
    if 0 == strcmp(freqbandname, 'bp1-45Hz') 
        [MEG_bpfreq_condition1] = adi_bpfilter_statistics(MEG_condition1, freqbandname, condition1, run)
        [MEG_bpfreq_condition2] = adi_bpfilter_statistics(MEG_condition2, freqbandname, condition2, run)
        [EEG_bpfreq_condition1] = adi_bpfilter_statistics(EEG_condition1.cleanEEG_interp, freqbandname, condition1, run)
        [EEG_bpfreq_condition2] = adi_bpfilter_statistics(EEG_condition2.cleanEEG_interp, freqbandname, condition2, run)
        clear MEG_condition1 MEG_condition2
    else
        % MEG muss von 1:45Hz gefilert werden wie EEG
        cfg = [];
        cfg.trials        = 'all'; 
        cfg.feedback      = 'yes';
        cfg.lpfilter      = 'yes';
        cfg.lpfreq        = 45;
        [MEG_bpfreq_condition1] = ft_preprocessing(cfg, MEG_condition1);   
        [MEG_bpfreq_condition2] = ft_preprocessing(cfg, MEG_condition2);  
        clear MEG_condition1 MEG_condition2   
    end

    % selecttimepoints:
    
    if length(MEG_bpfreq_condition1.trial{1,1}) == 3053
        for k = 1:length(MEG_bpfreq_condition1.trial)
            MEG_bpfreq_condition1.trial{1,k}(:,1) = [];
            MEG_bpfreq_condition1.time{1,k}(:,1) = [];
        end
       
    end
    if length(MEG_bpfreq_condition2.trial{1,1}) == 3053
        for k = 1:length(MEG_bpfreq_condition2.trial)
            MEG_bpfreq_condition2.trial{1,k}(:,1) = [];
            MEG_bpfreq_condition2.time{1,k}(:,1) = [];
        end
    end        
            
    cfg = [];
    cfg.latency = [-0.5 1.2];
    MEG_bpfreq_condition1 = ft_selectdata(cfg, MEG_bpfreq_condition1);
    MEG_bpfreq_condition2 = ft_selectdata(cfg, MEG_bpfreq_condition2);

    EEG_condition1 = ft_selectdata(cfg, EEG_condition1.cleanEEG_interp);
    EEG_condition2 = ft_selectdata(cfg, EEG_condition2.cleanEEG_interp);
     
    
    
    cfg = [];
    [condition1_combined_MEG_EEG] = ft_appenddata (cfg, MEG_bpfreq_condition1, EEG_condition1)
    clear EEG_condition1 MEG_bpfreq_condition1

    [condition2_combined_MEG_EEG] = ft_appenddata (cfg, MEG_bpfreq_condition2, EEG_condition2)
    clear  EEG_condition2  MEG_bpfreq_condition2
    
    [Condition1vs2] = adi_crossvalidation_EEG (condition1_combined_MEG_EEG, condition2_combined_MEG_EEG, freqbandname, NameCond1vs2, latency, ['500_', run], outPath)
    adi_figureTPRcrossval (Condition1vs2, time, condition1, condition2, NameCond1vs2, outPath, freqbandname, ['500_', run])

end



function [Condition1vs2] = adi_crossvalidation_EEG (condition1run, condition2run, freqbandname, NameCond1vs2, latency, run, outPath)

cfg             = [];
cfg.parameter   = 'trial';
cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
cfg.vartrllength = 2;
tCondition1     = ft_timelockanalysis(cfg, condition1run);    
tCondition2     = ft_timelockanalysis(cfg, condition2run); 

cfg         = [];
cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
cfg.statistic = {'accuracy', 'binomial', 'contingency'};
cfg.design  = [ones(size(condition1run.trial,2),1); 2*ones(size(condition2run.trial,2),1)]';
cfg.resample = 'true';
%     cfg.mva = {dml.standardizer dml.enet 'family', 'binomial', 'alpha', 0.3};

Condition1vs2 = [];
Condition1vs2.Accuracy=[];
Condition1vs2.Binominal=[];
Condition1vs2.Latency=latency;
Condition1vs2.stats='5f-crossvalidation';

    for i=1:length(latency)  
        cfg.latency = [latency(1,i) latency(2,i)]; 
        stat = ft_timelockstatistics(cfg, tCondition1, tCondition2);
        Condition1vs2.Accuracy(1,i) = stat.statistic.accuracy;
        Condition1vs2.Binominal(1,i) = stat.statistic.binomial;
        Condition1vs2.Contingency{1,i} = stat.statistic.contingency;
    end

    Condition1vs2.latency = latency;
    Condition1vs2.design = cfg.design;
    save ([outPath, NameCond1vs2, run, '_', freqbandname, '.mat'], 'Condition1vs2'); 
end

function [data_bpfreq] = adi_bpfilter_statistics(data, bpname, condition, run)


switch bpname
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
end

cfg=[];
cfg.trials        = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') 
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

[data_bpfreq] = ft_preprocessing(cfg, data);
end

