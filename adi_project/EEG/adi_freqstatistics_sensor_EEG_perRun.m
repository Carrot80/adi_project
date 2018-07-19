function adi_freqstatistics_sensorEEG_perRun (inPath, outPath, freqbandname, condition1, condition2, latency)

time =  latency(1,:);
NameCond1vs2 = ([condition1, '_vs_', condition2]);

%run 1:
condition1run1 = ([inPath, condition1, '500_1', '.mat']);
condition2run1 = ([inPath, condition2, '500_1', '.mat']);

if ~exist ([outPath, NameCond1vs2, '500_1_', freqbandname, '.mat'], 'file') && exist(condition1run1, 'file') && exist(condition2run1, 'file') 
    if 0 == strcmp(freqbandname, 'bp1-45Hz') 
        [data_bpfreq_condition1] = adi_bpfilter_EEG_statistics(inPath, freqbandname, condition1, num2str(1))
        [data_bpfreq_condition2] = adi_bpfilter_EEG_statistics(inPath, freqbandname, condition2, num2str(1))
    else 
        load(condition1run1)
        data_bpfreq_condition1 = cleanEEG_interp;
        clear cleanEEG_interp
        load(condition2run1)
        data_bpfreq_condition2 = cleanEEG_interp;
        clear cleanEEG_interp
    end
    [Condition1vs2] = adi_crossvalidation_EEG (data_bpfreq_condition1, data_bpfreq_condition2, freqbandname, NameCond1vs2, latency, '500_1', outPath)
    adi_figureTPRcrossval (Condition1vs2, time, condition1, condition2, NameCond1vs2, outPath, freqbandname, '500_1')
    clear Condition1vs2 data_bpfreq_condition1  data_bpfreq_condition2
end
  

%% run2:
condition1run2 = ([inPath, condition1, '500_2', '.mat']);
condition2run2 = ([inPath, condition2, '500_2', '.mat']);
if ~exist ([outPath, NameCond1vs2, '500_2_', freqbandname, '.mat'], 'file') && exist(condition1run2, 'file') && exist(condition2run2, 'file')
    if 0 == strcmp(freqbandname, 'bp1-45Hz') 
        [data_bpfreq_condition1] = adi_bpfilter_EEG_statistics(inPath, freqbandname, condition1, num2str(1))
        [data_bpfreq_condition2] = adi_bpfilter_EEG_statistics(inPath, freqbandname, condition2, num2str(1))
    else 
        load(condition1run2)
        data_bpfreq_condition1 = cleanEEG_interp;
        clear cleanEEG_interp
        load(condition2run2)
        data_bpfreq_condition2 = cleanEEG_interp;
        clear cleanEEG_interp
    end
    [Condition1vs2] = adi_crossvalidation_EEG (data_bpfreq_condition1, data_bpfreq_condition2, freqbandname, NameCond1vs2, latency, '500_2', outPath)
    adi_figureTPRcrossval (Condition1vs2, time, condition1, condition2, NameCond1vs2, outPath, freqbandname, '500_2')
    clear Condition1vs2 data_bpfreq_condition1  data_bpfreq_condition2
end

%run3:
condition1run3 = ([inPath, condition1, '500_3', '.mat']);
condition2run3 = ([inPath, condition2, '500_3', '.mat']);
if ~exist ([outPath, NameCond1vs2, '500_3_', freqbandname, '.mat'], 'file') && exist(condition1run3, 'file') && exist(condition2run3, 'file')
    if 0 == strcmp(freqbandname, 'bp1-45Hz') 
        [data_bpfreq_condition1] = adi_bpfilter_EEG_statistics(inPath, freqbandname, condition1, num2str(1))
        [data_bpfreq_condition2] = adi_bpfilter_EEG_statistics(inPath, freqbandname, condition2, num2str(1))
    else 
        load(condition1run3)
        data_bpfreq_condition1 = cleanEEG_interp;
        clear cleanEEG_interp
        load(condition2run3)
        data_bpfreq_condition2 = cleanEEG_interp;
        clear cleanEEG_interp
    end
    [Condition1vs2] = adi_crossvalidation_EEG (data_bpfreq_condition1, data_bpfreq_condition2, freqbandname, NameCond1vs2, latency, '500_3', outPath)
    adi_figureTPRcrossval (Condition1vs2, time, condition1, condition2, NameCond1vs2, outPath, freqbandname, '500_3')
end

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
cfg.channel = 'EEG';
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

function [data_bpfreq] = adi_bpfilter_EEG_statistics(inPath, bpname, condition, run)


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
if 1 == strcmp(bpname, 'delta') || 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

load([inPath, condition, '500_', run, '.mat']); 
[data_bpfreq] = ft_preprocessing(cfg, cleanEEG_interp);

end