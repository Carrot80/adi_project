function adi_crossvalidation_allRuns_perSubj_EEG (data_like, data_dislike, outPath, freqbandname, condition1, condition2, latency)

if ~exist([outPath, 'like_vs_dislike_allRuns_', freqbandname, '.fig'])

    time =  latency(1,:);
    NameCond1vs2 = ([condition1, '_vs_', condition2]);

    if ~exist (outPath, 'dir')
        mkdir(outPath)
    end

    if 0 == strcmp(freqbandname, 'bp1-45Hz') 
            [data_bpfreq_like] = adi_bpfilter_EEG_statistics(data_like, freqbandname)
            [data_bpfreq_dislike] = adi_bpfilter_EEG_statistics(data_dislike, freqbandname)
    else
        data_bpfreq_like = data_like;
        data_bpfreq_dislike = data_dislike;
    end

    [Condition1vs2] = adi_crossvalidation (data_bpfreq_like, data_bpfreq_dislike, freqbandname, NameCond1vs2, latency, outPath)
    adi_figureTPRcrossval (Condition1vs2, time, condition1, condition2, NameCond1vs2, outPath, freqbandname, '_allRuns')

end
end

function [Condition1vs2] = adi_crossvalidation (data_condition1, data_condition2, freqbandname, NameCond1vs2, latency, outPath)

cfg             = [];
cfg.parameter   = 'trial';
cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
cfg.vartrllength = 2;
tCondition1     = ft_timelockanalysis(cfg, data_condition1);    
tCondition2     = ft_timelockanalysis(cfg, data_condition2); 

cfg         = [];
cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
cfg.statistic = {'accuracy', 'binomial', 'contingency'};
cfg.design  = [ones(size(data_condition1.trial,2),1); 2*ones(size(data_condition2.trial,2),1)]';
cfg.resample = 'true';
% cfg.verbose = true;
% cfg.mva = {dml.standardizer dml.enet('family', 'binomial', 'alpha', 0.8)}; % 

Condition1vs2 = [];
Condition1vs2.Accuracy = [];
Condition1vs2.Binominal = [];
Condition1vs2.Latency = latency;
Condition1vs2.stats = '5f-crossvalidation';

    for i=1:length(latency)  
        cfg.latency = [latency(1,i) latency(2,i)]; 
        stat = ft_timelockstatistics(cfg, tCondition1, tCondition2);
        Condition1vs2.Accuracy(1,i) = stat.statistic.accuracy;
        Condition1vs2.Binominal(1,i) = stat.statistic.binomial;
        Condition1vs2.Contingency{1,i} = stat.statistic.contingency;
    end

    Condition1vs2.latency = latency;
    Condition1vs2.design = cfg.design;
    save ([outPath, NameCond1vs2, '_allRuns_', freqbandname, '.mat'], 'Condition1vs2'); 
end


function [data_bpfreq] = adi_bpfilter_EEG_statistics(data, bpname)


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

[data_bpfreq] = ft_preprocessing(cfg, data);

end