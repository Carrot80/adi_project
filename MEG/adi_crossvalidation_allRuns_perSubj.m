function adi_crossvalidation_allRuns_perSubj (inPath, outPath, freqbandname, condition1, condition2, latency)

time =  latency(1,:);
NameCond1vs2 = ([condition1, '_vs_', condition2]);

if ~exist (outPath, 'dir')
    mkdir(outPath)
end
condition1run1 = ([inPath, condition1, '_allRuns_', freqbandname, '.mat']);
condition2run1 = ([inPath, condition2, '_allRuns_', freqbandname, '.mat']);

if ~exist ([outPath, NameCond1vs2, '_allRuns_', freqbandname, '.mat'], 'file') && exist(condition1run1, 'file') && exist(condition2run1, 'file')
    [Condition1vs2] = adi_crossvalidation (condition1run1, condition2run1, freqbandname, NameCond1vs2, latency, outPath)
    adi_figureTPRcrossval (Condition1vs2, time, condition1, condition2, NameCond1vs2, outPath, freqbandname, '_allRuns')
end

end

function [Condition1vs2] = adi_crossvalidation (condition1run, condition2run, freqbandname, NameCond1vs2, latency, outPath)

condition1_data = load (condition1run);
condition2_data = load (condition2run);

cfg             = [];
cfg.parameter   = 'trial';
cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
cfg.channel     = 'MEG';
cfg.vartrllength = 2;
tCondition1     = ft_timelockanalysis(cfg, condition1_data.like_allRuns);    
tCondition2     = ft_timelockanalysis(cfg, condition2_data.dislike_allRuns); 

cfg         = [];
cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
cfg.channel = 'MEG';
cfg.statistic = {'accuracy', 'binomial', 'contingency'};
cfg.design  = [ones(size(condition1_data.like_allRuns.trial,2),1); 2*ones(size(condition2_data.dislike_allRuns.trial,2),1)]';
cfg.resample = 'true';
cfg.verbose = true;
cfg.mva = {dml.standardizer dml.enet('family', 'binomial', 'alpha', 0.8)}; % 

Condition1vs2 = [];
Condition1vs2.Accuracy=[];
Condition1vs2.Binominal=[];
Condition1vs2.Latency=latency;
Condition1vs2.stats='5f-crossvalidation';

    for i=1:length(latency)  
        cfg.latency = [latency(1,i) latency(2,i)]; 
        stat = ft_timelockstatistics(cfg, tCondition1, tCondition2);
        Condition1vs2.Accuracy(1,i)=stat.statistic.accuracy;
        Condition1vs2.Binominal(1,i)=stat.statistic.binomial;
        Condition1vs2.Contingency{1,i}=stat.statistic.contingency;
    end

    Condition1vs2.latency = latency;
    Condition1vs2.design = cfg.design;
    save ([outPath, NameCond1vs2, '_allRuns_', freqbandname, '.mat'], 'Condition1vs2'); 
end


