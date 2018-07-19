function adi_crossvalidation_group_MEG (inPath, outPath, latency, freq)

if ~exist([outPath, 'like_vs_dontcare_paired_group_', freq, '.mat'], 'file')
    load([inPath, 'like_appended_paired_', freq, '.mat']);
    load([inPath, 'dontcare_appended_paired_', freq, '.mat']);
    cfg              = [];
    cfg.parameter    = 'trial';
    cfg.keeptrials   = 'yes'; % classifiers operate on individual trials
    cfg.vartrllength = 2;
    tlike            = ft_timelockanalysis(cfg, like_appended_paired); 
    clear like_appended_paired
    tdontcare         = ft_timelockanalysis(cfg, dontcare_appended_paired); 
    clear dontcare_appended_paired
    
    time = latency(1,:);
    figure
    plot(tlike.time, tlike.avg(1:248,:))
    title('paired group like (MEG)')
    xlabel('time');
    ylabel('activation'); 
    axis tight;
    savefig([inPath, 'avg_like_', freq,'.fig'])
    fig = ([inPath, 'avg_like_', freq]);
    print('-dpng', fig);
    
    figure
    plot(tdontcare.time, tdontcare.avg(1:248,:))
    title('paired group dontcare (MEG)')
    xlabel('time');
    ylabel('activation'); 
    axis tight
    savefig([inPath, 'avg_dontcare_', freq,'.fig'])
    fig = ([inPath, 'avg_dontcare_', freq]);
    print('-dpng', fig);

    % tlike vs tdontcare
    [Condition1vs2] = adi_crossvalidation (tlike, tdontcare, 'like_vs_dontcare_paired', freq, latency, outPath)
    adi_figureTPRcrossval_groupstatistics (Condition1vs2, time, 'like', 'dontcare', 'like_vs_dontcare_paired', outPath, freq);
    clear tlike
end

%% tdislike vs tdontcare
     
if ~exist([outPath, 'dislike_vs_dontcare_paired_group_', freq, '.mat'], 'file')
    load([inPath, 'dislike_appended_paired_', freq, '.mat']);
    tdislike = ft_timelockanalysis(cfg, dislike_appended_paired); 
    clear dislike_appended_paired
    
    figure
    plot(tdislike.time, tdislike.avg(1:248,:))
    title('paired group dislike (MEG)')
    xlabel('time');
    ylabel('activation'); 
    axis tight;
    savefig([inPath, 'avg_dislike_', freq,'.fig'])
    fig = ([inPath, 'avg_dislike_', freq]);
    print('-dpng', fig);
    close all
    [Condition1vs2] = adi_crossvalidation (tdislike, tdontcare, 'dislike_vs_dontcare_paired', freq, latency, outPath);
    adi_figureTPRcrossval_groupstatistics (Condition1vs2, time, 'dislike', 'dontcare', 'dislike_vs_dontcare_paired', outPath, freq);
    
end  

end

function [Condition1vs2] = adi_crossvalidation (tCondition1, tCondition2, NameOfCondition, freq, latency, outPath)

    cfg         = [];
    cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
    cfg.channel = 'MEG';
    cfg.design  = [ones(size(tCondition1.trial,1),1); 2*ones(size(tCondition2.trial,1),1)]';
    cfg.resample = 'true';
    cfg.statistic = {'accuracy', 'binomial', 'contingency'};

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
  save ([outPath, NameOfCondition, '_group_', freq, '.mat'], 'Condition1vs2'); 

end