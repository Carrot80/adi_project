function adi_crossvalidation_group_MEG (inPath, outPath, latency, freq)

if ~exist([outPath, 'like_vs_dislike_group_', freq, '.mat'], 'file') % sollte Ergebnis noch nicht vorhanden sein, soll er Ausgangsdateien laden
    load([inPath, 'dislike_group_', freq, '.mat']);
    load([inPath, 'like_group_', freq, '.mat']);
    if exist('dislike_all_subj_red', 'var') % bei einer Ausgangsdatei hatte ich schon ft_appenddata angewandt, bei allen anderen hatte ich die trials der einzelnen Probanden untereinander angefügt (bietet besseren Überblick, man kann einzelne Probanden auseinanderhalten)
        size_dislike_all_subj_red = size(dislike_all_subj_red);
        cfg = [];
        dislike_all_subj_red_appended = ft_appenddata(cfg, dislike_all_subj_red(:,1), dislike_all_subj_red(:,2));

        for k = 3:size_dislike_all_subj_red (2)
            dislike_all_subj_red_appended = ft_appenddata(cfg, dislike_all_subj_red_appended, dislike_all_subj_red(:,k)) ;   
        end
        clear dislike_all_subj_red
    end
    if exist('like_all_subj_red', 'var')
        size_like_all_subj_red=size(like_all_subj_red);
        cfg = [];
        like_all_subj_red_appended = ft_appenddata(cfg, like_all_subj_red(:,1), like_all_subj_red(:,2));

        for k = 3:size_like_all_subj_red (2)
            like_all_subj_red_appended = ft_appenddata(cfg, like_all_subj_red_appended, like_all_subj_red(:,k)) ;   
        end
        clear like_all_subj_red
    end
    cfg              = [];
    cfg.parameter    = 'trial';
    cfg.keeptrials   = 'yes'; % classifiers operate on individual trials
    cfg.vartrllength = 2;
    tlike            = ft_timelockanalysis(cfg, like_all_subj_red_appended); 
    clear like_all_subj_red_appended
    tdislike         = ft_timelockanalysis(cfg, dislike_all_subj_red_appended); 
    clear dislike_all_subj_red_appended
    
    time =  latency(1,:);
    figure
    plot(tlike.time, tlike.avg(1:248,:))
    title('group MEG like')
    xlabel('time');
    ylabel('activation'); 
    axis tight
    savefig([inPath, 'avg_like_', freq,'.fig'])
    fig = ([inPath, 'avg_like_', freq]);
    print('-dpng', fig);
    
    figure
    plot(tdislike.time, tdislike.avg(1:248,:))
    title('group MEG dislike')
    xlabel('time');
    ylabel('activation'); 
    savefig([inPath, 'avg_dislike_', freq,'.fig'])
    fig = ([inPath, 'avg_dislike_', freq]);
    print('-dpng', fig);

    %% tlike vs tdislike

    [Condition1vs2] = kh_crossvalidation (tlike, tdislike, 'like_vs_dislike_group_30ms', freq, latency, outPath);
    adi_figureTPRcrossval (Condition1vs2, time, 'tlike', 'tdislike', 'like_vs_dislike_group_30ms', outPath, freq, []);
  
 end
end

function [Condition1vs2] = kh_crossvalidation (tCondition1, tCondition2, NameOfCondition, freq, latency, outPath)

    cfg         = [];
    cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
    cfg.channel = 'MEG';
    cfg.design  = [ones(size(tCondition1.trial,1),1); 2*ones(size(tCondition2.trial,1),1)]';
    cfg.resample = 'true';
    cfg.statistic = {'accuracy', 'binomial', 'contingency'};
    cfg.compact = 'true';
    cfg.mva = dml.one_against_one('mva', {dml.standardizer() dml.svm()}); 

    Condition1vs2 = [];
    Condition1vs2.Accuracy = [];
    Condition1vs2.Binominal = [];
    Condition1vs2.Latency = latency;
    Condition1vs2.stats = '5f-crossvalidation';

    for i=1:length(latency)  
        cfg.latency = [latency(1,i) latency(2,i)]; 
        stat = ft_timelockstatistics(cfg, tCondition1, tCondition2);
        Condition1vs2.Accuracy(1,i)=stat.statistic.accuracy;
        Condition1vs2.Binominal(1,i)=stat.statistic.binomial;
        Condition1vs2.Contingency{1,i}=stat.statistic.contingency;
    end

  Condition1vs2.latency = latency;
  Condition1vs2.design = cfg.design;
  save ([outPath, NameOfCondition, freq, '.mat'], 'Condition1vs2'); 

end