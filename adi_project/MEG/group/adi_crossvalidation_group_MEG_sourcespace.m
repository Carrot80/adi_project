function adi_crossvalidation_group_MEG (group_data_like, group_data_dislike, outPath, latency, freq)

if ~exist([outPath, 'like_vs_dislike_group_', freq, '.mat'], 'file') % sollte Ergebnis noch nicht vorhanden sein, soll er Ausgangsdateien laden
    
    if 0 == strcmp(freq, '1_45Hz')
        [data_bpfreq_like] = adi_bpfilter_EEG_statistics(group_data_like, freq);
        [data_bpfreq_dislike] = adi_bpfilter_EEG_statistics(group_data_dislike, freq);
    else
        data_bpfreq_like = group_data_like;
        clear group_data_like
        data_bpfreq_dislike = group_data_dislike;
        clear group_data_dislike        
    end

    cfg = [];
    cfg.latency = [-0.5 1];
    data_bpfreq_like_red = ft_selectdata(cfg, data_bpfreq_like)
    data_bpfreq_dislike_red = ft_selectdata(cfg, data_bpfreq_dislike)
    clear data_bpfreq_like data_bpfreq_dislike
    
    cfg              = [];
    cfg.parameter    = 'trial';
    cfg.keeptrials   = 'yes'; % classifiers operate on individual trials
    cfg.vartrllength = 2;
    tlike            = ft_timelockanalysis(cfg, data_bpfreq_like_red); 
    clear data_bpfreq_like_red
    tdislike         = ft_timelockanalysis(cfg, data_bpfreq_dislike_red); 
    clear data_bpfreq_dislike_red
    
    time =  latency(1,:);
    figure
    plot(tlike.time, tlike.avg)
    title('group MEG like')
    xlabel('time');
    ylabel('activation'); 
    axis tight
    savefig([outPath, 'avg_like_', freq,'.fig'])
    fig = ([outPath, 'avg_like_', freq]);
    print('-dpng', fig);
    
    figure
    plot(tdislike.time, tdislike.avg)
    title('group MEG dislike')
    xlabel('time');
    ylabel('activation'); 
    savefig([outPath, 'avg_dislike_', freq,'.fig'])
    fig = ([outPath, 'avg_dislike_', freq]);
    print('-dpng', fig);

    %% tlike vs tdislike

    [Condition1vs2] = kh_crossvalidation (tlike, tdislike, 'like_vs_dislike_group', freq, latency, outPath);
    adi_figureTPRcrossval (Condition1vs2, time, 'tlike', 'tdislike', 'like_vs_dislike_group', outPath, freq, []);
  
 end
end

function [Condition1vs2] = kh_crossvalidation (tCondition1, tCondition2, NameOfCondition, freq, latency, outPath)

    cfg         = [];
    cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
    cfg.design  = [ones(size(tCondition1.trial,1),1); 2*ones(size(tCondition2.trial,1),1)]';
    cfg.resample = 'true';
    cfg.statistic = {'accuracy', 'binomial', 'contingency'};
    cfg.verbose = true;
%     cfg.mva = {dml.standardizer dml.enet('family', 'binomial', 'alpha', 0.3)}; % 
%   cfg.mva =  dml.one_against_rest('mva', {dml.standardizer() dml.svm()}) 
%     cfg.mva =  dml.one_against_one('mva', {dml.standardizer() dml.svm()});

    Condition1vs2 = [];
    Condition1vs2.Accuracy = [];
    Condition1vs2.Binominal = [];
    Condition1vs2.Latency = latency;
    Condition1vs2.stats = '5f-crossvalidation';
    if isfield(cfg, 'mva')
        Condition1vs2.mva = cfg.mva;
    end

    for i=1:length(latency)  
        cfg.latency = [latency(1,i) latency(2,i)]; 
        stat = ft_timelockstatistics(cfg, tCondition1, tCondition2);
        Condition1vs2.Accuracy(1,i) = stat.statistic.accuracy;
        Condition1vs2.Binominal(1,i) = stat.statistic.binomial;
        Condition1vs2.Contingency{1,i} = stat.statistic.contingency;
    end

  Condition1vs2.latency = latency;
  Condition1vs2.design = cfg.design;
  save ([outPath, NameOfCondition, '_', freq, '.mat'], 'Condition1vs2'); 

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

cfg = [];
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


