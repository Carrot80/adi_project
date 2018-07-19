
function adi_source_SVM (path2vol, path2data, mriPath, SubjOutPath, outPath_extdisc, freqbandname, like, dislike, latency, prestim, poststim)

time =  latency(1,:);
% if ~exist([SubjOutPath '\run1\SVM_result\virtSens\virtsens_like_vs_dislike_1_' freqbandname '.mat'])
      
    [Condition1vs2] = adi_crossvalidation (virtsens_like, virtsens_dislike, freqbandname, 'virtsens_like_vs_dislike', latency, num2str(1), SubjOutPath)
    adi_figureTPRcrossval_SVM (Condition1vs2, time, like, dislike, 'virtsens_like_vs_dislike', SubjOutPath, freqbandname, num2str(1))

    [Condition1vs2] = adi_crossvalidation (virtsens_ns_like, virtsens_ns_dislike, freqbandname, 'virtsens_ns_like_vs_dislike', latency, num2str(1), SubjOutPath)
    adi_figureTPRcrossval_SVM (Condition1vs2, time, like, dislike, 'virtsens_ns_like_vs_dislike', SubjOutPath, freqbandname, num2str(1))

    clear virtsens_like virtsens_ns_like virtsens_dislike virtsens_ns_dislike

end



function [Condition1vs2] = adi_crossvalidation (condition1run, condition2run, freqbandname, NameCond1vs2, latency, run, SubjOutPath)

    cfg             = [];
    cfg.parameter   = 'trial';
    cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
    % cfg.channel     = 'MEG';
    cfg.vartrllength = 2;
    tCondition1     = ft_timelockanalysis(cfg, condition1run);    
    tCondition2     = ft_timelockanalysis(cfg, condition2run); 

    cfg         = [];
    cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
    % cfg.channel = 'MEG';
    cfg.statistic = {'accuracy', 'binomial', 'contingency'};
    cfg.design  = [ones(size(condition1run.trial,2),1); 2*ones(size(condition2run.trial,2),1)]';
    cfg.resample = 'true';
    %     cfg.mva = {dml.standardizer dml.enet 'family', 'binomial', 'alpha', 0.3};

    Condition1vs2 = [];
    Condition1vs2.Accuracy=[];
    Condition1vs2.Binominal=[];
    Condition1vs2.Latency=latency;
    Condition1vs2.stats='5f-crossvalidation';

    for i = 1:length(latency)  
        cfg.latency = [latency(1,i) latency(2,i)]; 
        stat = ft_timelockstatistics(cfg, tCondition1, tCondition2);
        Condition1vs2.Accuracy(1,i) = stat.statistic.accuracy;
        Condition1vs2.Binominal(1,i) = stat.statistic.binomial;
        Condition1vs2.Contingency{1,i} = stat.statistic.contingency;
    end

    Condition1vs2.latency = latency;
    Condition1vs2.design = cfg.design;

    outPath = ([SubjOutPath 'run' run '\SVM_result\virtSens' ]);
    if ~exist(outPath, 'file')
        mkdir (outPath)
    end
    save ([outPath filesep NameCond1vs2 '_' run '_' freqbandname '.mat'], 'Condition1vs2'); 


end


function [] = adi_figureTPRcrossval_SVM (Condition1vs2, time, cond1, cond2, NameCond1vs2, outPath, freqbandname, run)
 
figure
plot(time, Condition1vs2.Accuracy);
hold on

indSig = []; 
indSig(1:length(Condition1vs2.Binominal)) = NaN;
sig = find(Condition1vs2.Binominal <= 0.1);
indSig(sig) = Condition1vs2.Binominal(sig);
plot(time, indSig,'r+'); 

title([NameCond1vs2, '_', freqbandname]);
xlabel('time');
ylabel('accuracy/p-value'); 
ylim([0 1]);

savefig([outPath, '\run' run '\SVM_result\virtSens\' NameCond1vs2, run, '_', freqbandname,'.fig'])
fig = ([outPath, '\run' run '\SVM_result\virtSens\' NameCond1vs2, run, '_', freqbandname]);
print('-dpng', fig); 

total_cond1 = Condition1vs2.Contingency{1,1}(1,1) + Condition1vs2.Contingency{1,1}(1,2);
total_cond2 = Condition1vs2.Contingency{1,1}(2,1) + Condition1vs2.Contingency{1,1}(2,2);

for i = 1:length(Condition1vs2.Contingency)
    tpr_cond1(i)= Condition1vs2.Contingency{1,i}(1,1)/total_cond1; % true positive rate
    tpr_cond2(i)= Condition1vs2.Contingency{1,i}(2,2)/total_cond2; % true positive rate
end

figure
plot(time, Condition1vs2.Accuracy, 'r');
hold on
plot(time, indSig,'r+');  
hold on
plot(time, tpr_cond1,'b'); 
hold on
plot(time, tpr_cond2,'k'); 
leg_cond1 = (['TPR_', cond1]);
leg_cond2 = (['TPR_', cond2]);
legend ({'total TPR', 'p-value of total TPR', leg_cond1, leg_cond2});
title([NameCond1vs2, '_', freqbandname]);
xlabel('time');
ylabel('accuracy/p-value'); 
ylim ([0 1]);
savefig([outPath '\run' run '\SVM_result\virtSens\' NameCond1vs2 '_' run '_' freqbandname '_TPR.fig']);
fig = ([outPath '\run' run '\SVM_result\virtSens\' NameCond1vs2 '_' run '_' freqbandname '_TPR']);
print('-dpng', fig); 
close all
    
end
