function adi_virt_sens_SVM(path2data, SubjOutPath, path_extdisc, data_condition1, data_condition2, name_condition1, name_condition2, latency, freqbandname, cfg_virtsens)

NameCond1vs2 = [cfg_virtsens '_' name_condition1 '_vs_' name_condition2 '_allRuns'];
if ~exist([SubjOutPath NameCond1vs2 '_' freqbandname '.mat'], 'file')  
    
    if 1 == isempty(data_condition1)
        load ([path_extdisc '\runs_appended\virtsens\' cfg_virtsens '_' name_condition1 '_allRuns_' freqbandname '.mat'], 'vs_allRuns');
        data_condition1 = vs_allRuns;
        clear vs_allRuns
    end
    
    if 1 == isempty(data_condition2)
          load ([path_extdisc '\runs_appended\virtsens\' cfg_virtsens '_' name_condition2 '_allRuns_' freqbandname '.mat'], 'vs_allRuns');
          data_condition2 = vs_allRuns;
          clear vs_allRuns
    end
    
    cond1 = [cfg_virtsens '_' name_condition1 '_allRuns'];
    cond2 = [cfg_virtsens '_' name_condition2 '_allRuns'];

    label =  cell(1, length(data_condition1.trial{1,2}(:,1)));
    for k = 1:length(data_condition1.trial{1,2}(:,1))
       label{k} = num2str(k);
    end
    
    data_condition1.label = label;
    data_condition2.label = label;
    time = latency(1,:);
    [Condition1vs2] = adi_crossvalidation (data_condition1, data_condition2, freqbandname, NameCond1vs2, latency, SubjOutPath);
    adi_figureTPRcrossval_SVM (Condition1vs2, time, cond1, cond2, NameCond1vs2, SubjOutPath, freqbandname)

    clear data_condition1 data_condition2
end

end

function [Condition1vs2] = adi_crossvalidation (data_condition1, data_condition2, freqbandname, NameCond1vs2, latency, SubjOutPath)

    cfg             = [];
    cfg.parameter   = 'trial';
    cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
    cfg.vartrllength = 2;
    
    tCondition1 = ft_timelockanalysis(cfg, data_condition1);
    tCondition2 = ft_timelockanalysis(cfg, data_condition2); 
    
    cfg         = [];
    cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
    cfg.statistic = {'accuracy', 'binomial', 'contingency'};
    cfg.design  = [ones(size(data_condition1.trial,2),1); 2*ones(size(data_condition2.trial,2),1)]';
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
    
    if ~exist(SubjOutPath, 'dir')
        mkdir(SubjOutPath)
    end
        
    
    save ([SubjOutPath NameCond1vs2 '_' freqbandname '.mat'], 'Condition1vs2'); 


end


function adi_figureTPRcrossval_SVM (Condition1vs2, time, cond1, cond2, NameCond1vs2, outPath, freqbandname)
 
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

savefig([outPath, '\' NameCond1vs2, '_', freqbandname,'.fig'])
fig = ([outPath, '\' NameCond1vs2, '_', freqbandname]);
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
savefig([outPath NameCond1vs2 '_' freqbandname '_TPR.fig']);
fig = ([outPath NameCond1vs2 '_' freqbandname '_TPR']);
print('-dpng', fig); 
close all
    
end



