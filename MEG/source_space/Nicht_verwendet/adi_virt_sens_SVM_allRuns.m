function [] = adi_virt_sens_SVM(outPath, path_extdisc, data_all_conditions, conditions, freqbandname, cfg_virtsens)


adi_crossvalidation (data_all_conditions, conditions, freqbandname, outPath, cfg_virtsens);


% 
% NameCond1vs2 = [cfg_virtsens '_' name_condition1 '_vs_' name_condition2 '_allRuns'];
% if ~exist([SubjOutPath NameCond1vs2 '_' freqbandname '.mat'], 'file')  
%     
%     
%     
%     if 1 == isempty(data_condition1)
%         load ([path_extdisc '\MEG\sourcespace\runs_appended\virtsens\' cfg_virtsens '_' name_condition1 '_allRuns_' freqbandname '.mat'], 'vs_allRuns');
%         data_condition1 = vs_allRuns;
%         clear vs_allRuns
%     end
%     
%     if 1 == isempty(data_condition2)
%           load ([path_extdisc '\MEG\sourcespace\runs_appended\virtsens\' cfg_virtsens '_' name_condition2 '_allRuns_' freqbandname '.mat'], 'vs_allRuns');
%           data_condition2 = vs_allRuns;
%           clear vs_allRuns
%     end
%     
%     if isempty (data_condition1) || isempty(data_condition2)
%         return
%     end
%     
%     cond1 = [cfg_virtsens '_' name_condition1 '_allRuns'];
%     cond2 = [cfg_virtsens '_' name_condition2 '_allRuns'];
% 
%     label =  cell(1, length(data_condition1.trial{1,2}(:,1)));
%     for k = 1:length(data_condition1.trial{1,2}(:,1))
%        label{k} = num2str(k);
%     end
%     
%     data_condition1.label = label;
%     data_condition2.label = label;
%     [Condition1vs2] = adi_crossvalidation (data_all_conditions, freqbandname, SubjOutPath);
%     [Condition1vs2] = adi_crossvalidation (data_all_conditions, freqbandname, NameCond1vs2, SubjOutPath);
%     adi_figureTPRcrossval_SVM (Condition1vs2, time, cond1, cond2, NameCond1vs2, SubjOutPath, freqbandname)
% 
%     clear data_condition1 data_condition2
% end

end

function [Condition1vs2] = adi_crossvalidation (data_all_conditions, conditions, freqbandname, outPath, cfg_virtsens)

outPath_SVM_results_dmlt = [outPath 'SVM_results_dmlt\'];
if ~exist(outPath_SVM_results_dmlt, 'dir')
    mkdir(outPath_SVM_results_dmlt)
end

if exist([outPath_SVM_results_dmlt conditions{1} '_vs_' conditions{2} '_' freqbandname '.mat'], 'file')
    return
end

if isempty (data_all_conditions)
    load  ([outPath '\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqbandname, '.mat'], 'vs_allRuns');
    data_all_conditions = vs_allRuns;
end

conditions = fields(vs_allRuns);
cfg             = [];
cfg.parameter   = 'trial';
cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
switch size(fields(data_all_conditions), 1)
    case 3
        tCondition.(conditions{1}) = ft_timelockanalysis(cfg, data_all_conditions.(conditions{1}));
        tCondition.(conditions{2}) = ft_timelockanalysis(cfg, data_all_conditions.(conditions{2}));   
        tCondition.(conditions{3}) = ft_timelockanalysis(cfg, data_all_conditions.(conditions{3}));  
    case 2
        tCondition.(conditions{1}) = ft_timelockanalysis(cfg, data_all_conditions.(conditions{1}));
        tCondition.(conditions{2}) = ft_timelockanalysis(cfg, data_all_conditions.(conditions{2})); 
end
    
cfg         = [];
cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
cfg.statistic = {'accuracy', 'binomial', 'contingency'};
cfg.resample = 'true';
%     cfg.mva = {dml.standardizer dml.enet 'family', 'binomial', 'alpha', 0.3};
cfg.mva = {dml.one_against_one('mva',dml.svm)};
Condition1vs2 = [];
Condition1vs2.Accuracy = [];
Condition1vs2.Binominal = [];
%     Condition1vs2.Latency=latency;
Condition1vs2.stats = '5f-crossvalidation';
latency = [tCondition.(conditions{1}).time; tCondition.(conditions{1}).time];

    
switch size(fields(data_all_conditions), 1)
    case 3 
         cfg.design  = [ones(size(tCondition.(conditions{1}).trial,1),1); 2*ones(size(tCondition.(conditions{2}).trial,1),1)]';
         for i = 1:size(tCondition.(conditions{1}).time(1,:),2)
             cfg.latency = [latency(1,i) latency(2,i)]; 
             stat = ft_timelockstatistics(cfg, tCondition.(conditions{1}), tCondition.(conditions{2})); 
             Condition1vs2.Accuracy(1,i) = stat.statistic.accuracy;
             Condition1vs2.Binominal(1,i) = stat.statistic.binomial;
             Condition1vs2.Contingency{1,i} = stat.statistic.contingency;
         end
         Condition1vs2.latency = latency;
         Condition1vs2.design = cfg.design;
         Condition1vs2.condition = {conditions{1},conditions{2}} ;
         save ([outPath_SVM_results_dmlt conditions{1} '_vs_' conditions{2} '_' freqbandname '.mat'], 'Condition1vs2'); 
         adi_figureTPRcrossval_SVM (Condition1vs2, latency, outPath_SVM_results_dmlt, freqbandname)
         
         cfg.design  = [ones(size(tCondition.(conditions{1}).trial,1),1); 2*ones(size(tCondition.(conditions{3}).trial,1),1)]';
         for i = 1:size(tCondition.(conditions{1}).time(1,:),2)
             stat = ft_timelockstatistics(cfg, tCondition.(conditions{1}), tCondition.(conditions{3})); 
             Condition1vs3.Accuracy(1,i) = stat.statistic.accuracy;
             Condition1vs3.Binominal(1,i) = stat.statistic.binomial;
             Condition1vs3.Contingency{1,i} = stat.statistic.contingency;
         end
         Condition1vs3.latency = latency;
         Condition1vs3.design = cfg.design;             
         Condition1vs3.condition = {conditions{1},conditions{3}} ;
         save ([outPath_SVM_results_dmlt conditions{1} '_vs_' conditions{3} '_' freqbandname '.mat'], 'Condition1vs3'); 
         adi_figureTPRcrossval_SVM (Condition1vs3, latency, outPath_SVM_results_dmlt, freqbandname)
         
         cfg.design  = [ones(size(tCondition.(conditions{2}).trial,1),1); 2*ones(size(tCondition.(conditions{3}).trial,1),1)]';
         for i = 1:size(tCondition.(conditions{1}).time(1,:),2)
             stat = ft_timelockstatistics(cfg, tCondition.(conditions{2}), tCondition.(conditions{3})); 
             Condition2vs3.Accuracy(1,i) = stat.statistic.accuracy;
             Condition2vs3.Binominal(1,i) = stat.statistic.binomial;
             Condition2vs3.Contingency{1,i} = stat.statistic.contingency;
         end
         Condition2vs3.latency = latency;
         Condition2vs3.design = cfg.design; 
         Condition2vs3.condition = {conditions{2},conditions{3}} ;
         save ([outPath_SVM_results_dmlt conditions{2} '_vs_' conditions{3} '_' freqbandname '.mat'], 'Condition2vs3'); 
         adi_figureTPRcrossval_SVM (Condition2vs3, latency, outPath_SVM_results_dmlt, freqbandname)
     case 2 
         cfg.design = [ones(size(tCondition.(conditions{1}).trial,1),1); 2*ones(size(tCondition.(conditions{2}).trial,1),1)]';
         for i = 1:size(tCondition.(conditions{1}).time(1,:),2)
             cfg.latency = [latency(1,i) latency(2,i)]; 
             stat = ft_timelockstatistics(cfg, tCondition.(conditions{1}), tCondition.(conditions{2})); 
             Condition1vs2.Accuracy(1,i) = stat.statistic.accuracy;
             Condition1vs2.Binominal(1,i) = stat.statistic.binomial;
             Condition1vs2.Contingency{1,i} = stat.statistic.contingency;
         end
         Condition1vs2.latency = latency;
         Condition1vs2.design = cfg.design;
         Condition1vs2.condition = {conditions{1},conditions{2}} ;
         save ([outPath_SVM_results_dmlt conditions{1} '_vs_' conditions{2} '_' freqbandname '.mat'], 'Condition1vs2'); 
         adi_figureTPRcrossval_SVM (Condition1vs2, latency, outPath_SVM_results_dmlt, freqbandname)
end


end


function adi_figureTPRcrossval_SVM (Conditions, latency, outPath, freqbandname)

figure
plot(latency(1,:), Conditions.Accuracy);
hold on

indSig = []; 
indSig(1:length(Conditions.Binominal)) = NaN;
sig = find(Conditions.Binominal <= 0.1);
indSig(sig) = Conditions.Binominal(sig);
plot(latency(1,:), indSig,'r+'); 
plot([latency(1) latency(end)],[0.5 0.5])
axis tight
ylim([0 1])
title([Conditions.condition{1} '_vs_'  Conditions.condition{2} '_' freqbandname]);
xlabel('time');
ylabel('accuracy/p-value'); 
savefig([outPath, '\' Conditions.condition{1} '_vs_'  Conditions.condition{2}, '_', freqbandname,'.fig'])
fig = ([outPath, '\'  Conditions.condition{1} '_vs_'  Conditions.condition{2}, '_', freqbandname]);
print('-dpng', fig); 

total_cond1 = Conditions.Contingency{1,1}(1,1) + Conditions.Contingency{1,1}(1,2);
total_cond2 = Conditions.Contingency{1,1}(2,1) + Conditions.Contingency{1,1}(2,2);

for i = 1:length(Conditions.Contingency)
    tpr_cond1(i)= Conditions.Contingency{1,i}(1,1)/total_cond1; % true positive rate
    tpr_cond2(i)= Conditions.Contingency{1,i}(2,2)/total_cond2; % true positive rate
end

figure
plot(latency(1,:), Conditions.Accuracy, 'r');
hold on
plot(latency(1,:), indSig,'r+');  
hold on
plot(latency(1,:), tpr_cond1,'b'); 
hold on
plot(latency(1,:), tpr_cond2,'k'); 
leg_cond1 = (['TPR_', Conditions.condition{1}]);
leg_cond2 = (['TPR_', Conditions.condition{2}]);
legend ({'total TPR', 'p-value of total TPR', leg_cond1, leg_cond2});
title([Conditions.condition{1} '_vs_'  Conditions.condition{2} , '_', freqbandname]);
xlabel('time');
ylabel('accuracy/p-value'); 
ylim ([0 1]);
savefig([outPath Conditions.condition{1} '_vs_'  Conditions.condition{2}  '_' freqbandname '_TPR.fig']);
fig = ([outPath Conditions.condition{1} '_vs_'  Conditions.condition{2}  '_' freqbandname '_TPR']);
print('-dpng', fig); 
close all
    
end



