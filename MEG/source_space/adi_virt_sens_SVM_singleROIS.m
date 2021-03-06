

function [Condition1vs2] = adi_virt_sens_SVM (outPath, path_extdisc, freqbandname, cfg_virtsens)

load ('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\atlas_source_indices.mat')
outPath_SVM_results_dmlt = [outPath 'SVM_results_dmlt\'];
if ~exist(outPath_SVM_results_dmlt, 'dir')
    mkdir(outPath_SVM_results_dmlt)
end

% if exist([outPath_SVM_results_dmlt conditions{1} '_vs_' conditions{2} '_' freqbandname '.mat'], 'file')
%     return
% end

load  ([path_extdisc 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqbandname, '.mat'], 'vs_allRuns');
data_all_conditions = vs_allRuns;
clear vs_allRuns
conditions = fields(data_all_conditions);

%% ROIs extrahieren:
for m=1:length(atlas_downsampled.sources_roi_numbers)
    temp(m,1)=atlas_downsampled.sources_roi_numbers{m,1};
    temp(m,2) = m;
end
ROIs_sorted = sortrows(temp);

for k = 1:size(conditions,1)
    for m=1:length(atlas_downsampled.tissuelabel)
        index_tissuelabel = find (ROIs_sorted(:,1)==m)
        index_sources = ROIs_sorted(index_tissuelabel,2);
        for j = 1:length(data_all_conditions.(conditions{1}).trial)
            vs_roi.trial{1,j} = data_all_conditions.(conditions{1}).trial{1,j}(index_sources,:);
        end
        vs_roi.time = data_all_conditions.(conditions{1}).time;
        vs_roi.label = cellstr(num2cell(index_sources));
        cfg = []; 
        avg = ft_timelockanalysis(cfg, vs_roi); 
        % hier SVM? evtl. zu wenig daten => auch m�glich, alle
        % frequenzb�nder in SVM einzubeziehen 
        % evtl. avg und figure f�r alle ROIs?
        % hier weitermachen
    end
end








%%
cfg             = [];
cfg.parameter   = 'trial';
cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
cfg.vartrllength = 2;
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
cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts f�hren zum gleichen ergebnis
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



