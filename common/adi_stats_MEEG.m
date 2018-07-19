function [] = adi_figure_comparison(EEGpath, MEGpath, MEG_EEG_path, outpath, freq)

EEGstats = load ([EEGpath 'like_vs_dislike_allRuns_' freq '.mat']);
MEGstats = load ([MEGpath 'like_vs_dislike_allRuns_' freq '.mat']);
MEG_EEG_stats = load ([MEG_EEG_path 'like_vs_dislike_allRuns_' freq '.mat']);

%% prepare data for Anova:

data(:,1) = EEGstats.Condition1vs2.Accuracy(51:149);
data(:,2) = MEG_EEG_stats.Condition1vs2.Accuracy(51:149);
data(:,3) = MEGstats.Condition1vs2.Accuracy(51:149);
group = {'EEG', 'MEG', 'MEEG'};

[p,tbl,stats] = anova1(data, group, 'on')
[results,means] = multcompare(stats,'CType','bonferroni')

p = kruskalwallis(data, group)
p = friedman(data)

data2(:,1) = Condition1vs2_EEG.Accuracy(51:149);
data2(:,2) = Condition1vs2_MEG.Accuracy(51:149);

% Pearson linear correlation coefficient 
rho = corr(data2) 
[rho,pval] = corr(data2)
rho = corr(data)
%
[rho,pval] = corr(data2,'Type','Kendall','Rows','complete') 
[rho,pval] = corr(data2,'Type','Spearman') 
% repeated measures ANOVA:
[p, table] = anova_rm(data, displayopt) 

figure
time = EEGstats.Condition1vs2.latency(1,:);
plot(time, EEGstats.Condition1vs2.Accuracy, 'r')
hold on
plot(time, MEGstats.Condition1vs2.Accuracy, 'b')
hold on
plot(time, MEG_EEG_stats.Condition1vs2.Accuracy, 'k')
title(['tlike_vs_tdislike_' freq]);
ylim([0 1])
xlabel('time');
ylabel ('accuracy'); 
legend('EEG', 'MEG', 'combined') 
savefig([outpath filesep 'comparison_of_EEG_MEG_MEEG_allRuns_' freq '.fig']);
fig = ([outpath filesep, 'comparison_of_EEG_MEG_MEEG_allRuns_' freq]);
print('-dpng', fig); 
close all


end