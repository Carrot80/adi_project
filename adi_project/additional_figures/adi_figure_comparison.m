function [] = adi_figure_comparison(EEGpath, MEGpath, MEG_EEG_path, outpath, freq)

EEGstats = load ([EEGpath 'like_vs_dislike_allRuns_' freq '.mat']);
MEGstats = load ([MEGpath 'like_vs_dislike_allRuns_' freq '.mat']);
MEG_EEG_stats = load ([MEG_EEG_path 'like_vs_dislike_allRuns_' freq '.mat']);


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