%%  
    fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(fieldtripPath);
    ListSubj(1:2) = [];
    filter     = '1_45Hz';


%% comparison of MEG and EEG
% erstellt Abbildung der einzelnen Probanden für MEG und EEG separat

tlike_vs_dislike_all = []
outPathMEG = (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\additional_figures\']);
for i = 1:length(ListSubj)
        subjpath_MEG_stats = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\1_95Hz\04_timelockstatistics_interp\']); 
        [tlike_vs_dislike_all] = figure_comparison_singleSubjects (tlike_vs_dislike_all, subjpath_MEG_stats, ListSubj(i).name)
end

figure_comparison (tlike_vs_dislike_all, outPath, 'delta')
figure_comparison (tlike_vs_dislike_all, outPath, 'theta')
figure_comparison (tlike_vs_dislike_all, outPath, 'alpha')
figure_comparison (tlike_vs_dislike_all, outPath, 'bp1_45Hz')
 

%% comparison of MEG and EEG


tlike_vs_dislike_all = [];
outPathEEG = (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\additional_figures_EEG\']);
for i = 1:length(ListSubj)
        subjpath_EEG_stats = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\1_45Hz\04_timelockstatistics_interp\']); 
        [tlike_vs_dislike_all] = figure_comparison_singleSubjects (tlike_vs_dislike_all, subjpath_EEG_stats, ListSubj(i).name)
end

figure_comparison (tlike_vs_dislike_all, outPathEEG, 'delta', 'EEG')
figure_comparison (tlike_vs_dislike_all, outPathEEG, 'theta', 'EEG')
figure_comparison (tlike_vs_dislike_all, outPathEEG, 'alpha', 'EEG')
figure_comparison (tlike_vs_dislike_all, outPathEEG, 'bp1_45Hz', 'EEG')

%% Abbildung für MEG, EEG MEG_EEG Vergleich pro Proband - alle runs

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_45Hz';

for i = 1 :length(ListSubj)
    EEGstatsPath = ([fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\04_timelockstatistics_interp\']);
    MEGstatsPath = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\04_timelockstatistics_interp\']);
    MEG_EEGstatsPath = ([fieldtripPath ListSubj(i).name '\combined_MEG_EEG_analysis\02c_statistics_allRuns\']);
    outpath = ([fieldtripPath ListSubj(i).name '\comparison_of_all\']);
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end
    adi_figure_comparison(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'delta');
    adi_figure_comparison(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'theta');
    adi_figure_comparison(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'bp1-45Hz');
    adi_figure_comparison(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'alpha');
    adi_figure_comparison(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'beta');
end


%% Welche bzw. wie viele Probanden zeigen Accuracy-Werte über 70%: 

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
EEG_stats = [];
MEG_stats = [];
MEEG_stats = [];
table80 = [];
outpath = (['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\comparison_MEG_EEG_MEEG\best_subjects\']);
[table80] = adi_comparison_subjects(table80, fieldtripPath, outpath, 'delta');
[table80] = adi_comparison_subjects(table80, fieldtripPath, outpath, 'theta');
[table80] = adi_comparison_subjects(table80, fieldtripPath, outpath, 'bp1-45Hz');
[table80] = adi_comparison_subjects(table80, fieldtripPath, outpath, 'alpha');
[table80] = adi_comparison_subjects(table80, fieldtripPath, outpath, 'beta');


    
% figure
% time = EEGstats.Condition1vs2.latency(1,:);
% plot(time, EEGstats.Condition1vs2.Accuracy, 'r')
% hold on
% plot(time, MEGstats.Condition1vs2.Accuracy, 'b')
% hold on
% plot(time, MEG_EEG_stats.Condition1vs2.Accuracy, 'k')
% title(['tlike_vs_tdislike_' freq]);
% ylim([0 1])
% xlabel('time');
% ylabel ('accuracy'); 
% legend('EEG', 'MEG', 'combined') 
% savefig([outpath filesep 'comparison_of_EEG_MEG_MEEG_allRuns_' freq '.fig']);
% fig = ([outpath filesep, 'comparison_of_EEG_MEG_MEEG_allRuns_' freq]);
% print('-dpng', fig); 
% close all    

%% statistik: Unterscheidet sich die Accuracy zwischen  MEG, EEG und combined MEG/EEG pro Proband - was ist besser: MEG, EEG oder kombiniert 
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_45Hz';

for i = 1 :length(ListSubj)
    EEGstatsPath = ([fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\04_timelockstatistics_interp\']);
    MEGstatsPath = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\04_timelockstatistics_interp\']);
    MEG_EEGstatsPath = ([fieldtripPath ListSubj(i).name '\combined_MEG_EEG_analysis\02c_statistics_allRuns\']);
    outpath = ([fieldtripPath ListSubj(i).name '\comparison_of_all\stats']);
    if ~exist(outpath, 'dir')
        mkdir(outpath)
    end
    adi_stats_MEEG(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'delta');
    adi_stats_MEEG(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'theta');
    adi_stats_MEEG(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'bp1-45Hz');
    adi_stats_MEEG(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'alpha');
    adi_stats_MEEG(EEGstatsPath, MEGstatsPath, MEG_EEGstatsPath, outpath, 'beta');
end


