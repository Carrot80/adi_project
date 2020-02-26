%% data preparation: see adi_main_new_data_ab_30_08_19


%% Favoriten separieren und abspeichern pro Proband: 
clear
path2subjects = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
path2file = 'MEG_analysis\noisereduced\1_95Hz\02_interpolated\';
subjects_dir = dir(path2subjects);
subjects_dir(1:2,:) = [];
adi_favouriteball_analysis(subjects_dir, path2file)


%%  Favoriten avg across subjects - analysieren: 
clear
path2subjects = 'E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data_ab_11-10-2019\favorite_balldesigns\';
subjects_dir = dir(path2subjects);
subjects_dir(1:2,:) = [];
adi_favouriteball_analysis_across_subj(subjects_dir)

%% MVPA mit Favoritenanalyse

clear
path2subjects = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';

subjects_dir = dir(path2subjects);
subjects_dir(1:2,:) = [];

adi_leave_out_exemplar_singlsubj_favorites(path2subjects, subjects_dir)

%% read number of trials per subject

subjectpath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
dir_subjectpath = dir(subjectpath);
dir_subjectpath(1:3,:) = [];

for ii = 2:length(dir_subjectpath)
    
    path2datafile = [dir_subjectpath(ii).folder filesep dir_subjectpath(ii).name '\MEG_analysis\noisereduced\1_95Hz\mvpa\favorite_balldesigns\trl_no_favorites_vs_no_favorites.mat'];
    load (path2datafile)
    num_trials(ii,:) = cell2mat(perf.number_of_trials);
    
end

%% performance MPVA Favoritenanalyse - mean for subjects

clear
path2subjects = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subjects_dir = dir(path2subjects);
subjects_dir(1:2,:) = [];
path2file = 'MEG_analysis\noisereduced\1_95Hz\mvpa\favorite_balldesigns\perf_favorites_vs_no_favorites.mat';

adi_leave_out_exemplar_mean_subj_favorites(subjects_dir, path2file)

%% regression

