%% 4.12.2018
% neue Daten


%% main settings:
clear
% brainstormPath     = 'E:\adidas\Daten_Franzi\Database12-2018\BS_TrainingGroup\brainstorm_db\';
brainstormPath  = 'W:\neurochirurgie\science\Franzi\Database\BS_TrainingGroup\brainstorm_db\';
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '0.5_95Hz';
    
 %% Export Brainstorm-Files to Fieldtrip:   
    
 for i = 1%:length(ListSubj)    
     path_export_bst2ft   = ([fieldtripPath ListSubj(i).name '\MEG_EEG_input\noisereduced\' filter '\02_Export_Bst2Ft\']);
     if ~exist(path_export_bst2ft, 'dir')
         mkdir(path_export_bst2ft)
     end  
    adi_select_files_newData (brainstormPath, ListSubj(i).name, path_export_bst2ft)
 end
 
 
 
%% füge trigger und balldesign hinzu und entferne bad channels:

clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter  = '0.5_95Hz';

trigger.balldesign_short = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.balldesign = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggercodes = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime        = [102, 104, 106, 108, 101, 103, 105, 107, 109];

for i = 27:length(ListSubj)  
    path2cleanfile = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\01_clean\']);
    path2retval = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\' filter '\02_Export_Bst2Ft\']);
    adi_artifact_cleaning_newData(path2retval, path2cleanfile, ListSubj(i).name, filter, trigger)    
end

%% rejectvisual zusätzliche Säuberung einzelner Runs: für neue Daten nicht genutzt
% nl_adi_19, nl_adi21        
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];

    for i = 18% : length(ListSubj)  
        
        path2cleanfile = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\1_95Hz\01_clean\']);
        path_interpolated = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\']);
        adi_rejectvisual_MEG_extra (path2cleanfile, path_interpolated, ListSubj(i).name)
    end
       
%% für neue Daten nicht genutzt:
% nl_adi_33 => A248 rauswerfen
% nl_adi_19 => dontcare500_1: ab 1.28 sekunden Artefakt => Zeitintervall begrenzen:

clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
file = ([fieldtripPath, ListSubj(16).name, '\MEG_analysis\noisereduced\1_95Hz\01_clean\dontcare500_1.mat']);
load (file)  
cfgn                = [];
cfgn.parameter      = 'trial';
cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
cfgn.vartrllength   = 2;
tcleanMEG         = ft_timelockanalysis(cfgn,cleanMEG);  
figure
plot(tcleanMEG.time, tcleanMEG.avg(1:248,:)) % MEG
for k = 1:length(cleanMEG.trial)
cleanMEG.trial{k}(:,2310:end) = [];
cleanMEG.time{k}(:,2310:end) = [];    
end

tcleanMEG = ft_timelockanalysis(cfgn,cleanMEG);  
figure
plot(tcleanMEG.time, tcleanMEG.avg(1:248,:)) % MEG
    
% nl_adi_21:
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
file = ([fieldtripPath, 'nl_adi_21\MEG_analysis\noisereduced\1_95Hz\01_clean\dislike500_1.mat']);
load (file)  
cfgn                = [];
cfgn.parameter      = 'trial';
cfgn.keeptrials     = 'yes'; % classifiers operate on individual trials
cfgn.vartrllength   = 2;
tcleanMEG         = ft_timelockanalysis(cfgn,cleanMEG);  
figure
plot(tcleanMEG.time, tcleanMEG.avg(1:248,:)) % MEG
cfg = [];
cfg.method = 'summary';%'trial'% 'summary' %, 'trial'
cfg.channel = 'MEG'; % MEG
cfg.keepchannel = 'nan';
cfg.latency = [-0.5 0];
cfg.megscale = 1;
cfg.eegscale = 5e-8;
cfg.interactive = 'yes';
cfg.alim     = 1e-12; 
[cleanMEG]       = ft_rejectvisual(cfg, cleanMEG); 
tcleanMEG         = ft_timelockanalysis(cfgn,cleanMEG);  
figure
plot(tcleanMEG.time, tcleanMEG.avg(1:248,:)) % MEG
[cleanMEG]       = ft_rejectvisual(cfg, cleanMEG); 
 tcleanMEG         = ft_timelockanalysis(cfgn,cleanMEG);  
figure
plot(tcleanMEG.time, tcleanMEG.avg(1:248,:)) % MEG

% 
 
%% interpolate missing channels:
clear
brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\';
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '0.5_95Hz';

for i = 26%:length(ListSubj)  
    path2cleanfile = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\' filter '\01_clean\']);
    pathInterpolated = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\' filter '\02_interpolated\']);

    adi_interpolate_MEG (path2cleanfile, pathInterpolated)
end



%% überprüfe, ob 50ms trigger fälschlicherweise zugeordet wurden:

clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter  = '1_95Hz';


trigger.balldesign_short = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.balldesign = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggercodes = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime        = [102, 104, 106, 108, 101, 103, 105, 107, 109];

for i=1:length(ListSubj)
    path2cleanfile = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\0.5_95Hz\01_clean\'];
    adi_check_for_50ms_trigger(path2cleanfile, trigger, ListSubj(i).name)
end

%% export warped anatomy:
clear
brainstormPath     = 'E:\adidas\Daten_Franzi\Database12-2018\BS_TrainingGroup\brainstorm_db\';
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
MEG_filter     = '1_95Hz';

 for i = 3%:length(ListSubj)    

    outpath_T1warped = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
    if ~exist(outpath_T1warped, 'dir')
        mkdir(outpath_T1warped)
    end
    T1warped_Bst   =  load([brainstormPath ListSubj(i).name '\subjectimage_T1_warped.mat']);
    ftMri = out_fieldtrip_mri(T1warped_Bst, 'anatomy');
    save([outpath_T1warped 'T1warped'], 'ftMri')

 end
 
 %% read and realign mri   
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];

for i = 16:length(ListSubj) 
    mniPath = [fieldtripPath ListSubj(i).name '\MEG_EEG_input\T1_warped\'];
    outPath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\0.5_95Hz\vol\'];
    if ~exist(outPath)
        mkdir (outPath)
    end
    adi_prepare_hdm(mniPath, outPath, ListSubj(i).name)
end



%% compute lcmv Beamformer and virtual sensors for each subject- 
% calculations of spatial filter, multiplies it with avgdata a

clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 26%:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep 'MEG\sourcespace\']);
     condition = {'like', 'dislike', 'dontcare'};
   
     adi_virtsens_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'bp1_45Hz', condition)
   
   
end

%%

%% 29.5.2019: sensordaten für einzelne Bälle zusammenfassen mit neuen Verknüpfungsregeln

clear
balldesign = 'ggv';
path2excelsheet = 'U:\My Documents\MATLAB\eigene_Skripte\adi_project\excelsheets\subjects.xlsx';
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
[delete_run] = adi_boolean_operations(path2excelsheet, subject_list, balldesign);

freq = 'bp1_45Hz';

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];

% [grandavg] = adi_grandavg_sensorspace(subject_list, grandavg, delete_run);

load('U:\My Documents\MATLAB\atlas_source_indices.mat');
sensordata_all_subj = struct('trial', [], 'time', [], 'response_label', [], 'balldesign_short', [], 'run', [], 'subject', [],'label', [], 'fsample', [] ,'grad', [], 'cfg', []);

for i=1:length(subject_list)
    path2sensordata = ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\' subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\' ];
    [sensordata_all_subj] = adi_single_ball_sensordata(sensordata_all_subj, path2sensordata, subject_list(i).name, pattern, trigger, freq, balldesign, delete_run, i);
end

save(['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data\' balldesign '_ball\sensordata_all_subj.mat'], 'sensordata_all_subj')

% nachträgliche Interpolation von Kanälen, da es nicht überall geklappt
% hat:

[sensordata_all_subj_int] = adi_additionalMEG_interpolation(sensordata_all_subj, []);

% z-transformation:
[sensordata_all_subj_zscore] = adi_ztrans_sensorspace(sensordata_all_subj_int);

% clearvars sensordata_all_subj

% statistic
% time=[0.5 0.8];
% adi_stats_sensorspace(sensordata_all_subj_zscore, time) % hier noch vervollstänidgen
% 

session = sensordata_all_subj_zscore(1);
   for k = 2:length(sensordata_all_subj_zscore)
       session.trial = cat(2,session.trial, sensordata_all_subj_zscore(k).trial);
       session.time = cat(2,session.time, sensordata_all_subj_zscore(k).time);
       session.response_label = cat(2,session.response_label, sensordata_all_subj_zscore(k).response_label);
       session.balldesign_short = cat(2,session.balldesign_short, sensordata_all_subj_zscore(k).balldesign_short);
       session.subject = cat(2,session.subject, sensordata_all_subj_zscore(k).subject);
   end
   

for k=1:length(session.response_label)
    switch session.response_label{k}
        case 'like' % 'Volley'
            session.labels(k) = 1;
        case 'Neu_Like' % 'Volley'
            session.labels(k) = 1;
        case 'dislike' % 'Space'
            session.labels(k) = 2;
        case 'Neu_Dislike' % 'Space'
            session.labels(k) = 2;
    end
end

% wie häufig gibt es like und dislike?
n_subj_like = unique(session.subject(find(session.labels==1)));
n_subj_dislike = unique(session.subject(find(session.labels==2)));


% Daten aufteilen, für test data ca. 4 Probanden reservieren  

   % number of trial for learning (training data):
%    size_training_data = round(size(session.trial,2)*0.6)
   % trials 1:203 == training data
   % evtl. training: je 2 probanden mit 21 trials dislike, sowie 2 mit 21 trials like, 
   % Costfunktion plotten (Fehler), 
   % dann cv: jeweils 1 Proband pro Bedingung mit 21 trials  
   % dann langsam vortasten

% time = [0.05 0.16];
path2save = ['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data\' balldesign '\zscore\MVPA\';
time = [];   
initial_data = split_data(session); 


for kk = 1:length(initial_data.order_subj_final.like)
    switch kk
        case 1
         [added_data] = split_data(session, initial_data.like_training, initial_data.dislike_training, initial_data.order_subj_final.like{kk}, initial_data.order_subj_final.dislike{kk});   
%         adi_mvpa_sensorspace_kh_time(added_data.like_training, added_data.dislike_training, initial_data.like_cv, initial_data.dislike_cv, time, path2save)
        case {7, 8}
        [added_data] = split_data(session, added_data.like_training, added_data.dislike_training, initial_data.order_subj_final.like{kk}, []);   
%         adi_mvpa_sensorspace_kh_time(added_data.like_training, added_data.dislike_training, initial_data.like_cv, initial_data.dislike_cv, time, path2save)
        otherwise
        [added_data] = split_data(session, added_data.like_training, added_data.dislike_training, initial_data.order_subj_final.like{kk}, initial_data.order_subj_final.dislike{kk});   
%       
    end
end

%% 29.5.19 compute lda, logreg test and train:

adi_mvpa_sensorspace_lda_train_and_test(added_data.like_training, added_data.dislike_training, initial_data.like_cv, initial_data.dislike_cv, time)
adi_mvpa_sensorspace_logreg_train_and_test(added_data.like_training, added_data.dislike_training, initial_data.like_cv, initial_data.dislike_cv, time)

%% 29.5.19 compute lda, logreg crossvalidate:



cfg=[];
like = ft_appenddata(cfg, added_data.like_training, initial_data.like_cv);
dislike =  ft_appenddata(cfg, added_data.dislike_training, initial_data.dislike_cv);

adi_mvpa_sensorspace_lda_crossvalidate(added_data.like_training, added_data.dislike_training, initial_data.like_cv, initial_data.dislike_cv, time)
adi_mvpa_sensorspace_logreg_crossvalidate(added_data.like_training, added_data.dislike_training, initial_data.like_cv, initial_data.dislike_cv, time)

% original:
adi_mvpa_sensorspace(added_data.like_training, added_data.dislike_training, initial_data.like_cv, initial_data.dislike_cv, time) 
cfg=[];
like = ft_appenddata(cfg, added_data.like_training, initial_data.like_cv);
dislike =  ft_appenddata(cfg, added_data.dislike_training, initial_data.dislike_cv);
adi_mvpa_sensorspace(like, dislike, time) % original, d.h. unveränderte Mvpa-toolbox
 
%% searchlight: 5.6.19

[like, dislike] = split_data(session, 'cv')

comp1 = [0.05 0.15];
adi_mvpa_searchlight_sensorspace(like, dislike, comp1) 
comp2 = [0.15 0.5];
adi_mvpa_searchlight_sensorspace(like, dislike, comp2) 
comp3 = [0.5 0.8];
adi_mvpa_searchlight_sensorspace(like, dislike, comp3) 


    
% for p = 1:length(session.trial)
%     temp = session.trial{1,p};
%     session.data(p,:,:) = temp;
%     clear temp
% end
%    session = rmfield(session, 'trial');





%
[session_sorted] = regroup_session(session);
adi_sanity_check_sensorspace(session)
adi_stats_sensorspace(session) % hier noch vervollstänidgen
% 
group_path = ['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\svm\' balldesign filesep];
SVM_Guggenmos_adapted (session, group_path, freq)
% % SVM_Guggenmos_virtsens (session, group_path, freq)

%% 17.6. crossvalidation leave one out:

balldesign = 'gbf';
load (['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\new_data\' balldesign '_ball\zscore\' 'sensordata_all_subj_interp_zscore.mat'])

session = sensordata_all_subj_zscore(1);
   for k = 2:length(sensordata_all_subj_zscore)
       session.trial = cat(2,session.trial, sensordata_all_subj_zscore(k).trial);
       session.time = cat(2,session.time, sensordata_all_subj_zscore(k).time);
       session.response_label = cat(2,session.response_label, sensordata_all_subj_zscore(k).response_label);
       session.balldesign_short = cat(2,session.balldesign_short, sensordata_all_subj_zscore(k).balldesign_short);
       session.subject = cat(2,session.subject, sensordata_all_subj_zscore(k).subject);
   end
   

for k=1:length(session.response_label)
    switch session.response_label{k}
        case 'like' % 'Volley'
            session.labels(k) = 1;
        case 'Neu_Like' % 'Volley'
            session.labels(k) = 1;
        case 'dislike' % 'Space'
            session.labels(k) = 2;
        case 'Neu_Dislike' % 'Space'
            session.labels(k) = 2;
    end
end
clear sensordata_all_subj_zscore

adi_mvpa_sensorspace_cv_leave_subj_out(session) 



%% compute lcmv Beamformer and virtual sensors for each subject- 
% anderes Zeitintervall für Covariance!!


% calculations of spatial filter, multiplies it with avgdata 

clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 24%:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep 'MEG\sourcespace\']);
     condition = {'like', 'dislike', 'dontcare'};
   
     adi_virtsens_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'bp1_45Hz', condition)
   
   
end

%% single balldesign: select, z-normalize and compute SVM (from virtual sensors):

clear
path2excelsheet = 'U:\My Documents\MATLAB\eigene_Skripte\adi_project\excelsheets\subjects.xlsx';
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
[delete_run] = adi_boolean_operations(path2excelsheet, subject_list);

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];
% anpassen:
freq = 'bp1_45Hz';
balldesign = 'gbf';

load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'response', [], 'response_label', [], 'balldesign', [], 'subject', []);

for i=1:length(subject_list)
    
    path2virtsens = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' ];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' balldesign filesep];
    [virtsens_all_subj] = adi_vs_single_ball(virtsens_all_subj, path2virtsens, dir_out, subject_list(i).name, pattern, trigger, freq, balldesign, i, delete_run);
end

% z-transformation:
[virtsens_all_subj_zscore] = adi_ztrans_sensorspace(virtsens_all_subj); %evtl. function noch umbenennen
clearvars virtsens_all_subj

for k=1:length(virtsens_all_subj_zscore)
    for p=1:length(virtsens_all_subj_zscore(k).trial)
        virtsens_all_subj_zscore(k).subject_name{p} = virtsens_all_subj_zscore(k).subject;
        virtsens_all_subj_zscore(k).subject_num{p} = virtsens_all_subj_zscore(k).subject(end-1:end);
    end
end

[virtsens_all_subj_zscore] = adi_additionalMEG_interpolation(virtsens_all_subj_zscore, 'virtsens');

session = virtsens_all_subj_zscore(1);
   for k = 2:length(virtsens_all_subj_zscore)
       session.trial = cat(2,session.trial, virtsens_all_subj_zscore(k).trial);
       session.time = cat(2,session.time, virtsens_all_subj_zscore(k).time);
       session.response = cat(2,session.response, virtsens_all_subj_zscore(k).response);
       session.response_label = cat(2,session.response_label, virtsens_all_subj_zscore(k).response_label);
       session.balldesign = cat(2,session.balldesign, virtsens_all_subj_zscore(k).balldesign);
       session.subject_name = cat(2,session.subject_name, virtsens_all_subj_zscore(k).subject_name);
       session.subject_num = cat(2,session.subject_num, virtsens_all_subj_zscore(k).subject_num);
   end
    
for p = 1:length(session.trial)
    temp = session.trial{1,p};
    session.data(p,:,:) = temp;
    clear temp 
end
%    session = rmfield(session, 'trial');

for k=1:length(session.response_label)
    switch session.response_label{k}
        case 'like' % 'Volley'
            session.labels(k) = 1;
        case 'Neu_Like' % 'Volley'
            session.labels(k) = 1;
        case 'dislike' % 'Space'
            session.labels(k) = 2;
        case 'Neu_Dislike' % 'Space'
            session.labels(k) = 2;
    end
end
 
[session] = regroup_session(session);    
% adi_sanity_check_sensorspace(sessions)
% group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\gbf_ball\svm\';
% SVM_Guggenmos_virtsens (session, group_path, freq)
% SVM_Guggenmos_adapted2 (session, group_path, freq)
% SVM_Guggenmos_adapted_sourcespace (session, group_path, freq) % ohne noise whitening, funktioniert aber schlechter

%% source statistic

group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\gbf_ball\stats\';
time = [0 1]; % in s
adi_stats_sourcespace(session, time, group_path, 'gbf') 

%% source statistic über mehrere Bälle zur Hypothesengenerierung

grouppath = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\';
% balldesigns = {'gbf'; 'gbv' ; 'ggs'; 'ggf'};
% balldesigns_names = {'yellow_blue_soccerball' 'yellow_blue_volley' 'grey_green_space' 'grey_green_soccerball'};
% balldesigns = {'rwf'; 'gbs' ; 'rws'; 'rwv'; 'ggv'};
% balldesigns_names = {'red_white_soccerball' 'yellow_blue_space' 'red_white_space' 'red_white_volley' 'grey_green_volley' };

balldesigns = {'ggv'; 'rwv' ; 'rwf'};
balldesigns_names = {'grey_green_volley' 'red_white_volley' 'red_white_soccerball'};


time=[0 1];
% sign_clusters = {'Angular_R'; 'Occipital_Mid_L'; 'Occipital_Mid_R'; 'Parietal_Sup_R'; 'Parietal_Sup_L'; 'Precuneus_R'; 'Temporal_Mid_R'; 'Temporal_Inf_R'; 'Temporal_Sup_R'}; % nr. 66, 86, 90, 82

adi_stats_sourcespace_multiple_balldesigns(time, grouppath, balldesigns, balldesigns_names, []) 

%%
grouppath = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\';
balldesigns = {'rwf'; 'gbs' ; 'rws'; 'rwv'; 'ggv'};
balldesigns_names = {'red_white_soccerball' 'yellow_blue_space' 'red_white_space' 'red_white_volley' 'grey_green_volley' };

time = [0.003 0.25] ; % in s
sign_clusters = {'Angular_R'; 'Occipital_Mid_L'; 'Occipital_Mid_R'; 'Parietal_Sup_R'; 'Parietal_Sup_L'; 'Precuneus_R'; 'Temporal_Mid_R'; 'Temporal_Inf_R'; 'Temporal_Sup_R'}; % nr. 66, 86, 90, 82
adi_stats_sourcespace_multiple_balldesigns(time, grouppath, balldesigns, balldesigns_names, sign_clusters) 

%%

grouppath = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\';
balldesigns = {'rwf'; 'gbs' ; 'rws'; 'rwv'; 'ggv'};
balldesigns_names = {'red_white_soccerball' 'yellow_blue_space' 'red_white_space' 'red_white_volley' 'grey_green_volley' };

time = [0.003 0.25] ; % in s
sign_clusters = {'Angular_R'; 'Occipital_Mid_L'; 'Occipital_Mid_R'; 'Parietal_Sup_R'; 'Parietal_Sup_L'; 'Precuneus_R'; 'Temporal_Mid_R'; 'Temporal_Inf_R'; 'Temporal_Sup_R'}; % nr. 66, 86, 90, 82
adi_stats_sourcespace_multiple_balldesigns(time, grouppath, balldesigns, balldesigns_names, sign_clusters) 



%% source statistic über mehrere Bälle zur Hypothesengenerierung

grouppath = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\favourites\';
% like:
balldesigns = {'ggv'; 'rwv' ; 'rwf'};
balldesigns_names = {'grey_green_volley' 'red_white_volley' 'red_white_soccerball'};

%dislike
balldesigns = {'gbs'; 'ggs' ; 'rws'};
balldesigns_names = {'yellow_blue_space' 'grey_green_space' 'red_white_space'};
time=[0.0 1];
time=[0.06 0.3];
sign_clusters = {'Angular_R'; 'Angular_L'; 'Occipital_Mid_L'; 'Occipital_Mid_R'; 'Parietal_Sup_R'; ...
    'Parietal_Sup_L'; 'Precuneus_R'; 'Precuneus_L'; 'Temporal_Mid_R'; 'Temporal_Mid_L'; 'Temporal_Inf_R'; ...
    'Temporal_Inf_L'; 'Temporal_Sup_R'; 'Temporal_Sup_L'; 'Calcarine_R'; 'Calcarine_L'}; % nr. 66, 86, 90, 82

adi_stats_sourcespace_favourites(time, grouppath, balldesigns, balldesigns_names, [], 'dislike') 

%% rms grand avg für alle Bälle, um Zeitintervall der Komponenten festzulegen: (besser im sensor space)
% 21.3.2019
clear
path2balldesigns = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\';
balldesigns_short = {'gbf'; 'gbs'; 'gbv'; 'ggf'; 'ggv'; 'ggs' ; 'rwf';'rws'; 'rwv'};
balldesigns_long = {'yellow_blue_soccerball'; 'yellow_blue_space'; 'yellow_blue_volley'; 'grey_green_soccerball'; 'grey_green_volley'; 'grey_green_space' ; 'red_white_soccerball';'red_white_space'; 'red_white_volley'};

adi_grandavg_timeintervall(path2balldesigns, balldesigns_short, balldesigns_long)



%% wende Verknüpfungsregeln an, um trials einzelner Probanden rauszuwerfen
clear
path2excelsheet = 'U:\My Documents\MATLAB\eigene_Skripte\adi_project\excelsheets\subjects.xlsx';
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
[delete_run] = adi_boolean_operations(path2excelsheet, subject_list);

grandavg = []; 
[grandavg] = adi_grandavg_sensorspace(subject_list, grandavg, delete_run);
save('grandavg', 'grandavg')
% 
grandavg_like = grandavg.like(1).avg;

for k=2:length(grandavg.like)
    grandavg_like = (grandavg_like+grandavg.like(k).avg)/2;
end

grandavg_dislike = grandavg.dislike(1).avg;

for k=2:length(grandavg.dislike)
    grandavg_dislike = (grandavg_dislike+grandavg.dislike(k).avg)/2;
end

time=cleanMEG_interp.time{1,1};
figure
plot(time(1018:2500), grandavg_like(:,1018:2500));
axis tight
title('grand avg like')

figure
plot(time(1018:2600), grandavg_dislike(:,1018:2600));
axis tight
title('grand avg dislike')
% comp 1: 50-160ms; 
% comp 2: 160-500ms
% comp 3: 500-800ms

avg=(grandavg_like+grandavg_dislike)/2;
figure
plot(time(1018:2600), avg(:, 1018:2600));
title('avg_like_dislike')

grandavg_dislike(230,:)=[];
grandavg_like(230,:)=[];



find(grandavg_dislike(:, 1609)==max(grandavg_dislike))

[x,y]=find(grandavg_dislike(:, 1467)==max(grandavg_dislike(:, 1467)))







data_like.trial{1,1} = grandavg_like;
data_like.label = cleanMEG_interp.label;
data_like.time{1,1} = time;



figure
plot(time(1018:end), data_like.trial{1,1}(:, 1018:end))
axis tight
title('like');

% 
% poststim=1:round(1017.25*2);
% prestim=1:round(1017*1.5);
% 
% 
% x=1:3561;
% y=x/1017.25;
% z=(1:round(1017.25/2);
% 
% time=cleanMEG_interp.time{1,1};
% 
% 
% figure
% plot(time, grandavg_dislike)



%% sensor space grand avg
clear 
mainpath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\'; 
path2data = 'MEG_analysis\noisereduced\1_95Hz\02_interpolated\';
all_subject = dir(mainpath);
all_subject(1:2)=[];

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];

balldesign = 'ggs';
delete_subjects = {'nl_adi_20'; 'nl_adi_29'; 'nl_adi_34'};
for i=1:length(subject_list)
    ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
end
subject_list(find(ind),:) = [];
clear ind

delete_runs={'nl_adi_04', 3; 'nl_adi_05', [2];  'nl_adi_10', [1]; 'nl_adi_15', 2; 'nl_adi_17' 1; 'nl_adi_22', [2 3]; 'nl_adi_23', [1 2]; 'nl_adi_28', [1]; 'nl_adi_33', 1} ; 


grandavg.like=[];
grandavg.dislike = [];

delete_runs = struct();


load('U:\My Documents\MATLAB\atlas_source_indices.mat');
sensordata_single_ball_all_subj = struct('trial', [], 'time', [], 'response', [], 'response_label', [], 'balldesign', [], 'subject', []);

for i=1:length(subject_list)
    
    path2virtsens = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' ];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' balldesign filesep];
    [sensordata_single_ball_all_subj] = adi_sensor_single_ball(sensordata_single_ball_all_subj, path2virtsens, dir_out, subject_list(i).name, pattern, trigger, freq, balldesign, delete_runs, i);
end

for k=1:length(all_subject)
    
    [grandavg] = adi_grandavg_sensorspace([mainpath all_subject(k).name filesep path2data], grandavg);
    
    
    
end


grandavg_like = like_avg(1).avg;

for k=2:length(like_avg)
    grandavg_like= grandavg_like+like_avg(k).avg;
  
end
grandavg_like=grandavg_like/length(like_avg);

figure
plot(time, grandavg_like)
title('grandavg_like_gbf')


grandavg_dislike = dislike_avg(1).avg;

for k=2:length(dislike_avg)
    grandavg_dislike= grandavg_dislike+dislike_avg(k).avg;
  
end
grandavg_dislike=grandavg_dislike/length(dislike_avg);

figure
plot(time, grandavg_dislike)
title('grandavg_dislike_gbf')

grand_avg_both_cond = (grandavg_like+grandavg_dislike)/2;
figure
plot(time, grand_avg_both_cond)
title('grand_avg_both_cond')

figure
plot(time, grandavg_like)
hold on
plot(time, grandavg_dislike)








time=virtsens_all_subj_zscore(1).time{1,1};
mean_like= mean(like_avg(1).avg);

for k=2:length(like_avg)
    mean_like(k,:)= mean(like_avg(k).avg);
  
end

grand_avg_like=mean(mean_like);
%%




mean_dislike= mean(dislike_avg(1).avg);

for k=2:length(dislike_avg)
    mean_dislike(k,:)= mean(dislike_avg(k).avg);
  
end

grand_avg_dislike=mean(mean_dislike);


figure
plot(time, grand_avg_like)
hold on
plot(time, grand_avg_dislike)
legend({'like', 'dislike'});


figure
plot(time, mean_like)
hold on
plot(time, grand_avg_dislike)
legend({'like', 'dislike'});

adi_stats_mixed_models(session, path, balldesign)

%% 

load 'U:\My Documents\MATLAB\atlas_clusterstatistic_ohne_cerebellum_thalamus_basalganglien.mat';

roi_num = unique(cell2mat(atlas_cluster(:,1)));

for k=1:length(roi_num)
   
    ind_rois = find(cell2mat(atlas_cluster(:,1))==roi_num(k));
    for i = 1:length(data_like.trial)
%         data_like_roi{i}(k,:) = mean(data_like.trial{i}(ind_rois, :)); 
        data_dislike_roi{i}(k,:) = mean(data_dislike.trial{i}(ind_rois, :)); 
    end
    label(k,1) = atlas_cluster(ind_rois(1),2);
    clear ind_rois,
    
end

data_like.trial = data_like_roi;
data_dislike.trial = data_dislike_roi;
clear data_like_roi data_dislike_roi

data_like.label = label;
data_dislike.label = label;

% time-frequency analysis
cfg            = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';% oder: 'tfr'
cfg.taper        = 'hanning';
cfg.keeptrials = 'yes';
cfg.toi        = data_like.time{1}(1:2:end); % half of temporal resolution
cfg.foi        = [2:1:30];
cfg.t_ftimwin =  ones(length(cfg.foi),1).*0.5;
% cfg.t_ftimwin =  [repmat(1.5, 1, length(cfg.foi(cfg.foi<16)))./cfg.foi(cfg.foi<16) repmat(1, 1, length(cfg.foi(cfg.foi>=16)))./cfg.foi(cfg.foi>=16)];
% cfg.t_ftimwin = [2./(cfg.tapsmofrq(cfg.foi<16)*2),repmat(0.25,1,length(cfg.foi(cfg.foi>=16)))];

cfg.trials = 'all'; 
TFR_like = ft_freqanalysis(cfg, data_like);
TFR_dislike   = ft_freqanalysis(cfg, data_dislike);

cfg = [];
cfg.baseline     = [-0.5 0];
cfg.baselinetype = 'absolute';
cfg.maskstyle    = 'saturation';
cfg.channel      = TFR_like.label(66);

figure;
subplot(211); ft_singleplotTFR(cfg, TFR_like_1); title(['like_' TFR_like.label{66}]);
subplot(212); ft_singleplotTFR(cfg, TFR_dislike_1); title(['dislike_' TFR_like.label{66}]);

% 
% cfg = [];
% cfg.channel          = {'MEG'};
% cfg.latency          = 'all';
% cfg.frequency        = 10;
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'ft_statfun_indepsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan        = 2;
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 500;
% cfg.clustercritval = 0.05; 
% 
% load('E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\neighbours_rois.mat')
% 
% % prepare_neighbours determines what sensors may form clusters
% cfg_neighb.method    = 'distance';
% cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, data_like);
% 
% cfg.neighbours = neighbours;
% 
% design = zeros(1,size(TFR_like.powspctrm,1) + size(TFR_dislike.powspctrm,1));
% design(1,1:size(TFR_like.powspctrm,1)) = 1;
% design(1,(size(TFR_like.powspctrm,1)+1):(size(TFR_like.powspctrm,1)+...
% size(TFR_dislike.powspctrm,1))) = 2;
% 
% cfg.design           = design;
% cfg.ivar             = 1;
% 
% [stat] = ft_freqstatistics(cfg, TFR_like, TFR_dislike);



beta_like = data_like;
beta_dislike = data_dislike;

% create HGP as empty timelock structure with same dimensions as ERP, values will be filled in in the next steps
pb_1_30_like = rmfield(data_like,{'trial'});
pb_1_30_dislike   = rmfield(data_dislike,{'trial'});

% correct for the 1/f dropoff
freqcorr = reshape(TFR_like.freq.^2,[1 1 length(TFR_like.freq)]); %this vector accounts for the 1/f dropoff
% use repmat to create a matrix the same size as the TFR data
freqcorr_like = repmat(freqcorr,[size(TFR_like.powspctrm,1) size(TFR_like.powspctrm,2) 1 length(TFR_like.time)]);
freqcorr_dislike   = repmat(freqcorr,[size(TFR_dislike.powspctrm,1) size(TFR_dislike.powspctrm,2) 1 length(TFR_dislike.time)]);

% multiply data with freqcorr matrix and average over frequencies
pb_1_30_like.trial = squeeze(nanmean(TFR_like.powspctrm(:,:,:,:) .* freqcorr_like,3));
pb_1_30_dislike.trial   = squeeze(nanmean(TFR_dislike.powspctrm(:,:,:,:) .* freqcorr_dislike,3));

% calculate mean and variance
pb_1_30_like.avg = squeeze(nanmean(pb_1_30_like.trial,1));
pb_1_30_like.var = squeeze(nanvar(pb_1_30_like.trial,1));
pb_1_30_dislike.avg   = squeeze(nanmean(pb_1_30_dislike.trial,1));
pb_1_30_dislike.var   = squeeze(nanvar(pb_1_30_dislike.trial,1));


% baseline correction
cfg          = [];
cfg.baseline = [-.5 0];

pb_1_30_like_bl = ft_timelockbaseline(cfg,pb_1_30_like);
pb_1_30_dislike_bl   = ft_timelockbaseline(cfg,pb_1_30_dislike);

cfg           = [];
cfg.parameter = 'avg';
% cfg.xlim      = [-.3 .6];
cfg.channel   = TFR_like.label{78}; % other responsive channels: 'EEG PT_03-REF', 'EEG PT_04-REF', 'EEG IO_02-REF', 'EEG IO_04-REF', 'EEG SO_01-REF', 'EEG SO_02-REF''EEG SO_03-REF'

figure, ft_singleplotER(cfg, pb_1_30_like,pb_1_30_dislike)

cfg                  = [];
cfg.latency          = [0 .6];
cfg.parameter        = 'trial';
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.neighbours       = []; % no spatial information is exploited for statistical clustering
cfg.numrandomization = 500;
cfg.statistic        = 'indepsamplesT'; % idependent samples test for statistical testing on the single-trial level
cfg.channel          = 'all';
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.design           = [ones(1,size(ERP_object_bl.trial,1)), 2*ones(1,size(ERP_face_bl.trial,1))];

stats_ERP = ft_timelockstatistics(cfg,ERP_object_bl,ERP_face_bl);
stats_HGP = ft_timelockstatistics(cfg,HGP_object_bl,HGP_face_bl);

% look for significant channels in ERP stats
[chans time] = find(stats_ERP.mask)

% none of the channels contain significant differences between conditions

% look for significant channels in HGP stats
[chans time] = find(stats_HGP.mask);
chans = unique(chans)

% first, average over trials, otherwise we'll have problems with ft_singleplotER
cfg = [];
HGP_object_bl = ft_timelockanalysis(cfg, HGP_object_bl)
HGP_face_bl   = ft_timelockanalysis(cfg, HGP_face_bl)

% add statistical mask to data
time_idx = find(HGP_object_bl.time == stats_HGP.time(1)) : find(HGP_object_bl.time == stats_HGP.time(end)); % find indices of timepoints in data corresponding to timepoints in stats
HGP_object_bl.mask = false(size(HGP_object_bl.avg));
HGP_object_bl.mask(:,time_idx) = stats_HGP.mask;

% plot the ERP traces
cfg               = [];
cfg.parameter     = 'avg';
cfg.xlim          = [-.3 .6];
cfg.maskparameter = 'mask';

% loop over significant channels and plot
for ichan = 1:length(chans)
    cfg.channel = chans(ichan);
    figure, ft_singleplotER(cfg,HGP_object_bl,HGP_face_bl),
end

% time-frequency analysis
cfg            = [];
cfg.method     = 'tfr';
cfg.output     = 'pow';
cfg.keeptrials = 'yes'; % keep trials for statistics
cfg.foi        = [4:2:40 44:4:100 108:8:200]; % make the frequency spacing broader with higher frequencies
cfg.toi        = -.3:.02:.6;

cfg.trials = find(epoch_data_clean.trialinfo==3); % select 'object' trials
TFR_object = ft_freqanalysis(cfg, epoch_data_clean);

cfg.trials = find(epoch_data_clean.trialinfo==7); % select 'face' trials
TFR_face   = ft_freqanalysis(cfg, epoch_data_clean);

% baseline correction
cfg              = [];
cfg.baseline     = [-.3 .05];
cfg.baselinetype ='relchange';

TFR_object_bl = ft_freqbaseline(cfg, TFR_object);
TFR_face_bl   = ft_freqbaseline(cfg, TFR_face);

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.xlim      = [-.3 .6];
cfg.zlim      = [-10 10];
cfg.channel   = 'EEG IO_03-REF'; % other responsive channels: 'EEG PT_03-REF', 'EEG PT_04-REF', 'EEG IO_02-REF', 'EEG IO_04-REF', 'EEG SO_01-REF', 'EEG SO_02-REF''EEG SO_03-REF'

figure, ft_singleplotTFR(cfg,TFR_object_bl)
figure, ft_singleplotTFR(cfg,TFR_face_bl)

fg                  = [];
cfg.latency          = [0 .6];
cfg.parameter        = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.neighbours       = [];
cfg.numrandomization = 500;
cfg.statistic        = 'indepsamplesT';
cfg.channel          = 'all';
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';
cfg.clusteralpha     = 0.05;
cfg.design           = [ones(1, size(TFR_object_bl.powspctrm,1)), 2*ones(1,size(TFR_face_bl.powspctrm,1))];
cfg.ivar             = 1;

stats_TFR = ft_freqstatistics(cfg,TFR_object_bl,TFR_face_bl);
Finally, we will plot the masked t-values from significant channels of the statistical contrast.

% find significant channels
sig_tiles = find(stats_TFR.mask); % find significant time-frequency tiles
[chan freq time] = ind2sub(size(stats_TFR.mask),sig_tiles); % transform linear indices to subscript to extract significant channels, timepoints and frequencies
chan = unique(chan);

% plot TFRs
cfg               = [];
cfg.parameter     = 'stat';
cfg.maskparameter = 'mask';
cfg.maskalpha     = .4; % opacity value for non-significant parts
cfg.zlim          = [-4 4];
cfg.colorbar      = 'yes';

% loop over channel sets and plot
for ichan = 1:length(chan)
    cfg.channel = chan(ichan);
    figure, ft_singleplotTFR(cfg,stats_TFR),
end

%% statistische Analyse bei adi 05 und ggs ball; 

clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];
% anpassen:
freq = 'bp1_45Hz';
delete_runs={} ; 

balldesign = 'ggv';

load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'response', [], 'response_label', [], 'balldesign', [], 'subject', []);


for i=3%1:length(subject_list)
    
    path2virtsens = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' ];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' balldesign filesep];
    [virtsens_all_subj] = adi_vs_single_ball(virtsens_all_subj, path2virtsens, dir_out, subject_list(i).name, pattern, trigger, freq, balldesign, delete_runs, i);
end

virtsens_all_subj(1:i-1)=[];
% z-transformation:
[virtsens_all_subj_zscore] = adi_ztrans_sensorspace(virtsens_all_subj); %evtl. function noch umbenennen
clearvars virtsens_all_subj

[virtsens_all_subj_zscore] = adi_additionalMEG_interpolation(virtsens_all_subj_zscore, 'virtsens');

session = virtsens_all_subj_zscore(1);
   for k = 2:length(virtsens_all_subj_zscore)
       session.trial = cat(2,session.trial, virtsens_all_subj_zscore(k).trial);
       session.time = cat(2,session.time, virtsens_all_subj_zscore(k).time);
       session.response = cat(2,session.response, virtsens_all_subj_zscore(k).response);
       session.response_label = cat(2,session.response_label, virtsens_all_subj_zscore(k).response_label);
       session.balldesign = cat(2,session.balldesign, virtsens_all_subj_zscore(k).balldesign);
   end
   
   
session = virtsens_all_subj_zscore(1);
   for k = 2:length(virtsens_all_subj_zscore)
       session.trial = cat(2,session.trial, virtsens_all_subj_zscore(k).trial);
       session.time = cat(2,session.time, virtsens_all_subj_zscore(k).time);
       session.response = cat(2,session.response, virtsens_all_subj_zscore(k).response);
       session.response_label = cat(2,session.response_label, virtsens_all_subj_zscore(k).response_label);
       session.balldesign = cat(2,session.balldesign, virtsens_all_subj_zscore(k).balldesign);
   end
    
for p = 1:length(session.trial)
    temp = session.trial{1,p};
    session.data(p,:,:) = temp;
    clear temp
end
%    session = rmfield(session, 'trial');

for k=1:length(session.response_label)
    switch session.response_label{k}
        case 'like' % 'Volley'
            session.labels(k) = 1;
        case 'Neu_Like' % 'Volley'
            session.labels(k) = 1;
        case 'dislike' % 'Space'
            session.labels(k) = 2;
        case 'Neu_Dislike' % 'Space'
            session.labels(k) = 2;
    end
end
 
[session_ggv] = regroup_session(session);    

group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\rwf_ball\nl_adi_06\stats\';
time = [0 1]; % in s
% adi_stats_sourcespace(session_28, time, group_path, 'rwf') 





   %% select first run and compute svm
   
clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];
% anpassen:
freq = 'bp1_45Hz';
balldesign = 'gbf';
subject_list([1 2 16 18 24 25],:) = [];
delete_runs={'nl_adi_06', '2'; 'nl_adi_17', '1'; 'nl_adi_24', '2'; 'nl_adi_25', '1'};

% subject_list (2:29,:) = [];

load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'response', [], 'response_label', [], 'balldesign', []);

for i=1:length(subject_list)
    path2virtsens = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' ];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' balldesign filesep];
    [virtsens_all_subj] = adi_vs_single_ball(virtsens_all_subj, path2virtsens, dir_out, subject_list(i).name, pattern, trigger, freq, balldesign, delete_runs, i);
end

session = virtsens_all_subj(1);
   for k = 2:length(virtsens_all_subj)
       session.trial = cat(2,session.trial, virtsens_all_subj(k).trial);
       session.time = cat(2,session.time, virtsens_all_subj(k).time);
       session.response = cat(2,session.response, virtsens_all_subj(k).response);
       session.response_label = cat(2,session.response_label, virtsens_all_subj(k).response_label);
       session.balldesign = cat(2,session.balldesign, virtsens_all_subj(k).balldesign);
   end
    
   clearvars virtsens_all_subj
   
for p = 1:length(session.trial)
    temp = session.trial{1,p};
    session.data(p,:,:) = temp;
    clear temp
end
session = rmfield(session, 'trial');

for k=1:length(session.response_label)
    switch session.response_label{k}
        case 'like' % 'Volley'
            session.labels(k) = 1;
        case 'Neu_Like' % 'Volley'
            session.labels(k) = 1;
        case 'dislike' % 'Space'
            session.labels(k) = 2;
        case 'Neu_Dislike' % 'Space'
            session.labels(k) = 2;
    end
end
 
   group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\gbf_ball\svm\';
   SVM_Guggenmos_virtsens (session, group_path, freq)

  
%% select and compute SVM from single ball (from virtual sensors): Covariance Test

clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];
% anpassen:
freq = 'bp1_45Hz';
balldesign = 'gbf';
subject_list([1 2 16 18 24 25],:) = [];
delete_runs={'nl_adi_06', '2'; 'nl_adi_17', '1'; 'nl_adi_24', '2'; 'nl_adi_25', '1'};

% subject_list (2:29,:) = [];

load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'response', [], 'response_label', [], 'balldesign', []);

for i=1:length(subject_list)
    path2virtsens = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\Covariance_Test\virtsens\' ];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\Covariance_Test\virtsens\' balldesign filesep];
    [virtsens_all_subj] = adi_vs_single_ball(virtsens_all_subj, path2virtsens, dir_out, subject_list(i).name, pattern, trigger, freq, balldesign, delete_runs, i);
end

session = virtsens_all_subj(1);
   for k = 2:length(virtsens_all_subj)
       session.trial = cat(2,session.trial, virtsens_all_subj(k).trial);
       session.time = cat(2,session.time, virtsens_all_subj(k).time);
       session.response = cat(2,session.response, virtsens_all_subj(k).response);
       session.response_label = cat(2,session.response_label, virtsens_all_subj(k).response_label);
       session.balldesign = cat(2,session.balldesign, virtsens_all_subj(k).balldesign);
   end
    
   clearvars virtsens_all_subj
   
for p = 1:length(session.trial)
    temp = session.trial{1,p};
    session.data(p,:,:) = temp;
    clear temp
end
session = rmfield(session, 'trial');

for k=1:length(session.response_label)
    switch session.response_label{k}
        case 'like' % 'Volley'
            session.labels(k) = 1;
        case 'Neu_Like' % 'Volley'
            session.labels(k) = 1;
        case 'dislike' % 'Space'
            session.labels(k) = 2;
        case 'Neu_Dislike' % 'Space'
            session.labels(k) = 2;
    end
end
 
   group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\gbf_ball\Covariance_Test\svm\';
   SVM_Guggenmos_virtsens (session, group_path, freq)

%% compute virtual sensors mit noisenormalization for each subject and do a group analysis: like, dislike, dontcare  - 

clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];
freq = 'bp1_45Hz';
condition.like = 1;
condition.dislike = 2;
condition.dontcare = 3;
load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'balldesign', [], 'triggerlabel', [], 'label', [], 'fsample', [], 'response_label', []);

for i=3%:length(subject_list)
    
    vs_allRuns = struct();
    path2sensordata = [path2subj subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
    path2spatialfilter = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\spatialfilter\'];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtual_sensors\'];
    
    if ~exist([dir_out 'virtsens_ns\virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'file')
    if 1 == strcmp (subject_list(i).name, 'nl_adi_07')
        [vs_allRuns] = adi_regroup_subjects(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(2));
        [vs_allRuns] = adi_regroup_subjects(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(3));
    else
        [vs_allRuns] = adi_regroup_subjects(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(1));
        [vs_allRuns] = adi_regroup_subjects(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(2));
        [vs_allRuns] = adi_regroup_subjects(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(3)); 
    end
    run = fields(vs_allRuns);
    
    switch size(fields(vs_allRuns),1)
        case 3 
            vs_allRuns_appended.trial = [vs_allRuns.run1.trial vs_allRuns.run2.trial vs_allRuns.run3.trial];
            vs_allRuns_appended.time =  [vs_allRuns.run1.time vs_allRuns.run2.time vs_allRuns.run3.time];
            vs_allRuns_appended.balldesign =  [vs_allRuns.run1.balldesign vs_allRuns.run2.balldesign vs_allRuns.run3.balldesign];
            vs_allRuns_appended.triggerlabel =  [vs_allRuns.run1.triggerlabel vs_allRuns.run2.triggerlabel vs_allRuns.run3.triggerlabel];         
            vs_allRuns_appended.response_label =  [vs_allRuns.run1.response_label vs_allRuns.run2.response_label vs_allRuns.run3.response_label];         
        case 2
            if 1 == strcmp (subject_list(i).name, 'nl_adi_07')
                vs_allRuns_appended.trial = [vs_allRuns.run2.trial vs_allRuns.run3.trial];
                vs_allRuns_appended.time =  [vs_allRuns.run2.time vs_allRuns.run3.time];
                vs_allRuns_appended.balldesign =  [vs_allRuns.run2.balldesign vs_allRuns.run3.balldesign];
                vs_allRuns_appended.triggerlabel =  [vs_allRuns.run2.triggerlabel vs_allRuns.run3.triggerlabel];   
                vs_allRuns_appended.response_label =  [vs_allRuns.run2.response_label vs_allRuns.run3.response_label];            
            else
                vs_allRuns_appended.trial = [vs_allRuns.run1.trial vs_allRuns.run2.trial];
                vs_allRuns_appended.time =  [vs_allRuns.run1.time vs_allRuns.run2.time];
                vs_allRuns_appended.balldesign =  [vs_allRuns.run1.balldesign vs_allRuns.run2.balldesign];
                vs_allRuns_appended.triggerlabel =  [vs_allRuns.run1.triggerlabel vs_allRuns.run2.triggerlabel];   
                vs_allRuns_appended.response_label =  [vs_allRuns.run1.response_label vs_allRuns.run2.response_label];
            end
        end
        
    vs_allRuns_appended.label =  vs_allRuns.run2.label;
    vs_allRuns_appended.fsample =  vs_allRuns.run2.fsample;
    clearvars vs_allRuns
    
    % sanity check Nr. 3:
   
    path_avg = [dir_out '\sanity_check\avg\'];
    if ~exist ( path_avg, 'dir')
        mkdir (path_avg)
    end
    cfg = [];
    avg = ft_timelockanalysis(cfg, vs_allRuns_appended);
    figure;
    plot(avg.time, avg.avg)
    savefig([dir_out filesep 'sanity_check\avg\avg_virtsens_allRuns_' freq '.fig'])
    
    % Safe freq spectrum:
    path_spect = [dir_out '\sanity_check\freq_spectrum\'];
    if ~exist ( path_spect, 'dir')
        mkdir (path_spect)
    end
    
    sRate = vs_allRuns_appended.fsample;
    [FourRef,Fref]=fftBasic(avg.avg, round(sRate));
    figure;
    plot(FourRef, Fref)
    savefig([path_spect 'avg_virtsens_allRuns_freq_spectrum_' freq '.fig'])
    close
    
    pathAppended = [dir_out '\virtsens_ns\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end

    save ([pathAppended 'virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'vs_allRuns_appended', '-v7.3');
    
%     else
%         load([dir_out '\virtsens_ns\virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'vs_allRuns_appended')
    end
%     if ~isequal(length(vs_allRuns_appended.trial), length(vs_allRuns_appended.response_label))
%        error([subject_list(i).name ' : trial and resonse_label length do not match']) 
%     end
%     virtsens_all_subj(i) = vs_allRuns_appended;
%     clearvars vs_allRuns_appended
end

%     
%    session = virtsens_all_subj(1);
%    for k = 2:length(virtsens_all_subj)
%        session.trial = cat(2,session.trial, virtsens_all_subj(k).trial);
%        session.time = cat(2,session.time, virtsens_all_subj(k).time);
%        session.balldesign = cat(2,session.balldesign, virtsens_all_subj(k).balldesign);
%        session.triggerlabel = cat(2,session.triggerlabel, virtsens_all_subj(k).triggerlabel);
%        session.response_label = cat(2,session.response_label, virtsens_all_subj(k).response_label);
%    end
%     
%    clearvars virtsens_all_subj
%    
% for i = 1:length(session.trial)
%     temp = session.trial{1,i};
%     session.data(i,:,:) = temp;
%     clear temp
% end
%    session = rmfield(session, 'trial');
%    session.label = session.response_label;
%     
%    group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\pattern_comparison\comparison_response_labels\guggenmos_svm_results\';
%    svm_guggenmos (session, group_path, freq, condition, atlas_downsampled)
%     
    

%% compute svm per subject:

clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
% virtsens_all_subj = struct('fsample', [], 'trial', [], 'time', [], 'label', [], 'cfg', [], 'balldesign', [], 'response_label', [], 'triggerlabel', []);
filter = 'theta';

for i=3%:length(subject_list)
    subject_file_vs = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtual_sensors\virtsens_ns\virtsens_ns_all_conditions_allRuns_bp1_95Hz.mat'];
    [session] = adi_filter_virtsens(subject_file_vs, filter, subject_list(i).name);
   
%     virtsens_all_subj(i) = vs_filtered;
%     clearvars vs_allRuns_appended
end


%  session = virtsens_all_subj(1);
%    for k = 2:length(virtsens_all_subj)
%        session.trial = cat(2,session.trial, virtsens_all_subj(k).trial);
%        session.time = cat(2,session.time, virtsens_all_subj(k).time);
%        session.balldesign = cat(2,session.balldesign, virtsens_all_subj(k).balldesign);
%        session.triggerlabel = cat(2,session.triggerlabel, virtsens_all_subj(k).triggerlabel);
%        session.response_label = cat(2,session.response_label, virtsens_all_subj(k).response_label);
%    end
    
%    clearvars virtsens_all_subj
   
for p = 1:length(session.trial)
    temp = session.trial{1,p};
    session.data(p,:,:) = temp;
    clear temp
end
   session = rmfield(session, 'trial');
   session.label = session.response_label;
    
   subj_path = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtual_sensors\svm\'];
   SVM_Guggenmos_virtsens (session, subj_path, filter)
    
%% source analysis like vs. dislike single subject:

clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
filter = 'bp1_45Hz';

for i = 3%:length(ListSubj) 
     path2vol = [path2subj subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\vol\'];
     path2data = [path2subj subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
     path_T1warped_ft = ([path2subj, subject_list(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  subject_list(i).name filesep 'MEG\sourcespace\source_avg\']);
     condition = {'like', 'dislike'};
     adi_source_avg_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'bp1_45Hz', condition)
end

%%
clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
filter = 'bp1_45Hz';

like.run1=source_avg1.run1.like.avg.pow(source_avg1.run1.like.inside==1);
like.run2=source_avg2.run2.like.avg.pow(source_avg2.run2.like.inside==1);
like.run3=source_avg3.run3.like.avg.pow(source_avg3.run3.like.inside==1);
mean_like=(like.run1+like.run2+like.run3)/3;

dislike.run1=source_avg1.run1.dislike.avg.pow(source_avg1.run1.dislike.inside==1);
dislike.run2=source_avg2.run2.dislike.avg.pow(source_avg2.run2.dislike.inside==1);
dislike.run3=source_avg3.run3.dislike.avg.pow(source_avg3.run3.dislike.inside==1);
mean_dislike=(dislike.run1+dislike.run2+dislike.run3)/3;


avg_like = (source_avg1.run1.like.avg.pow(source_avg1.run1.like.inside==1)+source_avg2.run2.like.avg.pow(source_avg2.run2.like.inside==1)+source_avg3.run3.like.avg.pow(source_avg3.run3.like.inside==1))/3;
avg_dislike = (source_avg1.run1.dislike.avg.pow(source_avg1.run1.dislike.inside==1)+source_avg2.run2.dislike.avg.pow(source_avg2.run2.dislike.inside==1)+source_avg3.run3.dislike.avg.pow(source_avg3.run3.dislike.inside==1))/3;
% interpolate:
cfg=[];
cfg.funparameter
ft_sourceplot(cfg, source_avg1.run1.like);
sourcediff = avg_like-avg_dislike;
sum(avg_like)
sum(avg_dislike)
sum(sourcediff)

figure
plot(avg_data_like.time, like)

figure
plot(avg_data_like.time, dislike)

figure
plot(avg_data_like.time, diff)

sum(source_avg1.run1.like.avg.pow(source_avg1.run1.like.inside==1))

%% construct virtual channels of single balldesigns: gbf

clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];
% anpassen:
freq = 'bp1_95Hz';
balldesign = 'gbf';
nsubjects = [3:6 8:15 17:24 26 27 29 30];

condition.like = 1;
condition.dislike = 2;
condition.dontcare = 3;

load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'balldesign', [], 'triggerlabel', [], 'label', [], 'fsample', [], 'response_label', []);


for i=nsubjects
    
    vs_allRuns = struct();
    path2sensordata = [path2subj subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
    path2spatialfilter = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\spatialfilter\'];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtual_sensors\' balldesign filesep];
    
    if ~exist([dir_out 'virtsens_ns\virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'file')
    if 1 == strcmp (subject_list(i).name, 'nl_adi_07')
        [vs_allRuns] = adi_vs_single_ball(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(2), balldesign);
        [vs_allRuns] = adi_vs_single_ball(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(3), balldesign);
    else
        [vs_allRuns] = adi_vs_single_ball(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(1), balldesign);
        [vs_allRuns] = adi_vs_single_ball(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(2), balldesign);
        [vs_allRuns] = adi_vs_single_ball(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject_list, pattern, trigger, freq, num2str(3), balldesign); 
    end
    run = fields(vs_allRuns);
    
    switch size(fields(vs_allRuns),1)
        case 3 
            vs_allRuns_appended.trial = [vs_allRuns.run1.trial vs_allRuns.run2.trial vs_allRuns.run3.trial];
            vs_allRuns_appended.time =  [vs_allRuns.run1.time vs_allRuns.run2.time vs_allRuns.run3.time];
            vs_allRuns_appended.balldesign =  [vs_allRuns.run1.balldesign vs_allRuns.run2.balldesign vs_allRuns.run3.balldesign];
            vs_allRuns_appended.triggerlabel =  [vs_allRuns.run1.triggerlabel vs_allRuns.run2.triggerlabel vs_allRuns.run3.triggerlabel];         
            vs_allRuns_appended.response_label =  [vs_allRuns.run1.response_label vs_allRuns.run2.response_label vs_allRuns.run3.response_label];         
        case 2
            if 1 == strcmp (subject_list(i).name, 'nl_adi_07')
                vs_allRuns_appended.trial = [vs_allRuns.run2.trial vs_allRuns.run3.trial];
                vs_allRuns_appended.time =  [vs_allRuns.run2.time vs_allRuns.run3.time];
                vs_allRuns_appended.balldesign =  [vs_allRuns.run2.balldesign vs_allRuns.run3.balldesign];
                vs_allRuns_appended.triggerlabel =  [vs_allRuns.run2.triggerlabel vs_allRuns.run3.triggerlabel];   
                vs_allRuns_appended.response_label =  [vs_allRuns.run2.response_label vs_allRuns.run3.response_label];            
            else
                vs_allRuns_appended.trial = [vs_allRuns.run1.trial vs_allRuns.run2.trial];
                vs_allRuns_appended.time =  [vs_allRuns.run1.time vs_allRuns.run2.time];
                vs_allRuns_appended.balldesign =  [vs_allRuns.run1.balldesign vs_allRuns.run2.balldesign];
                vs_allRuns_appended.triggerlabel =  [vs_allRuns.run1.triggerlabel vs_allRuns.run2.triggerlabel];   
                vs_allRuns_appended.response_label =  [vs_allRuns.run1.response_label vs_allRuns.run2.response_label];
            end
        end
        
    vs_allRuns_appended.label =  vs_allRuns.run2.label;
    vs_allRuns_appended.fsample =  vs_allRuns.run2.fsample;
    clearvars vs_allRuns
    
    % sanity check Nr. 3:
   
    path_avg = [dir_out '\sanity_check\avg\'];
    if ~exist ( path_avg, 'dir')
        mkdir (path_avg)
    end
    cfg = [];
    avg = ft_timelockanalysis(cfg, vs_allRuns_appended);
    figure;
    plot(avg.time, avg.avg)
    savefig([dir_out filesep 'sanity_check\avg\avg_virtsens_allRuns_' freq '.fig'])
    
    % Safe freq spectrum:
    path_spect = [dir_out '\sanity_check\freq_spectrum\'];
    if ~exist ( path_spect, 'dir')
        mkdir (path_spect)
    end
    
    sRate = vs_allRuns_appended.fsample;
    [FourRef,Fref]=fftBasic(avg.avg, round(sRate));
    figure;
    plot(FourRef, Fref)
    savefig([path_spect 'avg_virtsens_allRuns_freq_spectrum_' freq '.fig'])
    close
    
    pathAppended = [dir_out '\virtsens_ns\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end

    save ([pathAppended 'virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'vs_allRuns_appended', '-v7.3');
    
%     else
%         load([dir_out '\virtsens_ns\virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'vs_allRuns_appended')
    end
%     if ~isequal(length(vs_allRuns_appended.trial), length(vs_allRuns_appended.response_label))
%        error([subject_list(i).name ' : trial and resonse_label length do not match']) 
%     end
    virtsens_all_subj(i) = vs_allRuns_appended;
    clearvars vs_allRuns_appended
end

virtsens_all_subj


%     
%    session = virtsens_all_subj(1);
%    for k = 2:length(virtsens_all_subj)
%        session.trial = cat(2,session.trial, virtsens_all_subj(k).trial);
%        session.time = cat(2,session.time, virtsens_all_subj(k).time);
%        session.balldesign = cat(2,session.balldesign, virtsens_all_subj(k).balldesign);
%        session.triggerlabel = cat(2,session.triggerlabel, virtsens_all_subj(k).triggerlabel);
%        session.response_label = cat(2,session.response_label, virtsens_all_subj(k).response_label);
%    end
%     
%    clearvars virtsens_all_subj
%    
% for i = 1:length(session.trial)
%     temp = session.trial{1,i};
%     session.data(i,:,:) = temp;
%     clear temp
% end
%    session = rmfield(session, 'trial');
%    session.label = session.response_label;
%     
%    group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\pattern_comparison\comparison_response_labels\guggenmos_svm_results\';
%    svm_guggenmos (session, group_path, freq, condition, atlas_downsampled)
%     
   
    %% construct virtual channels, nur erster run, Probanden 1-15, alle Bälle:

clear
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];

pattern.volley.labels = {'rwv'; 'ggv'; 'gbv'};
pattern.space.labels = {'rws'; 'ggs'; 'gbs'};
pattern.soccer.labels = {'rwf'; 'ggf'; 'gbf'};
trigger.labels_short    = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.labels          = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggerchannel  = [102, 104, 106, 108, 4196, 4198, 4200, 4202, 4204];
trigger.eprime          = [102, 104, 106, 108, 101, 103, 105, 107, 109];
% anpassen:
freq = 'bp1_95Hz';
run = 'Run1';
% nsubjects = [3 5 6 8:15 17:18];
nsubjects = 1:18;
condition.like = 1;
condition.dislike = 2;
condition.dontcare = 3;

load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'balldesign', [], 'triggerlabel', [], 'label', [], 'fsample', [], 'response_label', []);


for i=nsubjects
    
    path2sensordata = [path2subj subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
    path2spatialfilter = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\spatialfilter\'];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtual_sensors\' run filesep  ];
    
    if ~exist([dir_out 'virtsens_ns\virtsens_ns_all_conditions_', run '_' freq, '.mat'], 'file')
        [vs_Run1] = adi_regroup_subjects([], path2sensordata, path2spatialfilter, dir_out, subject_list(i).name, pattern, trigger, freq, num2str(1));
        vs_Run1.label =  vs_Run1.run1.label;
        vs_Run1.fsample =  vs_Run1.run1.fsample;
    
    % sanity check Nr. 3:
   
    path_avg = [dir_out '\sanity_check\avg\'];
    if ~exist ( path_avg, 'dir')
        mkdir (path_avg)
    end
    cfg = [];
    avg = ft_timelockanalysis(cfg, vs_Run1);
    figure;
    plot(avg.time, avg.avg)
    savefig([dir_out filesep 'sanity_check\avg\avg_virtsens_allRuns_' freq '.fig'])
    
    % Safe freq spectrum:
    path_spect = [dir_out '\sanity_check\freq_spectrum\'];
    if ~exist ( path_spect, 'dir')
        mkdir (path_spect)
    end
    
    sRate = vs_Run1.fsample;
    [FourRef,Fref]=fftBasic(avg.avg, round(sRate));
    figure;
    plot(FourRef, Fref)
    savefig([path_spect 'avg_virtsens_Run1_freq_spectrum_' freq '.fig'])
    close
    
    pathAppended = [dir_out '\virtsens_ns\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end

    save ([pathAppended 'virtsens_ns_all_conditions_Run1_', freq, '.mat'], 'vs_allRuns', '-v7.3');
    
%     else
%         load([dir_out '\virtsens_ns\virtsens_ns_all_conditions_allRuns_', freq, '.mat'], 'vs_allRuns_appended')
    end
%     if ~isequal(length(vs_allRuns_appended.trial), length(vs_allRuns_appended.response_label))
%        error([subject_list(i).name ' : trial and resonse_label length do not match']) 
%     end
    virtsens_all_subj(i) = vs_Run1;
end

virtsens_all_subj


%     
%    session = virtsens_all_subj(1);
%    for k = 2:length(virtsens_all_subj)
%        session.trial = cat(2,session.trial, virtsens_all_subj(k).trial);
%        session.time = cat(2,session.time, virtsens_all_subj(k).time);
%        session.balldesign = cat(2,session.balldesign, virtsens_all_subj(k).balldesign);
%        session.triggerlabel = cat(2,session.triggerlabel, virtsens_all_subj(k).triggerlabel);
%        session.response_label = cat(2,session.response_label, virtsens_all_subj(k).response_label);
%    end
%     
%    clearvars virtsens_all_subj
%    
% for i = 1:length(session.trial)
%     temp = session.trial{1,i};
%     session.data(i,:,:) = temp;
%     clear temp
% end
%    session = rmfield(session, 'trial');
%    session.label = session.response_label;
%     
%    group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\pattern_comparison\comparison_response_labels\guggenmos_svm_results\';
%    svm_guggenmos (session, group_path, freq, condition, atlas_downsampled)
%     
    %%
    
    %% Berechnung der mittleren Abweichung zwischen Beginn des Stimulus und Beginn des IOD:
  
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
path2data = 'MEG_analysis\noisereduced\1_95Hz\01_clean\';
deviance = struct();

for i=1:length(ListSubj) 
    [deviance] = adi_deviance_trig_IOD ([fieldtripPath ListSubj(i).name filesep path2data], ListSubj(i).name, deviance, i);  
  
end

for i=1:length(deviance)
    for k=1:length(deviance(i).file)
        mean_diff_file(k) = mean(deviance(i).file(k).diff);
    end 
    mean_diff_subject(i) = mean(mean_diff_file);
end

deviance_skalar = mean(mean_diff_subject);

fs = 1017.25;
deviance_msec = 1000*deviance_skalar/fs;
