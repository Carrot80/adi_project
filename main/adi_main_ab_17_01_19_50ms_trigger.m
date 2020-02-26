%% 50ms trigger, ab 2.9.2019
% neue Daten


%% main settings:
clear
brainstormPath  = 'W:\neurochirurgie\science\Franzi\Database\BS_TrainingGroup\brainstorm_db\';
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_95Hz';
    
 %% Export Brainstorm-Files to Fieldtrip:   

clear
brainstormPath  = 'W:\neurochirurgie\science\Franzi\Database\BS_TrainingGroup\brainstorm_db\';
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_95Hz';
   
 for i = 29:length(ListSubj)    
     path_export_bst2ft   = ([fieldtripPath ListSubj(i).name '\MEG_EEG_input\noisereduced\' filter '\02_Export_Bst2Ft\']);
     if ~exist(path_export_bst2ft, 'dir')
         mkdir(path_export_bst2ft)
     end  
    adi_select_files_newData (brainstormPath, ListSubj(i).name, path_export_bst2ft)
 end
 
 
 
%% füge trigger und balldesign hinzu und entferne bad channels:

clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
brainstormPath  = 'W:\neurochirurgie\science\Franzi\Database\BS_TrainingGroup\brainstorm_db\';

ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter  = '1_95Hz';

trigger.balldesign_short = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
trigger.balldesign = {'yellow_blue_volley', 'red_white_space', 'red_white_soccerball', 'grey_green_volley', 'yellow_blue_space', 'yellow_blue_soccerball', 'red_white_volley', 'grey_green_space', 'grey_green_soccerball'}; 
trigger.triggercodes = [4206, 4208, 4210, 4212, 110, 112, 114, 116, 118];
trigger.eprime        = [111, 113, 115, 117, 110, 112, 114, 116, 118];

for ii = 1:length(ListSubj)  
    path2cleanfile = ([fieldtripPath, ListSubj(ii).name, '\MEG_analysis\noisereduced\', filter, '\01_clean_50ms\']);
    path2ft_data = ([fieldtripPath, ListSubj(ii).name, '\MEG_EEG_input\noisereduced\' filter '\02_Export_Bst2Ft\']);
    brainstormPath  = ['W:\neurochirurgie\science\Franzi\Database\BS_TrainingGroup\brainstorm_db\' ListSubj(ii).name filesep];
    adi_artifact_cleaning_newData_50ms(path2ft_data, brainstormPath, path2cleanfile, ListSubj(ii).name, filter, trigger)    
end

 %% adi_04 =>> channels umsortieren
clear
inpath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\nl_adi_04\MEG_analysis\noisereduced\1_95Hz\01_clean_50ms\alte_Kanalsortierung\';
path_new_labels = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\nl_adi_05\MEG_analysis\noisereduced\1_95Hz\01_clean\Neu_Like500_1.mat';
outpath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\nl_adi_04\MEG_analysis\noisereduced\1_95Hz\01_clean_50ms\';

kh_sort_channels(inpath, path_new_labels, outpath)


%% rejectvisual zusätzliche Säuberung einzelner Runs
% nl_adi_19, nl_adi21, nl_adi20        
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];

    for ii = 17%16%18
        
        path2cleanfile = ([fieldtripPath, ListSubj(ii).name, '\MEG_analysis\noisereduced\1_95Hz\01_clean_50ms\']);
        adi_rejectvisual_MEG_extra_50ms (path2cleanfile, ListSubj(ii).name)
    end
       
%% interpolate missing channels:
clear
brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\';
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_95Hz';

for ii = 1:length(ListSubj)  
    path2cleanfile = ([fieldtripPath ListSubj(ii).name '\MEG_analysis\noisereduced\' filter '\01_clean_50ms\']);
    pathInterpolated = ([fieldtripPath ListSubj(ii).name '\MEG_analysis\noisereduced\' filter '\02_interpolated_50ms\']);

    adi_interpolate_MEG (path2cleanfile, pathInterpolated)
end

%% grand avg per subject:
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
path2save = 'MEG_analysis\noisereduced\1_95Hz\grandavg_50ms\';
path2inputfile = 'MEG_analysis\noisereduced\1_95Hz\02_interpolated_50ms\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
grandavg = []; 
[grandavg, trl_count] = adi_grandavg_subject(subject_list, grandavg, [], path2save, path2inputfile);

%% grand avg
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
folderpath = 'MEG_analysis\noisereduced\1_95Hz\grandavg_50ms\';
adi_grandavg_group(subject_list, folderpath);

%% überprüfe, ob 50ms trigger fälschlicherweise zugeordet wurden: nicht genutzt

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
%% svm für sensor daten Ball gbf und gbv

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

% balldesign = 'gbf';
% subject_list([1 2 16 18 24 25],:) = [];
% delete_runs={'nl_adi_06', '2'; 'nl_adi_17', '1'; 'nl_adi_24', '2'; 'nl_adi_25', '1'}; 

balldesign = 'gbv';
subject_list([12 19 25],:) = [];
delete_runs={'nl_adi_05', '1'; 'nl_adi_05', '2'; 'nl_adi_07', '1'; 'nl_adi_09', '3'; 'nl_adi_10', '1';  'nl_adi_13', '1'; 'nl_adi_17', '1'; 'nl_adi_25', '1'; 'nl_adi_33', '1'; 'nl_adi_34', '1'}; 


load('U:\My Documents\MATLAB\atlas_source_indices.mat');
sensordata_all_subj = struct('trial', [], 'time', [], 'response_label', [], 'balldesign_short', [],'label', [], 'fsample', [] ,'grad', [], 'cfg', []);

for i=1:length(subject_list)
    path2sensordata = ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\' subject_list(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\' ];
    [sensordata_all_subj] = adi_single_ball_sensordata(sensordata_all_subj, path2sensordata, subject_list(i).name, pattern, trigger, freq, balldesign, delete_runs, i);
end

% nachträgliche Interpolation von Kanälen, da es nicht überall geklappt
% hat:

[sensordata_all_subj] = adi_additionalMEG_interpolation(sensordata_all_subj);

% z-transformation:
[sensordata_all_subj_zscore] = adi_ztrans_sensorspace(sensordata_all_subj);

% clearvars sensordata_all_subj

% statistic
time=[0.5 0.8];
adi_stats_sensorspace(sensordata_all_subj_zscore, time) % hier noch vervollstänidgen
% 

session = sensordata_all_subj_zscore(1);
   for k = 2:length(sensordata_all_subj_zscore)
       session.trial = cat(2,session.trial, sensordata_all_subj_zscore(k).trial);
       session.time = cat(2,session.time, sensordata_all_subj_zscore(k).time);
       session.response_label = cat(2,session.response_label, sensordata_all_subj_zscore(k).response_label);
       session.balldesign_short = cat(2,session.balldesign_short, sensordata_all_subj_zscore(k).balldesign_short);
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
% 
[session] = regroup_session(session);
adi_sanity_check_sensorspace(session)
adi_stats_sensorspace(session) % hier noch vervollstänidgen
% 
group_path = ['E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\svm\' balldesign filesep];
SVM_Guggenmos_adapted (session, group_path, freq)
% % SVM_Guggenmos_virtsens (session, group_path, freq)



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

% balldesign = 'gbf';
% subject_list([1 2 16 18 24 25],:) = [];
% delete_runs={'nl_adi_06', 2; 'nl_adi_17', 1; 'nl_adi_24', 2; 'nl_adi_25', 1};
% 
% balldesign = 'gbv';
% subject_list([12 19 25],:) = [];
% delete_runs={'nl_adi_05', [1 2]; 'nl_adi_07', 1; 'nl_adi_09', 3; 'nl_adi_10', 1;  'nl_adi_13', 1; 'nl_adi_17', 1; 'nl_adi_25', 1; 'nl_adi_33', 1; 'nl_adi_34', 1}; 

% balldesign = 'rwf';
% subject_list([2 16 17],:) = [];
% delete_runs={'nl_adi_8', [1 2]; 'nl_adi_10', 1;  'nl_adi_11', 3; 'nl_adi_12', [1 2]; 'nl_adi_15' [1 2]; 'nl_adi_18', [1 2]; 'nl_adi_21', [1 2]; 'nl_adi_22', 1; 'nl_adi_25', 1; 'nl_adi_26', 1; 'nl_adi_28', 1}; 
% 
balldesign = 'ggf';
delete_subjects = {'nl_adi_04'; 'nl_adi_07'; 'nl_adi_12'; 'nl_adi_18'; 'nl_adi_20'};
for i=1:length(subject_list)
    ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
end
subject_list(find(ind),:) = [];
clear ind

delete_runs={'nl_adi_17', 1; 'nl_adi_19', 1;  'nl_adi_21', [2 3]; 'nl_adi_23', 3; 'nl_adi_28' [1 2]; 'nl_adi_29', 1; 'nl_adi_34', [2 3]}; % adi28: run2 sollte auch gelöscht werden, existiert aber nicht in den Daten
% 
% balldesign = 'ggs';
% delete_subjects = {'nl_adi_20'; 'nl_adi_29'; 'nl_adi_34'};
% for i=1:length(subject_list)
%     ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
% end
% subject_list(find(ind),:) = [];
% clear ind
% 
% delete_runs={'nl_adi_04', 3; 'nl_adi_05', [2];  'nl_adi_10', [1]; 'nl_adi_15', 2; 'nl_adi_17' 1; 'nl_adi_22', [2 3]; 'nl_adi_23', [1 2]; 'nl_adi_28', [1]; 'nl_adi_33', 1} ; 


load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'response', [], 'response_label', [], 'balldesign', [], 'subject', []);


for i=1:length(subject_list)
    
    path2virtsens = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' ];
    dir_out = ['E:\adidas\fieldtrip_Auswertung\single_subjects\' subject_list(i).name '\MEG\sourcespace\virtsens\' balldesign filesep];
    [virtsens_all_subj] = adi_vs_single_ball(virtsens_all_subj, path2virtsens, dir_out, subject_list(i).name, pattern, trigger, freq, balldesign, delete_runs, i);
end

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

group_path = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\ggs_ball\stats\';
time = [0 1]; % in s
adi_stats_sourcespace(session, time, group_path, 'ggs') 




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
    