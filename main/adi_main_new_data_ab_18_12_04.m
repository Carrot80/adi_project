%% 4.12.2018
% neue Daten


%% main settings:
clear
brainstormPath     = 'E:\adidas\Daten_Franzi\Database12-2018\BS_TrainingGroup\brainstorm_db\';
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '0.5_95Hz';
    
 %% Export Brainstorm-Files to Fieldtrip:   
    
 for i = 27:length(ListSubj)    
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
filter  = '0.5_95Hz';


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
MEG_filter     = '0.5_95Hz';

 for i = 16:20%length(ListSubj)    

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



%% SVM based on virtual sensors per Subject and run - lcmv Beamformer
% calculations of spatial filter, multiplies it with avgdata and runs Support
% Vector machine
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 3:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\0.5_95Hz\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\0.5_95Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep 'MEG\sourcespace\spatialfilter\']);
     condition = {'like', 'dislike', 'dontcare'};
   
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'beta', condition)
   
   
end


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
freq = 'beta';
condition.like = 1;
condition.dislike = 2;
condition.dontcare = 3;
load('U:\My Documents\MATLAB\atlas_source_indices.mat');
virtsens_all_subj = struct('trial', [], 'time', [], 'balldesign', [], 'triggerlabel', [], 'label', [], 'fsample', [], 'response_label', []);

for i=1:length(subject_list)
    
    vs_allRuns = struct();
    path2sensordata = [path2subj subject_list(i).name '\MEG_analysis\noisereduced\0.5_95Hz\02_interpolated\'];
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
    



