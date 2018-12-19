%% umfasst Konversion der Dateien aus Brainstorm-Datenbank vom Brainstorm-Format zum Fieldtrip-Format und Abpeicherung der Fieldtrip-Dateien in Pfad der Fieldtrop-Auswertung;  
% Datum der Erstellung: 16.5.2018
%% main settings:

    fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(fieldtripPath);
    ListSubj(1:2) = [];
    filter     = '1_45Hz';
    
    
    %% sensor level statistics per Run:
    
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
for i = 1 : length(ListSubj) 
        pathStatistics_perRun_combined_MEEG = ([fieldtripPath, ListSubj(i).name, '\combined_MEG_EEG_analysis\02c_statistics_perRun\']);
        
        if ~exist(pathStatistics_perRun_combined_MEEG, 'dir')
            mkdir (pathStatistics_perRun_combined_MEEG)
        end
        pathInterpolated_EEG = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\1_45Hz\02_interpolated\']);
        pathInterpolated_MEG = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\']);
        adi_freqstatistics_sensorMEG_EEG_perRun_kuerzer (pathInterpolated_EEG, pathInterpolated_MEG, pathStatistics_perRun_combined_MEEG, 'bp1-45Hz', 'like', 'dislike', latency)
        adi_freqstatistics_sensorMEG_EEG_perRun_kuerzer (pathInterpolated_EEG, pathInterpolated_MEG, pathStatistics_perRun_combined_MEEG, 'delta', 'like', 'dislike', latency)
        adi_freqstatistics_sensorMEG_EEG_perRun_kuerzer (pathInterpolated_EEG, pathInterpolated_MEG, pathStatistics_perRun_combined_MEEG, 'theta', 'like', 'dislike', latency)
        adi_freqstatistics_sensorMEG_EEG_perRun_kuerzer (pathInterpolated_EEG, pathInterpolated_MEG, pathStatistics_perRun_combined_MEEG, 'alpha', 'like', 'dislike', latency)
        adi_freqstatistics_sensorMEG_EEG_perRun_kuerzer (pathInterpolated_EEG, pathInterpolated_MEG, pathStatistics_perRun_combined_MEEG, 'beta', 'like', 'dislike', latency)       
        adi_freqstatistics_sensorMEG_EEG_perRun_kuerzer (pathInterpolated_EEG, pathInterpolated_MEG, pathStatistics_perRun_combined_MEEG, 'low_gamma', 'like', 'dislike', latency)
end
   
    %% statistics all Runs (appended data) per subject
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
  for i = 1 : length(ListSubj) 
      pathStats = ([fieldtripPath, ListSubj(i).name, '\combined_MEG_EEG_analysis\02c_statistics_allRuns\']);
      pathInterpolated_EEG = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\1_45Hz\02_interpolated\']);
      pathInterpolated_MEG = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\']);
      
      [data_combined_MEG_EEG_dislike, data_combined_MEG_EEG_like] = adi_appenddata_MEEG(pathInterpolated_EEG, 'EEG', pathInterpolated_MEG, 'MEG', ListSubj(i).name)
      adi_crossvalidation_allRuns_perSubj_MEG_EEG (data_combined_MEG_EEG_dislike, data_combined_MEG_EEG_like, pathStatistics, 'bp1-45Hz', 'like', 'dislike', latency)
      adi_crossvalidation_allRuns_perSubj_MEG_EEG (data_combined_MEG_EEG_dislike, data_combined_MEG_EEG_like, pathStatistics, 'delta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj_MEG_EEG (data_combined_MEG_EEG_dislike, data_combined_MEG_EEG_like, pathStatistics, 'theta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj_MEG_EEG (data_combined_MEG_EEG_dislike, data_combined_MEG_EEG_like, pathStatistics, 'alpha', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj_MEG_EEG (data_combined_MEG_EEG_dislike, data_combined_MEG_EEG_like, pathStatistics, 'beta', 'like', 'dislike', latency)
      adi_crossvalidation_allRuns_perSubj_MEG_EEG (data_combined_MEG_EEG_dislike, data_combined_MEG_EEG_like, pathStatistics, 'low_gamma', 'like', 'dislike', latency) 

  end
    
    
    %% group_analysis: sensor level statistics
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_45Hz';

pathStatistics_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\combined_MEEG\SVM_crossvalidation\';
latency = [-0.5:0.01:0.97; -0.47:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
group_data_like = [];
group_data_dislike = [];
for i = 1:length(ListSubj)
    inPathEEG = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\1_45Hz\02_interpolated\']);
    inPathMEG = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\']);
    [group_data_like, group_data_dislike] = adi_append_group_MEEG (inPathEEG, 'EEG', inPathMEG, 'MEG', ListSubj(i).name, group_data_like, group_data_dislike, i)
end

size_like_all_subj = size(group_data_like);
cfg = [];
group_like_appended = ft_appenddata(cfg, group_data_like(:,1), group_data_like(:,2));

for k = 3:size_like_all_subj(2)
    group_like_appended = ft_appenddata(cfg, group_like_appended, group_data_like(:,k)) ;   
end
clear group_data_like

size_dislike_all_subj = size(group_data_dislike);
cfg = [];
group_dislike_appended = ft_appenddata(cfg, group_data_dislike(:,1), group_data_dislike(:,2));

for k = 3:size_dislike_all_subj(2)
    group_dislike_appended = ft_appenddata(cfg, group_dislike_appended, group_data_dislike(:,k)) ;   
end
clear group_data_dislike

adi_crossvalidation_group_MEG_EEG(group_like_appended, group_dislike_appended, pathStatistics_group, latency, '1_45Hz')
adi_crossvalidation_group_MEG_EEG(group_like_appended, group_dislike_appended, pathStatistics_group, latency, 'delta')
adi_crossvalidation_group_MEG_EEG(group_like_appended, group_dislike_appended, pathStatistics_group, latency, 'theta')
adi_crossvalidation_group_MEG_EEG(group_like_appended, group_dislike_appended, pathStatistics_group, latency, 'alpha')
adi_crossvalidation_group_MEG_EEG(group_like_appended, group_dislike_appended, pathStatistics_group, latency, 'beta') 
adi_crossvalidation_group_MEG_EEG(group_like_appended, group_dislike_appended, pathStatistics_group, latency, 'low_gamma')   
    
    
%% group analyis: sensor level statistics paired group

pathStatistics_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_level_statistics\02_statistics_paired_dontcare_EEG\';
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
% hier evtl. weiternachen: dontcare-trials müssen noch bereinigt werden 
[group_data_like, group_data_dislike, group_data_dontcare] = adi_append_group_paired_dontcare_EEG (fieldtripPath)

adi_crossvalidation_group_paired_EEG_sensor(group_data_like, group_data_dislike, pathStatistics_group, latency, '1_45Hz')
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'delta')
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'theta')
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'alpha')
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'beta') 
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'low_gamma')   



%% group analyis: sensor level statistics paired group

pathAppended_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\01_appended_data_paired\';
pathStatistics_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_level_statistics\02_statistics_paired\';
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, '1_95Hz')

%%

%% figure: erstellt Abbildung, in der SVM Klassifikation aller Probanden pro Frequenzband dargestellt wird 
    
 subjpath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
 path2fig = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\combined_MEEG\SVM_crossvalidation\figures\';
 time = -0.5:0.01:0.98; 
 adi_figure_comparison_singleSubjectsMEEG(subjpath, 'like_vs_dislike_allRuns', 'bp1-45Hz', time, path2fig)   
 adi_figure_comparison_singleSubjectsMEEG(subjpath, 'like_vs_dislike_allRuns', 'delta', time, path2fig)    
 adi_figure_comparison_singleSubjectsMEEG(subjpath, 'like_vs_dislike_allRuns', 'theta', time, path2fig) 
 adi_figure_comparison_singleSubjectsMEEG(subjpath, 'like_vs_dislike_allRuns', 'alpha', time, path2fig) 
 adi_figure_comparison_singleSubjectsMEEG(subjpath, 'like_vs_dislike_allRuns', 'beta', time, path2fig)  





%% export warped anatomy and segment mri:

 for i = 3:length(ListSubj)   
     
path_T1warped_ft    = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
if ~exist(path_T1warped_ft, 'dir')
    mkdir(path_T1warped_ft)
end
T1warped_Bst   =  load([brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\subjectimage_T1_warped.mat']);
ftMri = out_fieldtrip_mri(T1warped_Bst, 'anatomy');
%     adi_select_BstFiles_T1(path_bst_T1warped, path_T1warped_ft)
save([path_T1warped_ft, 'T1warped'], 'ftMri')

% segment mri:
cfg           = [];
cfg.output = {'gray','white','csf','skull','scalp'};
segmentedmri  = ft_volumesegment(cfg, ftMri);
segmentedmri.anatomy = ftMri.anatomy;

save([path_T1warped_ft, 'segmentedmri.mat'], 'segmentedmri')
clear segmentedmri ftMri path_T1warped_ft T1warped_Bst    

 end
 
 %% head model:

 for i = 1:length(ListSubj) 
    outPath = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\source_space\']);
    if ~exist (outPath, 'dir')
        mkdir(outPath)
    end
    fileName = [outPath 'vol.mat'];
    if ~exist(fileName, 'file')
        cfg = [];
        cfg.method = 'bemcp';
        vol = ft_prepare_headmodel(cfg, segmentedmri);
        disp(vol)
        save(fileName, 'vol');
    else
        load(fileName);
    end

    % source model
    fileName = [outpath 'sourcemodel.mat'];
    if ~exist(filename, 'file')
        sourcemodel = conn_prepare_sourcemodel(segmentedmri, vol, 1, 0, 1);
        save(filename, 'sourcemodel');
    else
        load(filename);
    end
 end
%% source level statistics