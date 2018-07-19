%% umfasst Konversion der Dateien aus Brainstorm-Datenbank vom Brainstorm-Format zum Fieldtrip-Format und Abpeicherung der Fieldtrip-Dateien in Pfad der Fieldtrop-Auswertung;  
% Datum der Erstellung: 3.4.2018
%% main settings:

    fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(fieldtripPath);
    ListSubj(1:2) = [];
    filter     = '1_45Hz';
    
 %% Export Brainstorm-Files to Fieldtrip:   
 brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\adi_visual_1_45Hz\adi_visuell_1_45Hz\';
 for i = 2%:length(ListSubj)    
     path_export_bst2ft  = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\EEG\', filter,'\02_Export_Bst2Ft\']);
     if ~exist(path_export_bst2ft, 'dir')
         mkdir (path_export_bst2ft)
     end
     if ~exist([path_export_bst2ft, 'dislike500_1.mat'], 'file')
        % Run1:
        files_dislike = ([brainstormPath, ListSubj(i).name, '\c_rfhp0_band\']);
        adi_select_files_EEG (files_dislike, path_export_bst2ft, num2str(1))
     end
     
        %run2:
    if ~exist([path_export_bst2ft, 'dislike500_2.mat'], 'file')  
       path_bst_files_run2   =  ([brainstormPath, ListSubj(i).name, '\c_rfhp0_band_02\']);
       adi_select_files_EEG (path_bst_files_run2, path_export_bst2ft, num2str(2))
    end  
%         %run3:

        path_bst_files_run3  =  ([brainstormPath, ListSubj(i).name, '\c_rfhp0_band_03\']);
        adi_select_files_EEG (path_bst_files_run3, path_export_bst2ft, num2str(3))        
 
 end
 %% enferne aus Speicherplatzgründen die MEG-Kanäle 
 for i = 2%:length(ListSubj)     
     path_export_bst2ft  = [fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\EEG\', filter,'\02_Export_Bst2Ft\'];
     path_EEG = [fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\EEG\', filter, filesep];
     adi_delete_MEG(path_export_bst2ft, path_EEG) 
 end
 
 %% entferne Kanäle und Trials, die Franzi in ihrer Auswertung entfernt hat

    PathFranzi = 'L:\Arbeit\Adidas\Daten_Franzi\Database\BS_TrainingGroup\brainstorm_db\';
    for i = 2:length(ListSubj)  
        path_EEG = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\EEG\', filter, filesep, '02_Export_Bst2Ft']);
        path2cleanfile = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\1_45Hz\01_clean\']);
        adi_artifact_cleaningEEG(PathFranzi, path_EEG, path2cleanfile, ListSubj(i).name)
    end
  
    % evtl. vereinfachen:
for m = 1:length(EEG.trial)
    [ind] = find(EEG.ChannelFlag_Bst{1,m} == -1);
            if ~isempty(ind)
               EEG.trial{1, m}(ind, :) = NaN;
            end
            clear ind
end
    
    %% interpolate missing channels:
    
    for i = 1 : length(ListSubj)  
        path2cleanfile = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\', filter, '\01_clean\']);
        pathInterpolated = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\', filter, '\02_interpolated\']);
        
        adi_interpolate_EEG (path2cleanfile, pathInterpolated)
    end
    
% funktioniert nicht richtig, besser nochmal vereinfachen:

 [neighbours] = EEG_neighbours (cleanEEG); 
        close all
        cfgn                = [];
        cfgn.method         = 'weighted';   % 'average', 'spline', 'slap' or 'nan' (default = 'weighted')
        % cfgn.missingchannel = chans_cell;   % cell-array, see FT_CHANNELSELECTION for details
        % cfgn.layout        = '4D248.lay'; % weglassen, da grad sonst aus layout-file aufgebaut wird
        cfgn.neighbours     = neighbours;   % bourhood structure, see also FT_PREPARE_NEIGHBOURS
        cfgn.senstype     = 'EEG';

        cleanEEG_interp = cleanEEG;
        
for p = 1:length(cleanEEG.trial) 
    cfgn.trials  = p; 
    [ind] = find(cleanEEG.ChannelFlag_Bst{1,p} == -1);
    label = cleanEEG.label(ind)
    cfgn.badchannel    = label;
    [cleanEEG_interp_temp] = ft_channelrepair(cfgn, cleanEEG); 
    cleanEEG_interp.trial{1, p} = cleanEEG_interp_temp.trial{1,1};
    cleanEEG_interp.time{1,p} = cleanEEG_interp_temp.time{1,1};
    clear cleanEEG_interp_temp
end


    
    %% sensor level statistics per Run:
    
    latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
    for i = 1 : length(ListSubj) 
        pathInterpolated = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\', filter, '\02_interpolated\']);
        pathStatistics_perRun = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\', filter, '\02c_statistics_perRun\']);
        if ~exist(pathStatistics_perRun, 'dir')
            mkdir (pathStatistics_perRun)
        end
        adi_freqstatistics_sensor_EEG_perRun (pathInterpolated, pathStatistics_perRun, 'bp1-45Hz', 'like', 'dislike', latency)
        adi_freqstatistics_sensor_EEG_perRun (pathInterpolated, pathStatistics_perRun, 'delta', 'like', 'dislike', latency)
        adi_freqstatistics_sensor_EEG_perRun (pathInterpolated, pathStatistics_perRun, 'theta', 'like', 'dislike', latency)
        adi_freqstatistics_sensor_EEG_perRun (pathInterpolated, pathStatistics_perRun, 'alpha', 'like', 'dislike', latency)
        adi_freqstatistics_sensor_EEG_perRun (pathInterpolated, pathStatistics_perRun, 'beta', 'like', 'dislike', latency)       
        adi_freqstatistics_sensor_EEG_perRun (pathInterpolated, pathStatistics_perRun, 'low_gamma', 'like', 'dislike', latency)
    end
    
    
    %% statistics all Runs (appended data) per subject
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
  for i = 1 : length(ListSubj) 
      pathInterpolated = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\', filter, '\02_interpolated\']);
      pathStatistics = ([fieldtripPath, ListSubj(i).name, '\EEG_analysis\', filter, '\04_timelockstatistics_interp\']);
      [data_dislike, data_like] = adi_appenddata_EEG(pathInterpolated, filter);
      adi_crossvalidation_allRuns_perSubj_EEG (data_like, data_dislike, pathStatistics, 'bp1-45Hz', 'like', 'dislike', latency)
      adi_crossvalidation_allRuns_perSubj_EEG (data_like, data_dislike, pathStatistics, 'delta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj_EEG (data_like, data_dislike, pathStatistics, 'theta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj_EEG (data_like, data_dislike, pathStatistics, 'alpha', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj_EEG (data_like, data_dislike, pathStatistics, 'beta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj_EEG (data_like, data_dislike, pathStatistics, 'low_gamma', 'like', 'dislike', latency) 
  end
    
    
    %% group_analysis: sensor level statistics

pathStatistics_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_level_statistics\02_statistics_EEG_enet_0.3\';
latency = [0:0.01:0.98; 0.02:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
[group_data_like, group_data_dislike] = adi_append_group_EEG (fieldtripPath)
% adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, '1_45Hz')
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'delta')
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'theta')
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'alpha')
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'beta') 
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'low_gamma')   
    
    
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


%% figure: erstellt Abbildung, in der SVM Klassifikation aller Probanden pro Frequenzband dargestellt wird 
    
 subjpath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
 path2fig = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\EEG\SVM_crossvalidation\figures\';
 time = -0.5:0.01:0.98; 
 adi_figure_comparison_singleSubjectsEEG(subjpath, 'like_vs_dislike_allRuns', 'bp1-45Hz', time, path2fig)   
 adi_figure_comparison_singleSubjectsEEG(subjpath, 'like_vs_dislike_allRuns', 'delta', time, path2fig)    
 adi_figure_comparison_singleSubjectsEEG(subjpath, 'like_vs_dislike_allRuns', 'theta', time, path2fig) 
 adi_figure_comparison_singleSubjectsEEG(subjpath, 'like_vs_dislike_allRuns', 'alpha', time, path2fig) 
 adi_figure_comparison_singleSubjectsEEG(subjpath, 'like_vs_dislike_allRuns', 'beta', time, path2fig)  



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
 
 %% prepare head model:
% read and realign mri   

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];

for i = 13:length(ListSubj) 
    mniPath = [fieldtripPath ListSubj(i).name '\MEG_EEG_input\T1_warped\'];
    outPath = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\05_source_space\hdm\'];
    path2data = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\01_clean\'];
    
    if ~exist(outPath)
        mkdir (outPath)
    end
    adi_prepare_hdm_EEG(mniPath, path2data, outPath, ListSubj(i).name)
end

 
 
 %% wird nicht genutzt, kann gelöscht werden:
 
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_45Hz';
 
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
%% check sensor position

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
    path2vol = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\05_source_space\hdm\'];
    path2data = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\01_clean\'];
    outPath_run1 = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\05_source_space\run1\'];
    outPath_run2 = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\05_source_space\run2\'];
    outPath_run3 = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\05_source_space\run3\'];
    
    if exist([path2data 'dislike500_1.mat'], 'file')
        adi_check_sensor_position_EEG(path2vol, path2data, outPath_run1, num2str(1))
    end
    if exist([path2data 'dislike500_2.mat'], 'file')
        adi_check_sensor_position_EEG(path2vol, path2data, outPath_run2, num2str(2))
    end
    if exist([path2data 'dislike500_3.mat'], 'file')
        adi_check_sensor_position_EEG(path2vol, path2data, outPath_run3, num2str(3))   
    end
    

end


%% lcmv beamformer per run
% calculations of spatial filter

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\05_source_space\hdm\'];
     path2data = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
     outPath = [fieldtripPath ListSubj(i).name '\EEG_analysis\1_45Hz\05_source_space\'];
     outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\source_analysis\single_subjects\'  ListSubj(i).name filesep 'EEG\']);
     
     adi_compute_spatialfilter_EEG (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike')
     adi_compute_spatialfilter_EEG (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'delta', 'like', 'dislike')
     adi_compute_spatialfilter_EEG (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'theta', 'like', 'dislike')
     adi_compute_spatialfilter_EEG (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'alpha', 'like', 'dislike')
     adi_compute_spatialfilter_EEG (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'beta', 'like', 'dislike')
     adi_compute_spatialfilter_EEG (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'low_gamma', 'like', 'dislike')
   
end



















