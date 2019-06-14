%% umfasst Konversion der Dateien aus Brainstorm-Datenbank vom Brainstorm-Format zum Fieldtrip-Format und Abpeicherung der Fieldtrip-Dateien in Pfad der Fieldtrop-Auswertung;  
% Datum der Erstellung: 3.4.2018 bis 27.09.2018
%% main settings:
    
    brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\';
    fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(fieldtripPath);
    ListSubj(1:2) = [];
    filter     = '1_95Hz';
    
 %% Export Brainstorm-Files to Fieldtrip:   
 
 for i = 2%:length(ListSubj)    
     path_export_bst2ft      = strcat(fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter,'\02_Export_Bst2Ft\');
     if ~exist(strcat(fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter, '\02_Export_Bst2Ft\dislike500_3.mat'), 'file')
        
        %run1:
        path_bst_files_run1          =  strcat(brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\c_rfhp0_notch_band\');
        adi_select_files (path_bst_files_run1, path_export_bst2ft, num2str(1))

        %run2:
        if 1 == strcmp(ListSubj(i).name, 'nl_adi_06')
            path_bst_files_run2      =  strcat(brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\c_rfhp0_notch_band_02\');
        else
            path_bst_files_run2          =  strcat(brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\c_rfhp0_notch_02_band\');
        end
        adi_select_files (path_bst_files_run2, path_export_bst2ft, num2str(2))
         
%         %run3:
        if 1 == strcmp(ListSubj(i).name, 'nl_adi_05')
            path_bst_files_run3      =  ([brainstormPath 'adi_visual_1_45Hz\adi_visuell_1_45Hz\' ListSubj(i).name '\c_rfhp0_band_03\']);
        else
            path_bst_files_run3          =  strcat(brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\c_rfhp0_notch_03_band\');
        end
        adi_select_files (path_bst_files_run3, path_export_bst2ft, num2str(3))        
%      end    
 end
    
    %% entferne Kanäle und Trials, die in den bereinigten alten Dateien gespeichert waren (mittels "Auslesen_Der_artifakte.m")
    % pathFranzi: entferne von Franzi entfernte Sensoren und Trials in
    % eigener Datei
    PathFranzi = 'W:\neurochirurgie\science\Franzi\brainstorm_db\adi_visual_Franzi\adi_visual_Franzi\';
    for i = 5%length(ListSubj)  
        path2cleanfile = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\01_clean\');
        adi_artifact_cleaningMEG(PathFranzi, fieldtripPath, path2cleanfile, ListSubj(i).name, filter) 
        
    end
    % nl_adi_11: EEG entfernen, sofern noch existent:
    path_clean = strcat(fieldtripPath, 'nl_adi_11\MEG_analysis\noisereduced\', filter, '\01_clean\');
    list_clean = dir(fullfile(strcat(path_clean, '*mat')));
    for p =1:length(list_clean)
        load([path_clean, list_clean(p).name])
        for q = 1:length(cleanMEG.trial)
            cleanMEG.trial{1,q}(272:end,:) = [];
        end
        cleanMEG.label(272:end,:) = [];
        save([path_clean, list_clean(p).name], 'cleanMEG') 
        clear cleanMEG
    end
    
    
    
        %% rejectvisual zusätzliche Säuberung einzelner Runs:
    for i = 2% : length(ListSubj)  
       
        path2cleanfile = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\01_clean\']);
        adi_rejectvisual_MEG_extra (path2cleanfile)
    end
    
    
    %% interpolate missing channels:
        
    brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\';
    fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(fieldtripPath);
    ListSubj(1:2) = [];
    filter     = '1_95Hz';
    
    for i = 2 %length(ListSubj)  
        path2cleanfile = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\' filter '\01_clean\']);
        pathInterpolated = ([fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\' filter '\02_interpolated\']);
        
        adi_interpolate_MEG (path2cleanfile, pathInterpolated)
    end
 
%  path_nl_04 = load(strcat(fieldtripPath,'nl_adi_04', '\MEG_analysis\noisereduced\', filter, '\02_interpolated\dislike500_2.mat'));  
%  load path_nl_04;
%  [neighbours] = MEG_neighbours (cleanMEG_interp); 
%  cfgn                = [];
%  cfgn.method         = 'weighted';
%  cfgn.badchannel    = {'A247'; 'A248'};
%  cfgn.trials         = 13;  
%  cfgn.neighbours     = neighbours;   % bourhood structure, see also FT_PREPARE_NEIGHBOURS
%  cfgn.senstype     = 'MEG';
%  [cleanMEG_interp2] = ft_channelrepair(cfgn, cleanMEG_interp);  
%  cleanMEG_interp.trial{1,13} = cleanMEG_interp2.trial{1,1}
%  save path_nl_04
%      
    %% bandpass-filter per Run:
    
    for i = 1 : length(ListSubj) 
        pathInterpolated = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02_interpolated\');
        pathFreq = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02b_bpfreq\');
        if ~exist (pathFreq, 'dir')
            mkdir(pathFreq)
        end
        adi_bpfilter(pathInterpolated, pathFreq, 'bp1-45Hz')
%         adi_bpfilter(pathInterpolated, pathFreq, 'theta')
%         adi_bpfilter(pathInterpolated, pathFreq, 'alpha')
%         adi_bpfilter(pathInterpolated, pathFreq, 'delta')
%         adi_bpfilter(pathInterpolated, pathFreq, 'beta')
%         adi_bpfilter(pathInterpolated, pathFreq, 'low_gamma')
%         adi_bpfilter(pathInterpolated, pathFreq, 'high_gamma')
        % noch verifizieren, dass Filtern geklappt hat
    end
    
    %% sensor level statistics per Run:
    
    latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
    for i = 1 %: length(ListSubj) 
        pathFreq = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02b_bpfreq\');
        pathInterpolated = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02_interpolated\');
        pathStatistics_perRun = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02c_statistics_perRun\');
        if ~exist(pathStatistics_perRun, 'dir')
            mkdir (pathStatistics_perRun)
        end
        
%         adi_statistics_sensor (pathInterpolated, pathStatistics_perRun, '1_95Hz', 'like', 'dislike', latency)
%         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, 'bp1-45Hz', 'like', 'dislike', latency)
        adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, 'delta', 'like', 'dislike', latency)
%         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, 'theta', 'like', 'dislike', latency)
% %         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, theta, 'like', 'dontcare', latency)
% %         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, theta, 'dislike', 'dontcare', latency)
%         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, 'alpha', 'like', 'dislike', latency)
% %         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, alpha, 'like', 'dontcare', latency)
% %         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, alpha, 'dislike', 'dontcare', latency)
%         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, 'beta', 'like', 'dislike', latency)
% %         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, beta, 'like', 'dontcare', latency)
% %         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, beta, 'dislike', 'dontcare', latency)       
%         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, 'low_gamma', 'like', 'dislike', latency)
% %         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, low_gamma, 'like', 'dontcare', latency)
% %         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, low_gamma, 'dislike', 'dontcare', latency)
%         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, 'high_gamma', 'like', 'dislike', latency)
%         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, high_gamma, 'like', 'dontcare', latency)
%         adi_freqstatistics_sensor (pathFreq, pathStatistics_perRun, high_gamma, 'dislike', 'dontcare', latency)
    end
    
    %% appenddata (runs per subject):
    filter ='bp1-45Hz';
    
    for i = 3 : 4%length(ListSubj)  
        pathInterpolated = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02_interpolated\');
        pathFreq = strcat(fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02b_bpfreq\');
        pathAppended = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\03_appended_data\');
        adi_appenddata_MEG(pathInterpolated, pathAppended, filter);
        adi_appenddata_freq(pathFreq, pathAppended, 'bp1-45Hz')
%         adi_appenddata_freq(pathFreq, pathAppended, 'delta')
%         adi_appenddata_freq(pathFreq, pathAppended, 'theta')
%         adi_appenddata_freq(pathFreq, pathAppended, 'alpha')
%         adi_appenddata_freq(pathFreq, pathAppended, 'beta')
%         adi_appenddata_freq(pathFreq, pathAppended, 'low_gamma')
%         adi_appenddata_freq(pathFreq, pathAppended, 'high_gamma')

    end
    
    
    %% avg pro Frequenzband und Proband:
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_95Hz';
 
 for i = 1 : length(ListSubj)  
     pathAppended = ([fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\03_appended_data\']);    
     list_path_appended_dislike = dir(fullfile(pathAppended, 'dislike*.mat'))
     for p = 1:length(list_path_appended_dislike)
         adi_avg_MEG(pathAppended, list_path_appended_dislike(p).name, 'dislike')
     end
     
      list_path_appended_like = dir(fullfile(pathAppended, 'like*.mat'))
      for p = 1:length(list_path_appended_like)
         adi_avg_MEG(pathAppended, list_path_appended_like(p).name, 'like')
      end 
      
      list_path_appended_dontcare = dir(fullfile(pathAppended, 'dontcare*.mat'))
      for p = 1:length(list_path_appended_dontcare)
         adi_avg_MEG(pathAppended, list_path_appended_dontcare(p).name, 'dontcare')
      end 
      
 end
    
%% statistics appended data
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
  for i = 1 : length(ListSubj) 
      pathAppended = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\03_appended_data\');
      pathStatistics = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\04_timelockstatistics_interp\');
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'bp1-45Hz', 'like', 'dislike', latency)
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, '1_95Hz', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'delta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'theta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'alpha', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'beta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'low_gamma', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'high_gamma', 'like', 'dislike', latency) 
  end

%% Test: statistics appended data mit dml.standardizer alpha = 0.8 
latency = [0:0.01:0.6; 0.03:0.01:0.63]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
  for i = 1% : length(ListSubj) 
      pathAppended = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\03_appended_data\');
      pathStatistics = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\06_timelockstatistics_dml_0.8\');
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'bp1-45Hz', 'like', 'dislike', latency)
%       adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, '1_95Hz', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'delta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'theta', 'like', 'dislike', latency) 
      adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'alpha', 'like', 'dislike', latency) 
%       adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'beta', 'like', 'dislike', latency) 
%       adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'low_gamma', 'like', 'dislike', latency) 
%       adi_crossvalidation_allRuns_perSubj (pathAppended, pathStatistics, 'high_gamma', 'like', 'dislike', latency) 
  end
  
  
  
  
  
%% group_analysis: append all subjects
pathAppended_allSubj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\01_appended_data\';

adi_append_group_MEG (fieldtripPath, pathAppended_allSubj, filter)
adi_append_group_MEG (fieldtripPath, pathAppended_allSubj, 'bp1-45Hz')
adi_append_group_MEG (fieldtripPath, pathAppended_allSubj, 'delta')
adi_append_group_MEG (fieldtripPath, pathAppended_allSubj, 'theta')
adi_append_group_MEG (fieldtripPath, pathAppended_allSubj, 'alpha')
adi_append_group_MEG (fieldtripPath, pathAppended_allSubj, 'beta')
adi_append_group_MEG (fieldtripPath, pathAppended_allSubj, 'low_gamma')
adi_append_group_MEG (fieldtripPath, pathAppended_allSubj, 'high_gamma')

%% group files are too big for cache: data need to be reduced
pathAppended_allSubj_red = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\01_appended_data_red\';
pathAppended_allSubj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\01_appended_data\';

adi_reduce_group_size_MEG(pathAppended_allSubj, pathAppended_allSubj_red)
 
%% append group: dontcare vs. like/dislike paired  

pathAppended_allSubj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\01_appended_data_paired\';
adi_append_group_paired_dontcare (fieldtripPath, pathAppended_allSubj, filter)
adi_append_group_paired_dontcare (fieldtripPath, pathAppended_allSubj, 'bp1-45Hz')
adi_append_group_paired_dontcare (fieldtripPath, pathAppended_allSubj, 'delta')
adi_append_group_paired_dontcare (fieldtripPath, pathAppended_allSubj, 'theta')
adi_append_group_paired_dontcare (fieldtripPath, pathAppended_allSubj, 'alpha')
adi_append_group_paired_dontcare (fieldtripPath, pathAppended_allSubj, 'beta')
adi_append_group_paired_dontcare (fieldtripPath, pathAppended_allSubj, 'low_gamma')
adi_append_group_paired_dontcare (fieldtripPath, pathAppended_allSubj, 'high_gamma')

%% group_analysis: sensor level statistics

pathAppended_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\01_appended_data_red\';
pathStatistics_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\SMV_crossvalidation_30ms\';
latency = [-0.5:0.01:0.97; -0.47:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
% adi_crossvalidation_group_MEG(pathAppended_group, pathStatistics_group, latency, '1_95Hz')
% adi_crossvalidation_group_MEG(pathAppended_group, pathStatistics_group, latency, 'bp1-45Hz')
% adi_crossvalidation_group_MEG(pathAppended_group, pathStatistics_group, latency, 'alpha')
% adi_crossvalidation_group_MEG(pathAppended_group, pathStatistics_group, latency, 'beta')
% adi_crossvalidation_group_MEG(pathAppended_group, pathStatistics_group, latency, 'theta')
adi_crossvalidation_group_MEG(pathAppended_group, pathStatistics_group, latency, 'delta')
adi_crossvalidation_group_MEG(pathAppended_group, pathStatistics_group, latency, 'low_gamma')
% adi_crossvalidation_group_MEG(pathAppended_group, pathStatistics_group, latency, 'high_gamma')

%% group analyis: sensor level statistics paired group

pathAppended_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\01_appended_data_paired\';
pathStatistics_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_level_statistics\02_statistics_paired\';
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erhöhen - gibt es Unterschiede?
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, latency, '1_95Hz')
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, latency, 'bp1-45Hz')
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, latency, 'delta')
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, latency, 'theta')
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, latency, 'alpha')
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, latency, 'beta')
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, latency, 'low_gamma')
adi_crossvalidation_group_paired_MEG_sensor(pathAppended_group, pathStatistics_group, latency, 'high_gamma')


%% figure: erstellt Abbildung, in der SVM Klassifikation aller Probanden pro Frequenzband dargestellt wird 
    
 subjpath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
 path2fig = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\SVM_crossvalidation\figures\';
 time = -0.5:0.01:0.98; 
 figure_comparison_singleSubjects(subjpath, 'like_vs_dislike_allRuns', 'bp1-45Hz', time, path2fig)   
 figure_comparison_singleSubjects(subjpath, 'like_vs_dislike_allRuns', 'delta', time, path2fig)    
 figure_comparison_singleSubjects(subjpath, 'like_vs_dislike_allRuns', 'theta', time, path2fig) 
 figure_comparison_singleSubjects(subjpath, 'like_vs_dislike_allRuns', 'alpha', time, path2fig) 
 figure_comparison_singleSubjects(subjpath, 'like_vs_dislike_allRuns', 'beta', time, path2fig)    


%% export warped anatomy:

brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\';
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
MEG_filter     = '1_95Hz';

 for i = 1:length(ListSubj)    

    outpath_T1warped = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
    if ~exist(outpath_T1warped, 'dir')
        mkdir(outpath_T1warped)
    end
    T1warped_Bst   =  load([brainstormPath 'adi_visuell_' MEG_filter filesep 'adi_visuell_' MEG_filter filesep ListSubj(i).name '\subjectimage_T1_warped.mat']);
    ftMri = out_fieldtrip_mri(T1warped_Bst, 'anatomy');
    %     adi_select_BstFiles_T1(path_bst_T1warped, path_T1warped_ft)
    save([outpath_T1warped 'T1warped'], 'ftMri')

 end
 
 %% read and realign mri   

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];

for i = 1:length(ListSubj) 
    mniPath = [fieldtripPath ListSubj(i).name '\MEG_EEG_input\T1_warped\'];
    outPath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
    if ~exist(outPath)
        mkdir (outPath)
    end
    adi_prepare_hdm(mniPath, outPath, ListSubj(i).name)
end


%% load mri into fieldtrip
% dieser Abschnitt wurde nicht verwendet
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];

for i = 1:length(ListSubj)   
path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
load ([path_T1warped_ft 'T1warped.mat']);

% mri = ft_determine_coordsys(ftMri, 'interactive', 'yes');
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'spm';
mri_spm    = ft_volumerealign(cfg, ftMri);

cfg            = [];
cfg.resolution = 1;
cfg.dim        = [256 256 256];
mri_spm_rs     = ft_volumereslice(cfg, mri_spm);
transform_vox2spm = mri_spm_rs.transform;
% save the transformation matrix:
save([path_T1warped_ft 'transform_vox2spm'], 'transform_vox2spm');
% save the resliced anatomy in a FreeSurfer compatible format
cd path_T1warped_ft
cfg             = [];
cfg.filename    = ListSubj(i).name;
cfg.filetype    = 'nifti';
cfg.parameter   = 'anatomy';
ft_volumewrite(cfg, mri_spm_rs);

end


%% check sensor position

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
    path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
    path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\01_clean\'];
    outPath_run1 = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\run1\'];
    outPath_run2 = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\run2\'];
    outPath_run3 = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\run3\'];
    
    adi_check_sensor_position(path2vol, path2data, outPath_run1, num2str(1))
    adi_check_sensor_position(path2vol, path2data, outPath_run2, num2str(2))
    adi_check_sensor_position(path2vol, path2data, outPath_run3, num2str(3))

end

%% adi recoding of responses: da einige Probanden bei einigen Fußbällen sowohl like, dislike, und dontcare angegeben haben, werden diese
% Ball-Antworten recodiert und in die Kategorie "dontcare" zugeordnet.
% Beispiel: nl_adi_04 beurteilte den rwv-Ball mal positiv (like), mal
% negativ(dislike) => Antworten und Stimuli des rwv-Balls werden in die
% Kategorie 'dontcare' verschoben

fieldtrippath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';

subjects = {'nl_adi_04'; 'nl_adi_05'; 'nl_adi_08'; 'nl_adi_12'; 'nl_adi_17'};
ball_combinations2recode.nl_adi_04 = {'rwv', 'gbf'};
ball_combinations2recode.nl_adi_05 = {'rwv', 'ggs', 'gbf'};
ball_combinations2recode.nl_adi_08 = {'rwv', 'rws', 'rwf'};
ball_combinations2recode.nl_adi_12 = {'ggf', 'rwf'};
ball_combinations2recode.nl_adi_17 = {'ggf', 'ggv'};

% triggercodes = cell(9,2);

triggercodepattern = {'gbv', 'rws', 'rwf', 'ggv', 'gbs', 'gbf', 'rwv', 'ggs', 'ggf'}; 
triggercodes = {'102', '104', '106', '108', '4196', '4198', '4200', '4202', '4204'};

adi_recoding_of_responses(fieldtrippath, subjects, ball_combinations2recode, triggercodes, triggercodepattern, num2str(1));
adi_recoding_of_responses(fieldtrippath, subjects, ball_combinations2recode, triggercodes, triggercodepattern, num2str(2));
adi_recoding_of_responses(fieldtrippath, subjects, ball_combinations2recode, triggercodes, triggercodepattern, num2str(3));

%% SVM based on virtual sensors per Subject and run - lcmv Beamformer
% calculations of spatial filter, multiplies it with avgdata and runs Support
% Vector machine
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 10%:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep 'MEG\sourcespace\spatialfilter\']);
     condition = {'like', 'dislike', 'dontcare'};
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'bp1-45Hz', condition)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'delta', condition)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'theta', condition)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'alpha', condition)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'beta', condition)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'low_gamma', condition)
%    
   
end

%% append virtual sensor output (runs per subject):
% Gruppenvergleich nochmal durchführen (siehe unten)
% dmlt toolbox funktion nicht richtig, kann gelöscht werden
clear    
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
   
    for i = 10%:length(ListSubj)  
        path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
        path_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);
        outPath = [path_extdisc 'MEG\sourcespace\noROIs\'];
        condition = {'like', 'dislike', 'dontcare'};
        
        [virtsens_allRuns] = adi_append_virtSensMEG(path2data, path_extdisc, 'bp1-45Hz', 'virtsens', condition);
%         adi_virt_sens_SVM_allRuns (outPath, path_extdisc, virtsens_allRuns, condition, 'bp1-45Hz', 'virtsens')
        clear virtsens_allRuns
%         
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'bp1-45Hz', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (outPath, pathAppended, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'bp1-45Hz', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
% %  
%         [virtsens_allRuns] = adi_append_virtSensMEG(path2data, path_extdisc, 'delta', 'virtsens', condition);
%         adi_virt_sens_SVM_allRuns (outPath, path_extdisc, virtsens_allRuns, condition, 'delta', 'virtsens')
%         clear virtsens_allRuns
                
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'delta', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (outPath, pathAppended, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'delta', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
%       
%         [virtsens_allRuns] = adi_append_virtSensMEG(path2data, path_extdisc, 'theta', 'virtsens', condition);
%         adi_virt_sens_SVM_allRuns (outPath, path_extdisc, virtsens_allRuns, condition, 'theta', 'virtsens')
%         clear virtsens_allRuns
%         
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'theta', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (outPath, pathAppended, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'theta', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
        
%         [virtsens_allRuns] = adi_append_virtSensMEG(path2data, path_extdisc, 'alpha', 'virtsens', condition);
%         adi_virt_sens_SVM_allRuns (outPath, path_extdisc, virtsens_allRuns, condition, 'alpha', 'virtsens')
%         clear virtsens_allRuns 
        
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'alpha', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (outPath, pathAppended, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'alpha', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
        
%         [virtsens_allRuns] = adi_append_virtSensMEG(path2data, path_extdisc, 'beta', 'virtsens', condition);
%         adi_virt_sens_SVM_allRuns (outPath, path_extdisc, virtsens_allRuns, condition, 'beta', 'virtsens')
%         clear virtsens_allRuns
%         
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'beta', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (outPath, pathAppended, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'beta', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
        
%         [virtsens_allRuns] = adi_append_virtSensMEG(path2data, path_extdisc, 'low_gamma', 'virtsens', condition);
%         adi_virt_sens_SVM_allRuns (outPath, path_extdisc, virtsens_allRuns, condition, 'low_gamma', 'virtsens')
%         clear virtsens_allRuns
        
        
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'low_gamma', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (outPath, pathAppended, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'low_gamma', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
% %                       
    end
   
    
%% figure: erstellt Abbildung, in der SVM Klassifikation aller Probanden pro Frequenzband dargestellt wird 
    % uses fill function to create nice figures
 subjpath = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
 path2fig = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\source_space\MEG\SVM_crossvalidation\figures\';
 time = -0.4:0.01:0.97; 
 figure_comparison_singleSubjects2(subjpath, 'virtsens_like_vs_dislike_allRuns', 'bp1-45Hz', time, path2fig)   
 figure_comparison_singleSubjects2(subjpath, 'virtsens_like_vs_dislike_allRuns', 'delta', time, path2fig)    
 figure_comparison_singleSubjects2(subjpath, 'virtsens_like_vs_dislike_allRuns', 'theta', time, path2fig) 
 figure_comparison_singleSubjects2(subjpath, 'virtsens_like_vs_dislike_allRuns', 'alpha', time, path2fig) 
 figure_comparison_singleSubjects2(subjpath, 'virtsens_like_vs_dislike_allRuns', 'beta', time, path2fig)    


%% group analysis: append virtual sensors of all subjects and compute SVM:
% funktioniert nicht, out of memory
path2extdisc = 'E:\Adidas\fieldtrip_Auswertung\source_analysis\single_subjects\';
pathStatistics_group = 'E:\Adidas\fieldtrip_Auswertung\source_analysis\group_analysis\MEG\source_space\SVM_crossvalidation\';


[group_data_like, group_data_dislike] = adi_append_group_MEG_sourcespace (path2extdisc, 'bp1-45Hz', 'virtsens');
adi_crossvalidation_group_MEG_sourcespace(group_data_like, group_data_dislike, pathStatistics_group, latency, 'bp1-45Hz')
% hier weitermachen:
[group_data_like, group_data_dislike] = adi_append_group_MEG_sourcespace (path2extdisc, 'delta');
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'delta')

[group_data_like, group_data_dislike] = adi_append_group_MEG_sourcespace (path2extdisc, 'theta');
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'theta')

[group_data_like, group_data_dislike] = adi_append_group_MEG_sourcespace (path2extdisc, 'alpha');
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'alpha')

[group_data_like, group_data_dislike] = adi_append_group_MEG_sourcespace (path2extdisc, 'beta');
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'beta') 

[group_data_like, group_data_dislike] = adi_append_group_MEG_sourcespace (path2extdisc, 'low_gamma');
adi_crossvalidation_group_EEG(group_data_like, group_data_dislike, pathStatistics_group, latency, 'low_gamma')   


  %% keine ROIs werden extrahiert

clear    
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
   
    for i = 1 : length(ListSubj)  
        path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
        path_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);
        outPath = [path_extdisc 'MEG\sourcespace\singleROIs\'];
        
        adi_virt_sens_SVM_noROIs (outPath, path_extdisc, 'bp1-45Hz', 'virtsens')
        adi_virt_sens_SVM_noROIs (outPath, path_extdisc, 'delta', 'virtsens')
        adi_virt_sens_SVM_noROIs (outPath, path_extdisc, 'theta', 'virtsens')
        adi_virt_sens_SVM_noROIs (outPath, path_extdisc, 'alpha', 'virtsens')
        adi_virt_sens_SVM_noROIs (outPath, path_extdisc, 'beta', 'virtsens')
        adi_virt_sens_SVM_noROIs (outPath, path_extdisc, 'low_gamma', 'virtsens')                   
    end




    %% ähnlich wie oben, neu: einzelne ROIs werden extrahiert, funktioniert noch nciht
    % hier weitermachen
clear    
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
   
    for i = 3 : length(ListSubj)  
        path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
        path_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);
        outPath = [path_extdisc 'MEG\sourcespace\singleROIs\'];
        
        adi_virt_sens_SVM_singleROIS (outPath, path_extdisc, 'bp1-45Hz', 'virtsens')
        adi_virt_sens_SVM_singleROIS (outPath, path_extdisc, 'delta', 'virtsens')
        adi_virt_sens_SVM_singleROIS (outPath, path_extdisc, 'theta', 'virtsens')
        adi_virt_sens_SVM_singleROIS (outPath, path_extdisc, 'alpha', 'virtsens')
        adi_virt_sens_SVM_singleROIS (outPath, path_extdisc, 'beta', 'virtsens')
        adi_virt_sens_SVM_singleROIS (outPath, path_extdisc, 'low_gamma', 'virtsens')                   
    end




%% SANITY CHECK: teilweise sehen SVM-Ergebnisse (ohne ROIs) merkwürdig aus, deshalb diese function
clear 
subjmainpath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(subjmainpath);
ListSubj(1:2) = [];
path2atlas = ('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\atlas_source_indices.mat');

 for i = 1 : length(ListSubj) 
     
     subjpath = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep 'MEG\sourcespace\noRois\']); 
     explore_virtsens(subjpath, path2atlas, 'bp1-45Hz')
     
 end



%% lcmv Beamformer
% calculations of spatial filter, multiplies it with avgdata 
% extracts mean and max of each ROI
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep 'MEG\sourcespace\allROIs\']);
     adi_spatial_filter_perRun_ROIs_all90 (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_all90 (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'delta', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_all90 (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'theta', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_all90 (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'alpha', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_all90 (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'beta', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_all90 (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'low_gamma', 'like', 'dislike', 'dontcare')
   
end

%% passt Datenstruktur an LibSVM bzw. Python an udn führt Guggenmos Machine Learning durch (vergleich von oben)
% schließt an adi_spatial_filter_perRun_ROIs_all90.m an
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
 for i = 3:length(ListSubj) 
        outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
        path2data = [outPath_extdisc 'MEG\sourcespace\'];
        sessions = [];
        conditions = {'like'; 'dislike'; 'dontcare'};
        %freqbands = { 'delta'; 'theta'; 'alpha'; 'beta'; 'low_gamma'};
        freqbands = { 'bp1-45Hz'; 'delta'; 'theta'; 'alpha'; 'beta'; 'low_gamma'};
        [session] = adi_append_virtsens_all90ROIs_all_freqbands(sessions, path2data, outPath_extdisc, freqbands, conditions, 'mean');
%         adi_svm_matlab_all_freqbands (session, outPath_extdisc)% => noch nicht fertig
        Guggenmos_machineLearning_all90ROIs(session, outPath_extdisc) % => kommen leider merkwürdige Ergebnisse raus
            
 end
  
 %% passt Datenstruktur an LibSVM bzw. Python an udn führt Guggenmos Machine Learning durch (vergleich von oben)
  % funktioniert super 
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
condition = {'like'; 'dislike'; 'dontcare'};
 for i = 1:length(ListSubj) 
        outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
        path2data = [outPath_extdisc 'MEG\sourcespace\'];
        Guggenmos_machineLearning_noROIs(path2data, outPath_extdisc, 'bp1-45Hz', condition)
        Guggenmos_machineLearning_noROIs(path2data, outPath_extdisc, 'alpha', condition)
        Guggenmos_machineLearning_noROIs(path2data, outPath_extdisc, 'delta', condition)
        Guggenmos_machineLearning_noROIs(path2data, outPath_extdisc, 'theta', condition)
        Guggenmos_machineLearning_noROIs(path2data, outPath_extdisc, 'low_gamma', condition)
        Guggenmos_machineLearning_noROIs(path2data, outPath_extdisc, 'beta', condition)
 end
 
 
 %% figure: erstellt Abbildung, in der SVM Klassifikation aller Probanden pro Frequenzband dargestellt wird 
    % uses fill function to create nice figures
 subjpath = 'E:\adidas\fieldtrip_Auswertung\single_subjects\';
 path2fig = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\SVM_crossvalidation\all_1457_voxels\figures\';
 algo = {'svm'; 'WeiRD'; 'GNB'};
 figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'bp1-45Hz', path2fig, algo)   
 figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_', 'delta', path2fig, algo)    
 figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'theta', path2fig, algo) 
 figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'alpha', path2fig, algo) 
 figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'beta', path2fig, algo)    
 figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'low_gamma', path2fig, algo) 

%% hier weitermachen
subjpath = 'E:\adidas\fieldtrip_Auswertung\single_subjects\';
result_path = '\MEG\sourcespace\noRois\Guggenmos_decoding_results\';
 path2fig = 'E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\SVM_crossvalidation\all_1457_voxels\bars\';
 algo = {'svm'; 'WeiRD'; 'GNB'};
 perc_acc = ([50 70 80 90 96;  70 80 90 96 101]);
 time_beg = [50, 60, 290, 470, 870]; %(ms)
 time_end = [150, 90, 320, 560, 950];
 Subject_num = {'nl_adi_04'; 'nl_adi_05'; 'nl_adi_06'; 'nl_adi_07'; 'nl_adi_08'; 'nl_adi_09'; 'nl_adi_10'; 'nl_adi_11'; 'nl_adi_12'; 'nl_adi_13'; 'nl_adi_14'; 'nl_adi_15'};
histo_group_guggenmos(subjpath, result_path, 'result_decoding_' , 'bp1-45Hz', path2fig, algo, time_beg, time_end, perc_acc, Subject_num)   
histo_group_guggenmos(subjpath, result_path, 'result_decoding_' , 'delta', path2fig, algo, time_beg, time_end, perc_acc, Subject_num)   %  figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'theta', path2fig, algo) 
histo_group_guggenmos(subjpath, result_path, 'result_decoding_' , 'theta', path2fig, algo, time_beg, time_end, perc_acc, Subject_num)   %  figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'theta', path2fig, algo) 
histo_group_guggenmos(subjpath, result_path, 'result_decoding_' , 'alpha', path2fig, algo, time_beg, time_end, perc_acc, Subject_num)   %  figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'theta', path2fig, algo) 
histo_group_guggenmos(subjpath, result_path, 'result_decoding_' , 'beta', path2fig, algo, time_beg, time_end, perc_acc, Subject_num)   %  figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'theta', path2fig, algo) 
histo_group_guggenmos(subjpath, result_path, 'result_decoding_' , 'low_gamma', path2fig, algo, time_beg, time_end, perc_acc, Subject_num)   %  figure_comparison_singleSubjects_guggenmos_svm(subjpath, 'result_decoding_' , 'theta', path2fig, algo) 


%%

 %% passt Datenstruktur an LibSVM bzw. Python an udn führt Guggenmos Machine Learning durch (vergleich von oben)
  % funktioniert super 
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
condition = {'like'; 'dislike'; 'dontcare'};
load('U:\My Documents\MATLAB\atlas_source_indices_without_cerebellum.mat');


 for i = 1:length(ListSubj) 
        outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
        path2data = [outPath_extdisc 'MEG\sourcespace\'];
        Guggenmos_machineLearning_all_voxels_without_cerebellum(path2data, outPath_extdisc, 'bp1-45Hz', condition, atlas_indices_without_cerebellum)
        Guggenmos_machineLearning_all_voxels_without_cerebellum(path2data, outPath_extdisc, 'alpha', condition, atlas_indices_without_cerebellum)
        Guggenmos_machineLearning_all_voxels_without_cerebellum(path2data, outPath_extdisc, 'delta', condition, atlas_indices_without_cerebellum)
        Guggenmos_machineLearning_all_voxels_without_cerebellum(path2data, outPath_extdisc, 'theta', condition, atlas_indices_without_cerebellum)
        Guggenmos_machineLearning_all_voxels_without_cerebellum(path2data, outPath_extdisc, 'low_gamma', condition, atlas_indices_without_cerebellum)
        Guggenmos_machineLearning_all_voxels_without_cerebellum(path2data, outPath_extdisc, 'beta', condition, atlas_indices_without_cerebellum)
 end



  %% passt Datenstruktur an LibSVM bzw. Python an udn führt Guggenmos Machine Learning durch (vergleich von oben)
  %(wie oben, für alle frequenzen gleichzeitig)
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
condition = {'like'; 'dislike'; 'dontcare'};
 for i = 1:length(ListSubj) 
        outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
        path2data = [outPath_extdisc 'MEG\sourcespace\'];
        freqbands = { 'bp1-45Hz'; 'delta'; 'theta'; 'alpha'; 'beta'; 'low_gamma'};
        Guggenmos_svm_1457voxels_all_freqs_simultaniously(path2data, outPath_extdisc, freqbands, condition)
        
 end
 
 
   %% passt Datenstruktur an LibSVM bzw. Python an udn führt Guggenmos Machine Learning durch (vergleich von oben)
  %(wie oben, für alle frequenzen gleichzeitig)
  % tall array Versuch: funktoniert noch nicht
clear
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
condition = {'like'; 'dislike'; 'dontcare'};
 for i = 1:length(ListSubj) 
%         outPath_extdisc = (['D:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
%         path2data = [outPath_extdisc 'MEG\sourcespace\'];
        outPath_extdisc = 'D:\Kirsten\virtsens_adi_06\';
        path2data = 'D:\Kirsten\virtsens_adi_06\virtsens\MEG\sourcespace\';
        freqbands = { 'bp1-45Hz'; 'delta'; 'theta'; 'alpha'; 'beta'; 'low_gamma'};
        Guggenmos_svm_1457voxels_all_freqs_simultaniously_tall_array(path2data, outPath_extdisc, freqbands, condition)
        
 end
 %% wie oben, nur für matlab svm
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
condition = {'like'; 'dislike'; 'dontcare'};
 for i = 1:length(ListSubj) 
%         outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
%         path2data = [outPath_extdisc 'MEG\sourcespace\'];
        outPath_extdisc = 'D:\Kirsten\virtsens_adi_06\';
        path2data = 'D:\Kirsten\virtsens_adi_06\virtsens\MEG\sourcespace\';
        freqbands = { 'bp1-45Hz'; 'delta'; 'theta'; 'alpha'; 'beta'; 'low_gamma'};
        matlab_svm_1457voxels_all_freqs_simultaniously_all_array(path2data, outPath_extdisc, freqbands, condition)
        
 end
 
 
 
 %% noch schreibne
 

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
condition = {'like'; 'dislike'; 'dontcare'};
 for i = 3:length(ListSubj) 
        outPath_extdisc = (['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
        path2data = [outPath_extdisc 'MEG\sourcespace\'];
        adi_SVM_matlab(path2data, outPath_extdisc, 'bp1-45Hz', condition)
        adi_SVM_matlab(path2data, outPath_extdisc, 'alpha', condition)
        adi_SVM_matlab(path2data, outPath_extdisc, 'delta', condition)
        adi_SVM_matlab(path2data, outPath_extdisc, 'theta', condition)
        adi_SVM_matlab(path2data, outPath_extdisc, 'low_gamma', condition)
        adi_SVM_matlab(path2data, outPath_extdisc, 'beta', condition)
 end

 
 
 
 
 %% evtl. löschen
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
% bp = {'delta', 'theta', 'alpha', 'beta', 'low_gamma'};
bp = {'bp1-45Hz', 'theta', 'alpha', 'beta', 'low_gamma'};
cond = { 'like', 'dislike', 'dontcare'}; 
meanmax = 'mean';

 for i = 3 : length(ListSubj) 
        outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
        path2data = [outPath_extdisc 'MEG\sourcespace\'];
        Guggenmos_machineLearning_allfreqbands(path2data, bp, cond, meanmax)
        
 end

 
 
 
 
 
 
 %% source space for ROIs: computes virtual sensors for occipital, parietal and frontal ROIs

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
     outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep 'MEG\sourcespace\Clusters\']);
     adi_spatial_filter_perRun_ROIs_allclusters (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_allclusters (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'delta', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_allclusters (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'theta', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_allclusters (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'alpha', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_allclusters (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'beta', 'like', 'dislike', 'dontcare')
     adi_spatial_filter_perRun_ROIs_allclusters (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'low_gamma', 'like', 'dislike', 'dontcare')
  
   
end




%% passt Datenstruktur an LibSVM bzw. Python an

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
% latency = [-0.4:0.001:0.995; -0.395:0.001:1]; 
 for i = 9 : length(ListSubj) 
%         outPath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\allRuns\SVM_results\'];
        outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);  
        path2data = [outPath_extdisc 'MEG\sourcespace\'];
        Guggenmos_machineLearning2(path2data, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', 'dontcare', 'mean')
%         [virtsens_allRuns] = SVM_Guggenmos(path2data, outPath,
%         outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', 'dontcare',
%         'virtsens') % alte datei: evtl. löschen
        
 end
 

 %% aus virtsens aller 1457 Voxels werden ROIS extrahiert und Mittelwert oder Maximum pro ROI berechnet
 % alle Frequenzbänder werden in eine Matrix gesteckt und Guggenmos Machine Learning Skript durchgeführt
 % kommen merkwürdige Ergebnisse raus, deshalb unten fitcsvm
 
clear 
fieldtripPath      = 'E:\adidas\fieldtrip_Auswertung\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
load('U:\My Documents\MATLAB\atlas_ind.mat');
fieldnames_atlas = fields(ind_atlas);

for m = fieldnames_atlas(91:end)
    ind_atlas = rmfield(ind_atlas, m);
end

virtsens_all_subjects = struct('data', [], 'labels', [], 'time', []);
bp = {'bp1-45Hz', 'delta', 'theta', 'alpha', 'beta', 'low_gamma'};
condition = {'like'; 'dislike'; 'dontcare'};
 for i = 1:length(ListSubj) 
        outPath = ('E:\Adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\SVM_crossvalisation\all_90ROIs_all_freqs\');  
        path2data = ['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name '\MEG\sourcespace\'];
        virtsens_all_subjects = adi_extract_meanmax_virtsens_ROIs(virtsens_all_subjects, ind_atlas, path2data, outPath, bp, condition, 'max', i);
end           
 
  session = virtsens_all_subjects(1);
 for i = 2:length(virtsens_all_subjects)
    session.data = cat(1, session.data, virtsens_all_subjects(i).data);
    session.labels = cat (2, session.labels, virtsens_all_subjects(i).labels);
 end    
clear virtsens_all_subjects
Guggenmos_machineLearning_all90ROIs(session, outPath, 'all_freqs', 'max')

% matlab classification:
data = session.data;
ind_dontcare=find(session.labels ==3);
data(ind_dontcare,:,:)=[];
labels(ind_dontcare)=[];
classification_accuracy=nan(1,size(session.data,3));
for k = 1:size(session.data,3)
    
    X=data(:,:,k);
    SVMModel = fitcsvm(X,labels,'Standardize',true,'KernelFunction','linear',...
    'KernelScale','auto');
    CVSVMModel = crossval(SVMModel);
    classLoss = kfoldLoss(CVSVMModel, 'mode', 'average');
    classification_accuracy(k) = 1-classLoss;
    clear X SVMModel classLoss
end

save(['classification_accuracy_' 'max_6_freqs'], 'classification_accuracy')
figure
hold on
plot(session.time{1,1}, 100*classification_accuracy) 
legend('like vs dislike')
ylim([30 70])
xlabel('time [s]')
ylabel('predictive accuracy [%]')
 
 %% group anaylsis: mean of 90 ROIs: Kommen merkwürdige Ergebnisse raus
path2subjects      = 'E:\adidas\fieldtrip_Auswertung\single_subjects\';
outPath = ('E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\SVM_crossvalidation\all_90ROIs_all_freqs\');  
ListSubj = dir(path2subjects);
ListSubj(1:2) = [];
bp = {'bp1-45Hz', 'delta', 'theta', 'alpha', 'beta', 'low_gamma'};
cond = { 'like', 'dislike', 'dontcare'}; 
virtsens_all_subjects =  struct('data', [], 'labels', [], 'time', [], 'tissuelabel', [], 'ROIs', []);
 for i = 1 : length(ListSubj) 
        inPath = [path2subjects ListSubj(i).name '\MEG\sourcespace\allROIs\'];
        [virtsens_all_subjects] = adi_append_group_MEG_meanmax_ROIs(inPath, virtsens_all_subjects, bp, cond, 'mean', i);
 end

 if ~exist(outPath, 'dir')
     mkdir(outPath)
 end
 save ([outPath 'virtsens_mean_all_subjects'], 'virtsens_all_subjects', '-v7.3')
 session = virtsens_all_subjects(1);
 
for i = 2:length(virtsens_all_subjects)
    session.data = cat(1, session.data, virtsens_all_subjects(i).data);
    session.labels = cat (2, session.labels, virtsens_all_subjects(i).labels);
end    

clear virtsens_all_subjects
Guggenmos_machineLearning_all90ROIs(session, outPath, 'all_freqs', 'mean')
 



 %% group anaylsis: all clusters: 12.9.2018

clear
outPath = ('E:\adidas\fieldtrip_Auswertung\groups_analysis\source_space\MEG\SVM_crossvalidation\all_clusters_single_freqs\');  
fieldtripPath      = 'E:\adidas\fieldtrip_Auswertung\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
condition = {'like'; 'dislike'; 'dontcare'};
load('U:\My Documents\MATLAB\atlas_source_indices_all_clusters.mat');
virtsens_all_subjects =  struct('data', [], 'labels', [], 'time', []);
 for i = 1:length(ListSubj) 
        path2data = ['E:\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name  '\MEG\sourcespace\'];
        [virtsens_all_subjects] = adi_append_group_MEG_all_clusters(virtsens_all_subjects, path2data, 'bp1-45Hz', condition, atlas, i);
%       
 end
  
 session = virtsens_all_subjects(1);
 for i = 2:length(virtsens_all_subjects)
    session.data = cat(1, session.data, virtsens_all_subjects(i).data);
    session.labels = cat (2, session.labels, virtsens_all_subjects(i).labels);
end    
% Guggenmos_machineLearning_clusters(session, outPath, 'bp1-45Hz', condition, atlas)

% matlab classification:
data = session.data;
labels = session.labels;
ind_dontcare=find(session.labels ==3);
data(ind_dontcare,:,:)=[];
labels(ind_dontcare)=[];
classification_accuracy=nan(1,size(session.data,3));

for k = 1:size(data,3)
    
    X=data(:,:,k);
    SVMModel = fitcsvm(X,labels,'Standardize',true,'KernelFunction','linear',...
    'KernelScale','auto', 'OptimizeHyperparameters', 'none',  'ShrinkagePeriod', 0, 'verbose', 0);
    CVSVMModel = crossval(SVMModel);
    classLoss = kfoldLoss(CVSVMModel, 'mode', 'average');
    classification_accuracy(k) = 1-classLoss;
    clear X SVMModel classLoss
end

save(['classification_accuracy_' 'max_6_freqs'], 'classification_accuracy')
figure
hold on
plot(session.time{1,1}, classification_accuracy*100) 
legend('like vs dislike')
ylim([30 70])
xlabel('time [s]')
ylabel('predictive accuracy [%]')

% clear virtsens_all_subjects
 



