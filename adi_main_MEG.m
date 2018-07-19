%% umfasst Konversion der Dateien aus Brainstorm-Datenbank vom Brainstorm-Format zum Fieldtrip-Format und Abpeicherung der Fieldtrip-Dateien in Pfad der Fieldtrop-Auswertung;  
% Datum der Erstellung: 3.4.2018
%% main settings:
    
    brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\';
    fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(fieldtripPath);
    ListSubj(1:2) = [];
    filter     = '1_95Hz';
    
 %% Export Brainstorm-Files to Fieldtrip:   
 
 for i = 2%:length(ListSubj)    
     path_export_bst2ft      = strcat(fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter,'\02_Export_Bst2Ft\');
%      if ~exist(strcat(fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter, '\02_Export_Bst2Ft\dislike500_3.mat'), 'file')
%         
%         %run1:
%         path_bst_files_run1          =  strcat(brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\c_rfhp0_notch_band\');
%         adi_select_files (path_bst_files_run1, path_export_bst2ft, num2str(1))
% 
%         %run2:
%         if 1 == strcmp(ListSubj(i).name, 'nl_adi_06')
%             path_bst_files_run2      =  strcat(brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\c_rfhp0_notch_band_02\');
%         else
%             path_bst_files_run2          =  strcat(brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\c_rfhp0_notch_02_band\');
%         end
%         adi_select_files (path_bst_files_run2, path_export_bst2ft, num2str(2))
         
%         %run3:
        if 1 == strcmp(ListSubj(i).name, 'nl_adi_05')
            path_bst_files_run3      =  ([brainstormPath 'adi_visual_1_45Hz\adi_visuell_1_45Hz\' ListSubj(i).name '\c_rfhp0_band_03\']);
        else
            path_bst_files_run3          =  strcat(brainstormPath, 'adi_visuell_', filter, filesep, 'adi_visuell_', filter, filesep, ListSubj(i).name, '\c_rfhp0_notch_03_band\');
        end
        adi_select_files (path_bst_files_run3, path_export_bst2ft, num2str(3))        
%      end    
 end
    
    %% entferne Kan�le und Trials, die in den bereinigten alten Dateien gespeichert waren (mittels "Auslesen_Der_artifakte.m")
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
    
    
    
        %% rejectvisual zus�tzliche S�uberung einzelner Runs:
    for i = 8 : length(ListSubj)  
       
        path2cleanfile = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\01_clean\');
        adi_rejectvisual_MEG (path2cleanfile)
    end
    
    
    %% interpolate missing channels:
        
    brainstormPath     = 'W:\neurochirurgie\science\Kirsten\adidas\brainstorm_db\';
    fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(fieldtripPath);
    ListSubj(1:2) = [];
    filter     = '1_95Hz';
    
    for i = 14 %length(ListSubj)  
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
    
    latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erh�hen - gibt es Unterschiede?
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
    
    for i = 1 : length(ListSubj)  
        pathInterpolated = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02_interpolated\');
        pathFreq = strcat(fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\02b_bpfreq\');
        pathAppended = strcat (fieldtripPath, ListSubj(i).name, '\MEG_analysis\noisereduced\', filter, '\03_appended_data\');
%         adi_appenddata(pathInterpolated, pathAppended, filter);
        adi_appenddata_freq(pathFreq, pathAppended, 'bp1-45Hz')
        adi_appenddata_freq(pathFreq, pathAppended, 'delta')
        adi_appenddata_freq(pathFreq, pathAppended, 'theta')
        adi_appenddata_freq(pathFreq, pathAppended, 'alpha')
        adi_appenddata_freq(pathFreq, pathAppended, 'beta')
        adi_appenddata_freq(pathFreq, pathAppended, 'low_gamma')
        adi_appenddata_freq(pathFreq, pathAppended, 'high_gamma')

    end
    
    
    %% avg pro Frequenzband und Proband:
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
filter     = '1_95Hz';
 
 for i = 3 : length(ListSubj)  
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
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erh�hen - gibt es Unterschiede?
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
latency = [0:0.01:0.6; 0.03:0.01:0.63]; % Latentzeit probeweise auf 30ms erh�hen - gibt es Unterschiede?
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
latency = [-0.5:0.01:0.97; -0.47:0.01:1]; % Latentzeit probeweise auf 30ms erh�hen - gibt es Unterschiede?
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
latency = [-0.5:0.01:0.98; -0.48:0.01:1]; % Latentzeit probeweise auf 30ms erh�hen - gibt es Unterschiede?
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

%%

%% SVM based on virtual sensors per Subject and run - lcmv Beamformer
% calculations of spatial filter, multiplies it with avgdata and runs Support
% Vector machine

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = [1 2 5 9 14] %:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);
     latency = [-0.4:0.001:0.995; -0.395:0.001:1]; 
     
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'delta', 'like', 'dislike', latency)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'theta', 'like', 'dislike', latency)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'alpha', 'like', 'dislike', latency)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'beta', 'like', 'dislike', latency)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'low_gamma', 'like', 'dislike', latency)
     adi_spatial_filter_perRun (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', latency)
   
end


    %% append virtual sensor output (runs per subject):
    
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
   
    for i = [1 2 5 9 14]% 1 : length(ListSubj)  
        path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
        pathAppended = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);
        outPath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\allRuns\SVM_results\'];
        latency = [-0.4:0.01:0.97; -0.37:0.01:1];
        
%         [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'bp1-45Hz', 'virtsens');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'bp1-45Hz', latency, 'virtsens')
%         clear virtsens_dislike_allRuns virtsens_like_allRuns
%         
        [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'bp1-45Hz', 'virtsens_ns');
        adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'bp1-45Hz', latency, 'virtsens_ns')
        clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
 
%         [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'delta',  'virtsens');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'delta', latency, 'virtsens')
%         clear virtsens_dislike_allRuns virtsens_like_allRuns
%         
        [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'delta', 'virtsens_ns');
        adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'delta', latency, 'virtsens_ns')
        clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
%          
%         [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'theta', 'virtsens');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'theta', latency, 'virtsens')
%         clear virtsens_dislike_allRuns virtsens_like_allRuns
%         
        [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'theta', 'virtsens_ns');
        adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'theta', latency, 'virtsens_ns')
        clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
        
%         [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'alpha','virtsens');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'alpha', latency, 'virtsens')
%         clear virtsens_dislike_allRuns virtsens_like_allRuns
%         
        [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'alpha', 'virtsens_ns');
        adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'alpha', latency, 'virtsens_ns')
        clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
        
%         [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'beta', 'virtsens');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'beta', latency, 'virtsens')
%         clear virtsens_dislike_allRuns virtsens_like_allRuns
%         
        [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'beta', 'virtsens_ns');
        adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'beta', latency, 'virtsens_ns')
        clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
        
%         [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'low_gamma', 'virtsens');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'low_gamma', latency, 'virtsens')
%         clear virtsens_dislike_allRuns virtsens_like_allRuns
        
        [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'low_gamma', 'virtsens_ns');
        adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'low_gamma', latency, 'virtsens_ns')
        clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
%                       
    end
   
    
%% Beamformer ROI-Analysis 
    
    
    
    
    
    
    
    
    
    
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


%% calculates Beamformer for allRuns per subject and run SVM afterwards 
% noch nicht durchgef�hrt
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
    
   path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
   path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
   outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\source_analysis\single_subjects\'  ListSubj(i).name '\sourcespace\runs_appended\']);
   path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\'];
   outPath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\allRuns\source_reconstruction\SVM_results\'];
   latency = [-0.4:0.01:0.97; -0.37:0.01:1];
   
   if ~exist(outPath, 'dir')
      mkdir (outPath)
   end
    
   [virtsens_like, virtsens_dislike] = adi_source_reconstruction_allRuns(path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', latency, 'virtsens');
   adi_virt_sens_SVM_allRuns (outPath, path2data, outPath_extdisc, virtsens_like, virtsens_dislike, 'like', 'dislike', latency, 'bp1-45Hz', 'virtsens')
   clear virtsens_like virtsens_dislike
   
   [virtsens_like, virtsens_dislike] = adi_source_reconstruction_allRuns(path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'delta', 'like', 'dislike', latency, 'virtsens');
   adi_virt_sens_SVM_allRuns (outPath, path2data, outPath_extdisc, virtsens_like, virtsens_dislike, 'like', 'dislike', latency, 'delta', 'virtsens')
   clear virtsens_like virtsens_dislike
   
   [virtsens_like, virtsens_dislike] = adi_source_reconstruction_allRuns(path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'theta', 'like', 'dislike', latency, 'virtsens');
   adi_virt_sens_SVM_allRuns (outPath, path2data, outPath_extdisc, virtsens_like, virtsens_dislike, 'like', 'dislike', latency, 'theta', 'virtsens')
   clear virtsens_like virtsens_dislike
   
   [virtsens_like, virtsens_dislike] = adi_source_reconstruction_allRuns(path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'alpha', 'like', 'dislike', latency, 'virtsens');
   adi_virt_sens_SVM_allRuns (outPath, path2data, outPath_extdisc, virtsens_like, virtsens_dislike, 'like', 'dislike', latency, 'alpha', 'virtsens')
   clear virtsens_like virtsens_dislike
   
   [virtsens_like, virtsens_dislike] = adi_source_reconstruction_allRuns(path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'beta', 'like', 'dislike', latency, 'virtsens');
   adi_virt_sens_SVM_allRuns (outPath, path2data, outPath_extdisc, virtsens_like, virtsens_dislike, 'like', 'dislike', latency, 'beta', 'virtsens')
   clear virtsens_like virtsens_dislike
   
   [virtsens_like, virtsens_dislike] = adi_source_reconstruction_allRuns(path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'low_gamma', 'like', 'dislike', latency, 'virtsens');
   adi_virt_sens_SVM_allRuns (outPath, path2data, outPath_extdisc, virtsens_like, virtsens_dislike, 'like', 'dislike', latency, 'low_gamma', 'virtsens')
   clear virtsens_like virtsens_dislike
    
   [virtsens_like, virtsens_dislike] = adi_source_reconstruction_allRuns(path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'high_gamma', 'like', 'dislike', latency, 'virtsens');
   adi_virt_sens_SVM_allRuns (outPath, path2data, outPath_extdisc, virtsens_like, virtsens_dislike, 'like', 'dislike', latency, 'high_gamma', 'virtsens')
   clear virtsens_like virtsens_dislike
    
end




%% group analysis: append virtual sensors of all subjects and compute SVM:
% funktioniert nicht, out of memory
path2extdisc = 'L:\Arbeit\Adidas\fieldtrip_Auswertung\source_analysis\single_subjects\';
pathStatistics_group = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\source_space\MEG\SVM_crossvalidation\';
latency = [-0.4:0.01:0.97; -0.37:0.01:1]; 

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




%% adi recoding of responses: da einige Probanden bei einigen Fu�b�llen sowohl like, dislike, und dontcare angegeben haben, werden diese
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

% adi_recoding_of_responses(fieldtrippath, subjects, ball_combinations2recode, triggercodes, triggercodepattern, num2str(1));
% adi_recoding_of_responses(fieldtrippath, subjects, ball_combinations2recode, triggercodes, triggercodepattern, num2str(2));
adi_recoding_of_responses(fieldtrippath, subjects, ball_combinations2recode, triggercodes, triggercodepattern, num2str(3));

%%










%% SVM based on virtual sensors per Subject and run - lcmv Beamformer
% calculations of spatial filter, multiplies it with avgdata and runs Support
% Vector machine
% hier morgen weitermachen: 5 trials mitteln! ROI-Analyse
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep]);
     latency = [-0.4:0.001:0.995; -0.395:0.001:1]; 
     
     adi_spatial_filter_perRun_ROI (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'delta', 'like', 'dislike', latency)
     adi_spatial_filter_perRun_ROI (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'theta', 'like', 'dislike', latency)
     adi_spatial_filter_perRun_ROI (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'alpha', 'like', 'dislike', latency)
     adi_spatial_filter_perRun_ROI (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'beta', 'like', 'dislike', latency)
     adi_spatial_filter_perRun_ROI (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'low_gamma', 'like', 'dislike', latency)
     adi_spatial_filter_perRun_ROI (path2vol, path2data, path_T1warped_ft, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', latency)
   
end

