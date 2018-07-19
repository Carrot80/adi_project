
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


%% lcmv Beamformer
% calculations of spatial filter

%5 trials mitteln! ROI-Analyse
fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02b_bpfreq\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
%      Figoutpath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\source_analysis\single_subjects\'  ListSubj(i).name filesep]);
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
   
    for i = 1 : length(ListSubj)  
        path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02b_bpfreq\'];
        pathAppended = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\source_analysis\single_subjects\'  ListSubj(i).name filesep]);
        outPath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\allRuns\SVM_results\'];
        latency = [-0.4:0.01:0.97; -0.37:0.01:1];
        
        [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'bp1-45Hz', 'virtsens');
        adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'bp1-45Hz', latency, 'virtsens')
        clear virtsens_dislike_allRuns virtsens_like_allRuns
        
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'bp1-45Hz', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'bp1-45Hz', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
%  
        [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'delta',  'virtsens');
        adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'delta', latency, 'virtsens')
        clear virtsens_dislike_allRuns virtsens_like_allRuns
        
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'delta', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'delta', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
%          
        [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'theta', 'virtsens');
        adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'theta', latency, 'virtsens')
        clear virtsens_dislike_allRuns virtsens_like_allRuns
        
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'theta', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'theta', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
%         
        [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'alpha','virtsens');
        adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'alpha', latency, 'virtsens')
        clear virtsens_dislike_allRuns virtsens_like_allRuns
        
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'alpha', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'alpha', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
%         
        [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'beta', 'virtsens');
        adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'beta', latency, 'virtsens')
        clear virtsens_dislike_allRuns virtsens_like_allRuns
        
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'beta', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'beta', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
%         
        [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'low_gamma', 'virtsens');
        adi_virt_sens_SVM_allRuns (path2data, outPath, pathAppended, virtsens_like_allRuns, virtsens_dislike_allRuns, 'like', 'dislike', 'low_gamma', latency, 'virtsens')
        clear virtsens_dislike_allRuns virtsens_like_allRuns
        
%         [virtsens_ns_dislike_allRuns, virtsens_ns_like_allRuns] = adi_append_virtSensMEG(path2data, pathAppended, 'low_gamma', 'virtsens_ns');
%         adi_virt_sens_SVM_allRuns (path2data, outPath, virtsens_ns_like_allRuns, virtsens_ns_dislike_allRuns, 'like', 'dislike', 'low_gamma', latency, 'virtsens_ns')
%         clear virtsens_ns_dislike_allRuns virtsens_ns_like_allRuns
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
% noch nicht durchgeführt
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



