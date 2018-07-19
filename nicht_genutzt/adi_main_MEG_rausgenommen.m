%% SVM based on virtual sensors per Subject and run - lcmv Beamformer
% calculations of spatial filter, multiplies it with avgdata and runs Support
% Vector machine

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02b_bpfreq\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
     outPath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\source_analysis\single_subjects\'  ListSubj(i).name filesep]);
     latency = [-0.4:0.01:0.97; -0.37:0.01:1];
     
     adi_virt_sens_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.5])
     adi_virt_sens_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'delta', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.5])
     adi_virt_sens_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'theta', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.5])
     adi_virt_sens_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'alpha', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.5])
     adi_virt_sens_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'beta', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.5])
     adi_virt_sens_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'low_gamma', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.5])
   
end


%% SVM based on source reconstruction prestim/poststim-Ratio per Subject and run - Beamformer
% nicht verwendet
% funktioniert so nicht, da keine trials vorhanden
% calculations spatial filter, multiplies it with avgdata and runs Support
% Vector machine

fieldtripPath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(fieldtripPath);
ListSubj(1:2) = [];
for i = 1:length(ListSubj) 
     path2vol = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\vol\'];
     path2data = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\02b_bpfreq\'];
     path_T1warped_ft = ([fieldtripPath, ListSubj(i).name, '\MEG_EEG_input\T1_warped\']);
     outPath = [fieldtripPath ListSubj(i).name '\MEG_analysis\noisereduced\1_95Hz\05_source_space\'];
     outPath_extdisc = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\source_analysis\single_subjects\'  ListSubj(i).name filesep]);
     latency = [-0.4:0.01:0.97; -0.37:0.01:1];
     
     adi_sourceRatio_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'bp1-45Hz', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.50])
     adi_sourceRatio_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'delta', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.50])
     adi_sourceRatio_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'theta', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.50])
     adi_sourceRatio_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'alpha', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.50])
     adi_sourceRatio_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'beta', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.50])
     adi_sourceRatio_SVM_perRun (path2vol, path2data, path_T1warped_ft, outPath, outPath_extdisc, 'low_gamma', 'like', 'dislike', latency, [-0.5 -0.02], [0.02 0.50])
   
end