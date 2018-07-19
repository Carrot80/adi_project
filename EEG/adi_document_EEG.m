function main ()

%% main settings:
    
    submainpath     = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
    ListSubj = dir(submainpath);
    ListSubj(1:2) = [];
    filter     = '1_45Hz';
    
 %% Export Brainstorm-Files to Fieldtrip:   
 
%     for i = 1:length(ListSubj)    
%         path_export_bst2ft      = strcat(submainpath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter,'\02_Export_Bst2Ft\');
%         if ~exist(strcat(path_export_bst2ft, 'dislike500_1.mat'), 'file')
%         
%         %run1:
%         path_bst_files_run1          =  strcat(submainpath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter, '\01_Brainstorm_Files\c_rfhp0_band\');
%         adi_select_files (path_bst_files_run1, path_export_bst2ft, num2str(1))
% 
%         %run2:
%         path_bst_files_run2          =  strcat(submainpath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter, '\01_Brainstorm_Files\c_rfhp0_band_02\');
%         adi_select_files (path_bst_files_run2, path_export_bst2ft, num2str(2))
% 
%         %run3:
%         path_bst_files_run3          =  strcat(submainpath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter, '\01_Brainstorm_Files\c_rfhp0_band_03\');
%         adi_select_files (path_bst_files_run3, path_export_bst2ft, num2str(3))       
%         end
%     end
     
    %% rejectvisual
    for i = 3: length(ListSubj)  
        path_export_bst2ft = strcat(submainpath, ListSubj(i).name, '\MEG_EEG_input\noisereduced\', filter,'\02_Export_Bst2Ft\');
        path2cleanfile = strcat (submainpath, ListSubj(i).name, '\EEG_analysis\', filter, '\01_clean\');
        adi_rejectvisual_EEG (path_export_bst2ft, path2cleanfile)
    end

end


