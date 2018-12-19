function [] = add_triggercode_labels(path2subjects, subject_list, path_add_triggers, fn_add_triggers, path_trigger_sources, fn_trigger_sources, color, triggercodes, triggercode_labels)


for i = 1%[3 4 6 7 8 10 11 12 13 15]%:length(subject_list)

    switch i
        case {1, 2, 5, 9, 14}
            fileList_trigger_sources = dir(fullfile([path2subjects subject_list(i).name filesep path_trigger_sources '\alt\*.mat'])); 
            fileList_add_triggers = dir(fullfile([path2subjects subject_list(i).name filesep path_add_triggers '\alt\*.mat'])); 
        otherwise
            fileList_trigger_sources = dir(fullfile([path2subjects subject_list(i).name filesep path_trigger_sources '*.mat'])); 
            fileList_add_triggers = dir(fullfile([path2subjects subject_list(i).name filesep path_add_triggers '\*.mat'])); 
    end
              
    for k=4:length(fileList_add_triggers)
        
        disp(['loading source file ' fn_trigger_sources ' ' fileList_trigger_sources(k).name ' of subject ' subject_list(i).name])
        switch i
            case {1, 2, 5, 9, 14}
                file_triggersources = load([path2subjects subject_list(i).name filesep path_trigger_sources '\alt\' fileList_trigger_sources(k).name], fn_trigger_sources);
            otherwise
                file_triggersources = load([path2subjects subject_list(i).name filesep path_trigger_sources fileList_trigger_sources(k).name], fn_trigger_sources);
        end
        
        disp(['loading file to add triggers ' fn_add_triggers ' ' fileList_add_triggers(k).name ' of subject ' subject_list(i).name])
        switch i
            case {1, 2, 5, 9, 14}
                file_add_triggers = load([path2subjects subject_list(i).name filesep path_add_triggers '\alt\'  fileList_add_triggers(k).name], fn_add_triggers);
            otherwise
                 file_add_triggers = load([path2subjects subject_list(i).name filesep path_add_triggers filesep fileList_add_triggers(k).name], fn_add_triggers);
        end
        
        if isfield(file_add_triggers.(fn_add_triggers), 'sampleinfo') && isequal(file_add_triggers.(fn_add_triggers).sampleinfo, file_triggersources.(fn_trigger_sources).sampleinfo)
            
            file_add_triggers.(fn_add_triggers).trialinfo.triggerchannel = file_triggersources.(fn_trigger_sources).trialinfo.triggerchannel;
            file_add_triggers.(fn_add_triggers).trialinfo.responsechannel = file_triggersources.(fn_trigger_sources).trialinfo.responsechannel;
            file_add_triggers.(fn_add_triggers).trialinfo.triggerlabel = file_triggersources.(fn_trigger_sources).trialinfo.triggerlabel;
            file_add_triggers.(fn_add_triggers).trialinfo.response = file_triggersources.(fn_trigger_sources).trialinfo.response;
            file_add_triggers.(fn_add_triggers).trialinfo.response_label = file_triggersources.(fn_trigger_sources).trialinfo.response_label;
        elseif isfield(file_add_triggers.(fn_add_triggers), 'sampleinfo_orig') && isequal(file_add_triggers.(fn_add_triggers).sampleinfo_orig, file_triggersources.(fn_trigger_sources).sampleinfo_orig)
            file_add_triggers.(fn_add_triggers).trialinfo.triggerchannel = file_triggersources.(fn_trigger_sources).trialinfo.triggerchannel;
            file_add_triggers.(fn_add_triggers).trialinfo.responsechannel = file_triggersources.(fn_trigger_sources).trialinfo.responsechannel;
            file_add_triggers.(fn_add_triggers).trialinfo.triggerlabel = file_triggersources.(fn_trigger_sources).trialinfo.triggerlabel;
            file_add_triggers.(fn_add_triggers).trialinfo.response = file_triggersources.(fn_trigger_sources).trialinfo.response;
            file_add_triggers.(fn_add_triggers).trialinfo.response_label = file_triggersources.(fn_trigger_sources).trialinfo.response_label;

        else 
            warning('not adding trialinfo, as sampleinfo does not match')
            break
        end 
        cleanMEG_interp = file_add_triggers.(fn_add_triggers);
        switch i
            case {1, 2, 5, 9, 14}
                save ([path2subjects subject_list(i).name filesep path_add_triggers '\alt\' fileList_add_triggers(k).name], 'cleanMEG_interp')
            otherwise
                save ([path2subjects subject_list(i).name filesep path_add_triggers filesep fileList_add_triggers(k).name], 'cleanMEG_interp')
        end
        clearvars file_triggersources file_add_triggers cleanMEG_interp
        disp('done');
        
    end
        
end
    

end

