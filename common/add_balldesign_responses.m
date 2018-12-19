function add_balldesign_response(mainpath, path2data, filename, subject, triggercode_labels, triggercodes)

fileList = dir(fullfile([mainpath subject filesep path2data '*.mat']));   
    for k = 1:length(fileList)
        disp(['loading ' fileList(k).name ' of ' subject])
        data = load ([mainpath subject filesep path2data fileList(k).name]); 
        if isfield(data.(filename).trialinfo, 'sampleinfo')
            data.(filename).trialinfo.triggerlabel = data.(filename).trialinfo.sampleinfo(:,5)';
        end
        if ~isfield(data.(filename).trialinfo, 'triggerchannel') && exist([mainpath subject '\MEG_EEG_input\noisereduced\1_95Hz\02_Export_Bst2FT\' fileList(k).name], 'file')
            load([mainpath subject '\MEG_EEG_input\noisereduced\1_95Hz\02_Export_Bst2FT\' fileList(k).name])
            RetVal.trial(data.(filename).rejectedTrials)=[];
            for p = 1:length(data.(filename).trial)
                for j = 1:length(RetVal.trial)
                    if isequal(data.(filename).trial{p}(1,:), RetVal.trial{j}(1,:)) %&& isequal(data.(filename).trialinfo.triggerlabel(p), RetVal.trial{j}(272,find(ismember(RetVal.trial{j}(272,:), triggercodes),1)))
                       data.(filename).trialinfo.triggerchannel(p,:) = RetVal.trial{j}(272,:);
                    end
                end
            end
            
        end
        if ~isequal(size(data.(filename).trialinfo.triggerlabel,2), length(data.(filename).trial)) 
            error('unequal trial size')
        end
        data.(filename).trialinfo.trl_num_corrected_photodiode = zeros(1, length(data.(filename).trial));
        if isfield(data.(filename).trialinfo, 'triggerchannel') && ~isequal(size(data.(filename).trialinfo.triggerchannel,1), length(data.(filename).trial)) 
            warning('unequal trial size')
        elseif isfield(data.(filename).trialinfo, 'triggerchannel')  && isequal(size(data.(filename).trialinfo.triggerchannel,1), length(data.(filename).trial)) 
            for m = 1:length(data.(filename).trial)
                overlapping_trigger = find(data.(filename).trialinfo.triggerchannel(m,:) == data.(filename).trialinfo.triggerlabel(m)+512);
                if ~isempty(overlapping_trigger)
                    warning(['overlapping photodiode and stimulus trigger found in trialnumber ' num2str(m) ' in subject ' subject ' of ' fileList(k).name ' cleanMEG. Correcting...'])
                    stim_onset_sample = find(data.(filename).time{m}==0);
                    time_offset = stim_onset_sample - overlapping_trigger(1);
                    data.(filename).trialinfo.triggerchannel(m, overlapping_trigger(1))
                    data.(filename).trial{m}(end+1,:) = data.(filename).trialinfo.triggerchannel(m,:);
                    corrected_trial = nan(size(data.(filename).trial{m},1), size(data.(filename).trial{m},2));
                    corrected_trial(:,time_offset+1:size(data.(filename).trial{m},2)) = data.(filename).trial{m}(:,1:size(data.(filename).trial{m},2)-time_offset) ;  
                    data.(filename).trial{m} = corrected_trial;
                    data.(filename).trialinfo.trl_num_corrected_photodiode(m) = m;
                    data.(filename).trialinfo.(['trl_corrected_photodiode_' num2str(m)]) = corrected_trial;
                    data.(filename).trialinfo.triggerchannel(m,overlapping_trigger) = 512;
                    data.(filename).trialinfo.triggerchannel(m,time_offset+1:end) = data.(filename).trialinfo.triggerchannel(m,1:end-time_offset);
                    clearvars corrected_trial overlapping_trigger
                end
            end
            if isfield(data.(filename).trialinfo, 'triggerchannel') 
                stim_onset_sample = find(data.(filename).time{m}==0);
                indx = zeros(length(data.(filename).trial),1);
                for m = 1:length(data.(filename).trial)
                    temp = find(data.(filename).trialinfo.triggerchannel(m,1:stim_onset_sample-1) == 512);                
                    if ~isempty(temp)
                        indx(m)=temp(1);
                        time_offset =  stim_onset_sample - indx(m);
                        data.(filename).trial{m}(end+1,:) = data.(filename).trialinfo.triggerchannel(m,:);
                        corrected_trial = nan(size(data.(filename).trial{m},1), size(data.(filename).trial{m},2));
                        corrected_trial(:,time_offset+1:size(data.(filename).trial{m},2)) = data.(filename).trial{m}(:,1:size(data.(filename).trial{m},2)-time_offset) ;  
                        data.(filename).trial{m} = corrected_trial;
                        data.(filename).trialinfo.trl_num_corrected_photodiode(m) = m;
                        data.(filename).trialinfo.(['trl_corrected_photodiode_' num2str(m)]) = corrected_trial;
                        data.(filename).trialinfo.triggerchannel(m,time_offset+1:end) = data.(filename).trialinfo.triggerchannel(m,1:end-time_offset);
                        clearvars corrected_trial
                    end
                    clearvars temp 
                end
            end 
        end
        
        
        if 1 == contains(fileList(k).name, 'dislike')
                data.(filename).trialinfo.response_label(1:length(data.(filename).trial)) = {'dislike'};
                data.(filename).trialinfo.response(1:length(data.(filename).trial)) = 32;
        elseif 1 == contains(fileList(k).name, 'dontcare')
                data.(filename).trialinfo.response_label(1:length(data.(filename).trial)) = {'dontcare'};
                data.(filename).trialinfo.response(1:length(data.(filename).trial)) = 64;
        elseif strcmp(fileList(k).name(1:end-6), 'like500')
                data.(filename).trialinfo.response_label(1:length(data.(filename).trial)) = {'like'};
                data.(filename).trialinfo.response(1:length(data.(filename).trial)) = 16;
        end

        for m=1:length(data.(filename).trialinfo.triggerlabel)
            data.(filename).trialinfo.balldesign(m) = triggercode_labels(find(triggercodes == data.(filename).trialinfo.triggerlabel(m)));
        end

        save ([mainpath subject filesep path2data fileList(k).name], '-struct', 'data'); 
        clearvars data


    end
end