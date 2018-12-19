function adi_interpolate_MEG (inPath, outPath,trigger)

list = dir(fullfile([inPath, '*.mat'])); 
for k = 1:length(list)
%     if ~exist([outPath, list(k).name], 'file');
        disp(['loading ' list(k).name ' of clean MEG data'])
        load([inPath list(k).name]); 
        disp(['loading ' list(k).name ' of interpolated MEG data'])
        load([outPath list(k).name]); 
        if contains(list(k).name, 'dontcare')
            if 1==isfield(cleanMEG_interp, 'trialinfo') && 1==isfield(cleanMEG_interp.trialinfo, 'sampleinfo')
                ind_all = find(cleanMEG_interp.trialinfo.sampleinfo(:,5) == 0);  
                for p=1:length(ind_all)
                     ind_single = find(trigger.eprime == cleanMEG_interp.trialinfo.sampleinfo(ind_all(p),4));
                     cleanMEG_interp.trialinfo.sampleinfo(ind_all(p),5) = trigger.triggerchannel(ind_single);
                end
                if ~isempty(ind_all)
                    save ([outPath, list(k).name], 'cleanMEG_interp'); 
                end
            elseif 0==isfield(cleanMEG_interp, 'trialinfo')
                 cleanMEG_interp.trialinfo = cleanMEG.trialinfo;
                 save ([outPath, list(k).name], 'cleanMEG_interp'); 
            end    
        elseif ~isfield(cleanMEG_interp, 'trialinfo') 
            cleanMEG_interp.trialinfo = cleanMEG.trialinfo;
            save ([outPath, list(k).name], 'cleanMEG_interp'); 
        end
        if isfield(cleanMEG_interp.trialinfo, 'sampleinfo')
            disp(cleanMEG_interp.trialinfo.sampleinfo)
        else
            warning('cleanMEG_interp.trialinfo does not contain a sampleinfo');
        end
        clearvars cleanMEG_interp
        
        if contains(list(k).name, 'dontcare')
            if 1==isfield(cleanMEG, 'trialinfo') && 1==isfield(cleanMEG.trialinfo, 'sampleinfo')
                ind_all = find(cleanMEG.trialinfo.sampleinfo(:,5) == 0);  
                for p=1:length(ind_all)
                     ind_single = find(trigger.eprime == cleanMEG.trialinfo.sampleinfo(ind_all(p),4));
                     cleanMEG.trialinfo.sampleinfo(ind_all(p),5) = trigger.triggerchannel(ind_single);
                end
                save ([inPath, list(k).name], 'cleanMEG'); 
            elseif 0==isfield(cleanMEG, 'trialinfo')
                 cleanMEG.trialinfo = cleanMEG.trialinfo;
                 save ([inPath, list(k).name], 'cleanMEG'); 
            end
          
        elseif ~isfield(cleanMEG, 'trialinfo') 
            cleanMEG.trialinfo = cleanMEG.trialinfo;
            save ([inPath, list(k).name], 'cleanMEG_interp'); 
        end
        if isfield(cleanMEG.trialinfo, 'sampleinfo')
            disp(cleanMEG.trialinfo.sampleinfo)
        else
            warning('cleanMEG.trialinfo does not contain a sampleinfo');
        end
%         % für jedes Trial einzeln Sensoren mit 'NaN' finden:
%         for m = 1:length(cleanMEG.trial)
%             [pos(:,m)]       = isnan (cleanMEG.trial{1,m}(1:248,1));
%         end
% 
%         [badChans_allTrials] = find(sum(transpose(pos))== length(cleanMEG.trial)); 
%         
%         if ~isempty(badChans_allTrials)
%             cfgn.trials         = 'all';        % or a selection given as a 1xN vector (default = 'all')
%             cfgn.badchannel    = cleanMEG.label(badChans_allTrials'); % wichtig: semicolon, sonst funktioniert es nicht!
%             [cleanMEG_interp] = ft_channelrepair(cfgn, cleanMEG);
%         else
%             cleanMEG_interp = cleanMEG;
%         end
%         clear pos
%         for m = 1:length(cleanMEG_interp.trial)
%             [pos(:,m)]       = isnan (cleanMEG_interp.trial{1,m}(1:248,1));
%         end
% 
%         [badChans_singleTrials] = find(sum(transpose(pos)) < length(cleanMEG.trial) & sum(transpose(pos)) > 0);
%         if ~isempty(badChans_singleTrials)  
%             for m=1:length(badChans_singleTrials)
%                 indTrl = find(pos(badChans_singleTrials(m),:));
%                 cfgn.trials         = indTrl;  
%                 cfgn.badchannel    = cleanMEG.label(badChans_singleTrials(m)); % wichtig: semicolon, sonst funktioniert es nicht!
%                 cfgn.m    = cleanMEG.label(badChans_singleTrials(m)); % w
%                 [cleanMEG_interp_temp] = ft_channelrepair(cfgn, cleanMEG_interp);
%                 if length(cfgn.trials) == 1
%                     cleanMEG_interp.trial{1, cfgn.trials} = cleanMEG_interp_temp.trial;
%                 else
%                     cleanMEG_interp.trial(1, cfgn.trials) = cleanMEG_interp_temp.trial;
%                 end
%                 clear cleanMEG_interp_temp
%             end
%         end
%         clear pos
%         for m = 1:length(cleanMEG_interp.trial)
%             [pos(:,m)]       = isnan (cleanMEG_interp.trial{1,m}(1:248,1));
%         end
%         
%         if any(sum(pos))
%             disp('cleanMEG_interp still contains NaNs, using selfwritten function:')
%             [cleanMEG_interp] = interpolate_MEG(cleanMEG, cleanMEG_interp, pos, neighbours)         
%         end
%         clear pos      
%         for m = 1:length(cleanMEG_interp.trial)
%             [pos(:,m)]       = isnan (cleanMEG_interp.trial{1,m}(1:248,1));
%         end    
%         if any(sum(pos))
%             error('cleanMEG_interp still contains NaNs, check data')
%         end
%    
%         save ([outPath, list(k).name], 'cleanMEG_interp'); 
%         clearvars -except k list inPath outPath
    clear cleanMEG 

    end
end

function [cleanMEG_interp] = interpolate_MEG(cleanMEG, cleanMEG_interp, pos, neighbours)

 badChans_singleTrials = find(sum(transpose(pos)) < length(cleanMEG.trial) & sum(transpose(pos)) > 0);
 badChans_labels = cleanMEG.label(badChans_singleTrials);
 for m = 1:length(badChans_labels)
      bad_chan = char(badChans_labels(m));
      bad_chan_num = str2num(bad_chan(2:end));
      bad_chan_neighb = neighbours(bad_chan_num).neighblabel; 
      for n = 1:length(bad_chan_neighb)
          ind = find(strcmp(bad_chan_neighb{n}, cleanMEG.label));
          neighbour_index(n) = ind; % Zielvariable Nachbarn: hier stehen die indices der Nachbarn drin
          clear ind
      end
      trl_index = find(pos(badChans_singleTrials(m),:)); % Zielvariable Trialnummern
      for n = 1:length(trl_index)
          trl_interp = cleanMEG_interp.trial{1,trl_index(n)};
          activation_neighbours=trl_interp(neighbour_index, :);
          activation_neighbours(find(isnan(activation_neighbours(:,1))),:) = []; 
          interp = sum(activation_neighbours)./size(activation_neighbours,1);
          trl_interp(badChans_singleTrials(m),:) = interp;
          cleanMEG_interp.trial{1,trl_index(n)} = trl_interp;
          clear trl_interp activation_neighbours interp
      end
      clear bad_chan bad_chan_num bad_chan_neighb trl_index
 end

end
