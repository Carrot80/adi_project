function adi_interpolate_MEG (inPath, outPath)

list = dir(fullfile([inPath, '*.mat'])); 
for k = 1:length(list)
    if ~exist([outPath, list(k).name], 'file')
        load([inPath list(k).name]); 
        
        [neighbours] = MEG_neighbours (cleanMEG); 
        close all
        cfgn                = [];
        cfgn.method         = 'average';   % 'average', 'spline', 'slap' or 'nan' (default = 'weighted')
        % cfgn.missingchannel = chans_cell;   % cell-array, see FT_CHANNELSELECTION for details
        % cfgn.layout        = '4D248.lay'; % weglassen, da grad sonst aus layout-file aufgebaut wird
        cfgn.neighbours     = neighbours;   % bourhood structure, see also FT_PREPARE_NEIGHBOURS
        cfgn.senstype     = 'MEG';
        [cleanMEG_interp] = ft_channelrepair(cfgn, cleanMEG);
        
        if sum(cleanMEG_interp.trial{1,1}(225,:)) == 0
            cfgn.badchannel = {'A248'};
            [cleanMEG_interp] = ft_channelrepair(cfgn, cleanMEG_interp);
        end
        
        
        
        % für jedes Trial einzeln Sensoren mit 'NaN' finden:
        for m = 1:length(cleanMEG_interp.trial)
            [pos(:,m)]       = isnan (cleanMEG_interp.trial{1,m}(1:248,1));
        end

        [badChans_allTrials] = find(sum(transpose(pos))== length(cleanMEG.trial)); 
        
        if ~isempty(badChans_allTrials)
            cfgn.trials         = 'all';        % or a selection given as a 1xN vector (default = 'all')
            cfgn.badchannel    = cleanMEG.label(badChans_allTrials'); % wichtig: semicolon, sonst funktioniert es nicht!
            [cleanMEG_interp] = ft_channelrepair(cfgn, cleanMEG_interp);
        end
        clear pos
        for m = 1:length(cleanMEG_interp.trial)
            [pos(:,m)]       = isnan (cleanMEG_interp.trial{1,m}(1:248,1));
        end

        [badChans_singleTrials] = find(sum(transpose(pos)) < length(cleanMEG_interp.trial) & sum(transpose(pos)) > 0);
        if ~isempty(badChans_singleTrials)  
            for m=1:length(badChans_singleTrials)
                indTrl = find(pos(badChans_singleTrials(m),:));
                cfgn.trials         = indTrl;  
                cfgn.badchannel    = cleanMEG.label(badChans_singleTrials(m)); % wichtig: semicolon, sonst funktioniert es nicht!
                cfgn.m    = cleanMEG.label(badChans_singleTrials(m)); % w
                [cleanMEG_interp_temp] = ft_channelrepair(cfgn, cleanMEG_interp);
                if length(cfgn.trials) == 1
                    cleanMEG_interp.trial{1, cfgn.trials} = cleanMEG_interp_temp.trial;
                else
                    cleanMEG_interp.trial(1, cfgn.trials) = cleanMEG_interp_temp.trial;
                end
                clear cleanMEG_interp_temp
            end
        end
        clear pos
        for m = 1:length(cleanMEG_interp.trial)
            [pos(:,m)]       = isnan (cleanMEG_interp.trial{1,m}(1:248,1));
        end
        
        if any(sum(pos))
            disp('cleanMEG_interp still contains NaNs, using selfwritten function:')
            [cleanMEG_interp] = interpolate_MEG(cleanMEG, cleanMEG_interp, pos, neighbours)         
        end
        clear pos      
        for m = 1:length(cleanMEG_interp.trial)
            [pos(:,m)]       = isnan (cleanMEG_interp.trial{1,m}(1:248,1));
        end    
        if any(sum(pos))
            error('cleanMEG_interp still contains NaNs, check data')
        end
        cleanMEG_interp.trialinfo = cleanMEG.trialinfo;    
       
        %% save outfile:
        if ~exist(outPath, 'dir')
            mkdir(outPath)
        end
        save ([outPath, list(k).name], 'cleanMEG_interp'); 
        clearvars -except k list inPath outPath

    end
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
