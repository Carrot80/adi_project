function adi_interpolate_EEG (inPath, outPath)

list = dir(fullfile([inPath, '*.mat'])); 
for k = 1:length(list)
    if ~exist([outPath, list(k).name], 'file');
%         load([inPath list(k).name]); 
        [neighbours] = EEG_neighbours (cleanEEG); 
        close all
        cfgn                = [];
        cfgn.method         = 'weighted';   % 'average', 'spline', 'slap' or 'nan' (default = 'weighted')
        % cfgn.missingchannel = chans_cell;   % cell-array, see FT_CHANNELSELECTION for details
        % cfgn.layout        = '4D248.lay'; % weglassen, da grad sonst aus layout-file aufgebaut wird
        cfgn.neighbours     = neighbours;   % bourhood structure, see also FT_PREPARE_NEIGHBOURS
        cfgn.senstype     = 'EEG';

        % für jedes Trial einzeln Sensoren mit 'NaN' finden:
        for m = 1:length(cleanEEG.trial);
            [pos(:,m)]       = isnan (cleanEEG.trial{1,m}(1:64,1));
        end
        if any(pos) 
            [chan_nr_Idcs, trl_Idcs] = find(pos~=0); % how many elements of each row are non-zero
            badChannels_num = unique(chan_nr_Idcs);       % temp vector of vals
            rowIdcs = sort(chan_nr_Idcs);          % sorted input aligns with temp (lowest to highest)
            t = zeros(size(badChannels_num)); % vector for freqs
            % frequency for each value
            for i = 1:length(badChannels_num)
               t(i) = sum(rowIdcs == badChannels_num(i));
            end
            ind_length_trl = find(t == length(cleanEEG.trial));
            chans_cell_all_trl = cleanEEG.label(badChannels_num(ind_length_trl));
            % which sensor contains NaN in all trials
            if ~isempty (chans_cell_all_trl)
                cfgn.trials         = 'all';        % or a selection given as a 1xN vector (default = 'all')
                cfgn.badchannel    = chans_cell_all_trl; % wichtig: semicolon, sonst funktioniert es nicht!
                [cleanEEG_interp] = ft_channelrepair(cfgn, cleanEEG);  
                clear chans_cell_all_trl;
            else
                cleanEEG_interp = cleanEEG;
            end
            % find bad channels in single trials which are not present in
            % each trial
            ind = find(t < length(cleanEEG.trial));
            [chan_nr_Idcs, trl_Idcs] = find(pos~=0); %
            for p = 1:length(badChannels_num)
                 ind_badchan_all_trl = find(chan_nr_Idcs == badChannels_num(p));
                 chan_nr_Idcs(ind_badchan_all_trl) = [];
                 trl_Idcs(ind_badchan_all_trl) = [];
            end
   
            if ~isempty (chan_nr_Idcs)
                for i = 1:length(chan_nr_Idcs) 
                    cfgn.trials  = trl_Idcs(i); 
                    cfgn.badchannel    = cleanEEG_interp.label(chan_nr_Idcs(i));
                    [cleanEEG_interp_temp] = ft_channelrepair(cfgn, cleanEEG_interp);
                    cleanEEG_interp.trial{1, cfgn.trials} = cleanEEG_interp_temp.trial{1,1};
                    clear cleanEEG_interp_temp
                end
            end
        else
             cleanEEG_interp = cleanEEG;
        end
        save ([outPath, list(k).name], 'cleanEEG_interp'); 
        clearvars -except k list inPath outPath
    end
end
end
      
function [neighbours] = EEG_neighbours (data)

% ---> % using gradiometers specified in the data, not template:

        cfg_neigh = [];
        cfg_neigh.method        = 'template';               %'distance', 'triangulation' or 'template'
%         cfg_neigh.template      = 'bti248_neighb.mat';      % name of the template file, e.g. CTF275_neighb.mat
%         cfg_neigh.layout        = '4D248.lay';              % filename of the layout, see FT_PREPARE_LAYOUT
        cfg_neigh.channel       = {'all'};                  % channels for which neighbours should be found
        cfg_neigh.feedback      = 'yes';                    % 'yes' or 'no' (default = 'no')
        neighbours = ft_prepare_neighbours(cfg_neigh, data);
        
%         cfg =[];
%         cfg.neighbours = neighbours;
%         cfg.enableedit = 'yes';
%         cfg.senstype = 'MEG';
%         neighbours2= ft_neighbourplot (cfg, data_bpf);

%         sens = ft_read_sens('c,rfhp0.1Hz,n')
%         chan_str = [];
%         chan_str{1,1} = chans(1,1)
%         chan_str{1,2} = chans(1,2)
end
