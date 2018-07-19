      
function [neighbours] = kh_neighours (data)

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