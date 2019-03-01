function [data_out] = adi_interpolate_MEG (data_in, data_source)
    data_out = data_in;
    if 1 == strcmp(data_source, 'virtsens')
       load ('E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\neighbours_distance.mat'); 
    else
         [neighbours] = MEG_neighbours (data_in(1)); 
         close all
    end
    
for p=1:length(data_in)

    cfgn                = [];
    cfgn.method         = 'average';   % 'average', 'spline', 'slap' or 'nan' (default = 'weighted')
    % cfgn.missingchannel = chans_cell;   % cell-array, see FT_CHANNELSELECTION for details
    % cfgn.layout        = '4D248.lay'; % weglassen, da grad sonst aus layout-file aufgebaut wird
    cfgn.neighbours     = neighbours;   % bourhood structure, see also FT_PREPARE_NEIGHBOURS
    cfgn.senstype     = 'MEG';

    % für jedes Trial einzeln Sensoren mit 'NaN' finden:
    for m = 1:length(data_in(p).trial)
        [pos_nan(:,m)]       = isnan (data_in(p).trial{1,m}(1:248,100));
    end
    sum(pos_nan(:))
        
    if 0==isequal(sum(pos_nan(:)), 0)
        error(['NaNs in data of subject ' num2str(p)])
    end

    for m = 1:length(data_in(p).trial)
        [sum_(:,m)]     = sum(data_in(p).trial{1,m}(1:248,:)');
    end


    try
          for oo=1:size(sum_,2)
            ind=find(sum_(:,oo)==0);
            if ~isempty(ind)
                cfgn.trials         = oo;        % or a selection given as a 1xN vector (default = 'all')
                cfgn.badchannel    = data_in(p).label(ind'); 
                [interp] = ft_channelrepair(cfgn, data_in(p));
                data_out(p).trial(oo)=interp.trial;
                clearvars interp
            end
         end    
        
    catch
        
    end
    clearvars sum_ pos_nan 
end

clearvars data_in

end
