function [] = adi_check_sensor_position(path2vol, path2data, outPath, run)
    
%% align the electrodes - dies muss für jeden Run gemacht werden:
if ~exist([outPath, 'sensorposition_check.fig'], 'file')   
    load ([path2vol 'vol_EEG.mat'], 'vol');
    vol = ft_convert_units(vol, 'mm');
    load ([path2data, 'dislike500_', run, '.mat'], 'cleanEEG' );
    % electrodes
    elec = cleanEEG.elec;
    elec = ft_convert_units(elec, 'mm');   
    figure;  
    ft_plot_mesh(vol.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]); 
    hold on;
    ft_plot_sens(elec);    

    if ~exist(outPath, 'dir')
        mkdir(outPath)
    end
    savefig([outPath, 'sensorposition_check.fig']);   
    fig = ([outPath 'electrodealignment']);
    print('-dpng', fig); 
    close all
end

% cfg           = [];
% cfg.method    = 'interactive';
% cfg.elec      = elec;
% cfg.headshape = vol.bnd(1);
% elec_aligned  = ft_electroderealign(cfg);


  
    

end
