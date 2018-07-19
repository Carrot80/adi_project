function [] = adi_check_sensor_position(path2vol, path2data, outPath, run)
    
    load ([path2vol 'vol.mat']);
    vol = ft_convert_units(vol, 'cm');
    
    load ([path2data, 'dislike500_', run, '.mat']);
    grad = ft_convert_units(cleanMEG.grad, 'cm');
    ft_plot_sens(grad)
    ft_plot_vol(vol)
    hs = ft_read_headshape([outPath, 'hs_file']);
    hs = ft_convert_units(hs, 'cm');
    ft_plot_headshape(hs)
    savefig([outPath, 'sensorposition_check.fig']);   
    close all
end
% 
% ft_plot_vol(vol)
% ft_plot_sens(data_bpfreq.grad)
% hs = ft_read_headshape([outPath, 'hs_file']);
% hs = ft_convert_units(hs, 'cm');
% hs = rmfield(hs, 'fid'); % remove the fiducials -> these are stored in MRI-voxel coords
% ft_plot_headshape(hs)
% 
% figure
% ft_plot_sens(data_bpfreq.grad, 'unit', 'mm')
% ft_plot_headshape(hs, 'unit', 'mm')
% ft_plot_vol(vol, 'unit', 'mm')
% ft_plot_axes([], 'unit', 'mm');
% 

% 
% cfg = [];
% cfg.anaparameter = 'anatomy';
% cfg.funparameter = 'brain';
% cfg.location = [0 0 60];
% ft_sourceplot(cfg, mri_resliced)
% 
% 
% 
% %% figure 5, MRI scalp surface and polhemus headshape
% 
% cfg = [];
% cfg.tissue = 'scalp';
% cfg.method = 'isosurface';
% cfg.numvertices = 10000;
% scalp = ft_prepare_mesh(cfg, mri_segmented);
% 
% figure
% ft_plot_mesh(scalp, 'facecolor', 'skin')
% lighting phong
% camlight left
% camlight right
% material dull
% alpha 0.5
% ft_plot_headshape(hs, 'vertexcolor', 'k');
% 
% 
% 
% 
% % hier stimmt irgend etwas nicht:
% location = [0 0 60];
% figure
% ft_plot_ortho(mri_resliced.anatomy, 'transform', mri_resliced.transform, 'location', location, 'intersectmesh', vol.bnd)
% 
% figure
% ft_plot_montage(mri_resliced.anatomy, 'transform', mri_resliced.transform, 'intersectmesh', vol.bnd)
% 

% 
% 
% hold on
% plot3(grid.pos(:,1),grid.pos(:,2),grid.pos(:,3),'.');
% ft_plot_vol(vol);
% ft_plot_sens(data_bpfreq.grad);
