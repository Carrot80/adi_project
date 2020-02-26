function [] = check_headposition_within_sensor (ListSubj, fieldtripPath, path2folder)

for ii = 1:length(ListSubj)

    dir_hs_files = dir([ListSubj(ii).folder filesep ListSubj(ii).name filesep path2folder filesep 'hs_file_*']);
    
    for kk=1:length(dir_hs_files)
        
        if exist([dir_hs_files(kk).folder filesep 'Neu_Like500_' dir_hs_files(kk).name(end) '.mat'], 'file')
            shape = ft_read_headshape([dir_hs_files(kk).folder filesep dir_hs_files(kk).name], 'format', '4d_hs');
            load([dir_hs_files(kk).folder filesep 'Neu_Like500_' dir_hs_files(kk).name(end) '.mat'], 'cleanMEG')
            figure
            ft_plot_headshape(shape)
            hold on 
            plot3(cleanMEG.grad.chanpos(1:248,1), cleanMEG.grad.chanpos(1:248,2), cleanMEG.grad.chanpos(1:248,3), 'o')      
            ax = gca;
            set(gca, 'CameraPosition', [-0.242 -2.599 0.105]);  
            title([ListSubj(ii).name ' run_' dir_hs_files(kk).name(end) ])
            savefig ([dir_hs_files(kk).folder filesep dir_hs_files(kk).name '.fig'])
%             print(gcf, '-append','-dpsc2', ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\headpositions.ps']);
            close
            clear shape cleanMEG
        end
        
        
        
    end


 
end
 

 
 
 
 
 