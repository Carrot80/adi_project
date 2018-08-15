    

function [] = adi_explore_virtsens(subjpath, path2atlas, freqband)
  
load (path2atlas)
load ([subjpath 'runs_appended\virtsens\' 'virtsens_all_conditions_allRuns_' freqband '.mat'])

path2save = [subjpath 'sanity_check\explore_virtsens\'];
if ~exist(path2save, 'dir')
    mkdir(path2save)
end

cfg              = [];
cfg.parameter    = 'trial';
cfg.keeptrials   = 'yes'; % classifiers operate on individual trials
cfg.vartrllength = 2;
conditions = fields(vs_allRuns);

for i = 1:size(conditions,1)
    tcondition = ft_timelockanalysis(cfg, vs_allRuns.(conditions{i}));
    figure;
    plot(tcondition.time, tcondition.avg)
    title(['virtsens_avg ' conditions{i}])
    savefig([path2save, 'avg_virtsens ' conditions{i} '_' freqband '.fig'])
    close
    % find bad virtsens: für jeden time sample soll er Außreißer finden 
    for k = 1:size(tcondition.avg,2)
        for p = 1:size(tcondition.avg,1)
            z_scores(p,k) = (tcondition.avg(p,k) - mean(tcondition.avg(:,k)))/std(tcondition.avg(:,k));
        end
    end
    figure;
    plot(tcondition.time, z_scores)
    hold on
    plot([-0.6 1.1], [3 3], 'k-')
    title(['z_scores_virtsens_avg ' conditions{i}])
    savefig([path2save, 'z_scores_virtsens_avg ' conditions{i} '_' freqband '.fig'])
    [ind_extrem_values_row, ind_extrem_values_column] = find (z_scores >=3);
    extreme_voxels = unique(ind_extrem_values_row);
    
%     extrem_values = tcondition.avg(ind_extrem_values_row, ind_extrem_values_column);
%     figure    
%     plot(tcondition.time, tcondition.avg(extreme_voxels, :))
%     index_atlas2 = find(atlas_downsampled.tissue); % !!!!
    for k=1:size(extreme_voxels,1)
        extreme_values.(conditions{i}).roi{k,1} = atlas_downsampled.sources_roi_numbers{extreme_voxels(k),2};
        extreme_values.(conditions{i}).roi{k,2} = extreme_voxels(k);
    end
        
    save([path2save 'extreme_values' freqband], 'extreme_values')    
        
        
end



%%
% nonzero = find(indx_atlas);
% 
% for k = 1:length(extrem_values)
%     num(k) = indx_atlas(extrem_values(k));
%     
% end
% % 
% num = index_atlas2(960); %!!!
% [~,b] = find(indx_atlas == num);
% atlas_downsampled.tissuelabel(b)
% zscore(415)
% 
% 
% indx = find(template_grid.tissue==x); % where x is the number of your choice
%  

end

