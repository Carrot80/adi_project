    

function main ()

subjmainpath      = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
ListSubj = dir(subjmainpath);
ListSubj(1:2) = [];

 for i = 1 : length(ListSubj) 
     
     subjpath = (['L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\'  ListSubj(i).name filesep 'MEG\sourcespace\noRois\']); 
     explore_virtsens(subjpath)
     
 end
 
end


function [] = explore_virtsens(subjpath)
  
load ('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\atlas_source_indices.mat');
path_adi_04 = 'L:\Arbeit\Adidas\fieldtrip_Auswertung\single_subjects\nl_adi_04\MEG\sourcespace\noROIs\runs_appended\virtsens\';
load ([path_adi_04 'virtsens_all_conditions_allRuns_bp1-45Hz.mat'])

cfg              = [];
cfg.parameter    = 'trial';
cfg.keeptrials   = 'yes'; % classifiers operate on individual trials
cfg.vartrllength = 2;
tdontcare           = ft_timelockanalysis(cfg, vs_allRuns.dontcare);
tlike           = ft_timelockanalysis(cfg, vs_allRuns.like); 
tdislike           = ft_timelockanalysis(cfg, vs_allRuns.dislike); 
figure
plot(tdontcare.time, tdontcare.avg)
figure
plot(tlike.time, tlike.avg)
figure
plot(tdislike.time, tdislike.avg)
    
% find bad virtsens:
% für jeden time sample soll er Außreißer finden (maximum)
for k = 1:size(tlike.avg,2)
    for p = 1:size(tlike.avg,1)
        z_score(p,k) = (tlike.avg(p,k) - mean(tlike.avg(:,k)))/std(tlike.avg(:,k));
    end
end
    
% z_score_all = zscore(tlike.avg); the same
figure
plot(tdislike.time, z_score_all)
hold on
plot([-0.6 1.1], [2 2], 'k-')
[ind_extrem_values_row, ind_extrem_values_column] = find (z_score_all >=3);

find(ind_extrem_values_column)

extreme_voxels = unique(ind_extrem_values_row);
for i= length(unique(ind_extrem_values_row))
    
    ind(i,:)
end

extrem_values = tdislike.avg(ind_extrem_values_row, ind_extrem_values_column);
figure    
plot(tdislike.time, tdislike.avg(extreme_voxels, :))




%%
% nonzero = find(indx_atlas);
index_atlas2 = find(atlas_downsampled.tissue); % !!!!
for k = 1:length(extrem_values)
    num(k) = indx_atlas(extrem_values(k));
    
end

num = index_atlas2(960); %!!!
[~,b] = find(indx_atlas == num);
atlas_downsampled.tissuelabel(b)
zscore(415)


indx = find(template_grid.tissue==x); % where x is the number of your choice
 

end

