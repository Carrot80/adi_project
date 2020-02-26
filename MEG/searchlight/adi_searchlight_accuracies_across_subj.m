function [] = searchlight_accuracies_across_subj(subjects_dir, path2data, filename, comp)

for ii = 2:length(subjects_dir)

    load([subjects_dir(ii).folder filesep subjects_dir(ii).name filesep path2data filesep filename])

    lda_acc(ii,:) = perf.lda.mean_accuracy;
    lda_num_of_trials(ii) = perf.number_of_trials;

end

weights = lda_num_of_trials/sum(lda_num_of_trials);

for kk = 1:size(lda_acc,2)
    weighted_mean_lda(kk) = sum(lda_acc(:,kk) .* weights');
end


stat = [];
stat.time = comp(1);
stat.mvpa.perf = weighted_mean_lda;
stat.accuracy = weighted_mean_lda;
stat.label = perf.features;
stat.dimord = 'time';
title(' avg searchlight LDA') 
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
% cfg.marker = 'labels';
ft_topoplotER(cfg, stat);
savefig(['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\searchlight\' 'lda_all_subjects_equal_number_trials_x2' num2str(comp(1)) '_' num2str(comp(2)) 'ms_weighted.fig'])

%% 2. Abbildung : standardabweichung

std_lda_acc = std(lda_acc);

stat = [];
stat.time = comp(1);
stat.mvpa.perf = std_lda_acc;
stat.accuracy = std_lda_acc;
stat.label = perf.features;
stat.dimord = 'time';
title(' avg searchlight LDA') 
cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = '4D248_helmet.mat';            
%    cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
ft_topoplotER(cfg, stat);
savefig(['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\searchlight\' 'lda_all_subjects_equal_number_trials_x2' num2str(comp(1)) '_' num2str(comp(2)) 'ms_std.fig'])


max(lda_acc')

%%

% proband nr. 04 und 21 rausnehmen, da nur 1 run

end