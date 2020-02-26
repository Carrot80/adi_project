function check_number_trials_testset(subjects_dir, balldesign, path2file)
% modified: 13.2.20

fsample = 256.001;
time = -0.5:(1/fsample):1;
pp = 1;
for ii = 1:length(subjects_dir)  
    load([subjects_dir(ii).folder filesep subjects_dir(ii).name filesep path2file 'performance_real_data.mat'], 'perf')
    design_testset = {};
    for kk = 1:length(perf.CV)
        design_testset{kk} = perf.CV(kk).design;
    end
 ii
    for kk = 1:length(design_testset)
        if isfield(perf, design_testset{kk}) 
            trlno_testset(pp,1) = perf.(design_testset{kk}).number_of_trials.testset; 
            like_dislike_ratio_testset(pp,1) = numel(find(perf.CV(kk).labels_testset==1))./(numel(find(perf.CV(kk).labels_testset==1))+numel(find(perf.CV(kk).labels_testset==2)));
            if 1==isfield(perf.CV, 'trialnumer_likes_trainingsset')
                like_dislike_ratio_trainingset(pp,1) = perf.CV(kk).trialnumer_likes_trainingsset./(perf.CV(kk).trialnumer_likes_trainingsset+perf.CV(kk).trialnumer_dislikes_trainingsset);
            else
                like_dislike_ratio_trainingset(pp,1) = perf.CV(kk).trialnumber_likes_trainingsset./(perf.CV(kk).trialnumber_likes_trainingsset+perf.CV(kk).trialnumber_dislikes_trainingsset);
            end
            accuracy_prestim(pp,1) = nanmean(perf.(design_testset{kk}).lda.accuracy(1:129));
            accuracy_poststim_72ms(pp,1) = nanmean(perf.(design_testset{kk}).lda.accuracy(nearest(time, 0.072)));
        else
            trlno_testset(pp,1) = NaN; 
            like_dislike_ratio_testset(pp,1) = NaN;
            accuracy_prestim(pp,1) = NaN;
            like_dislike_ratio_trainingset(pp,1) = NaN;
            accuracy_poststim_72ms(pp,1) = NaN;
        end
        subject_no{pp,1} = subjects_dir(ii).name;
        design{pp,1} = design_testset{kk};        
        pp = pp+1;
    end
end

like_dislike_ratio = like_dislike_ratio_testset.(design_testset{1})';
for kk = 2:length(design_testset)
    like_dislike_ratio = cat(1, like_dislike_ratio, like_dislike_ratio_testset.(design_testset{kk})');
end

mean_accuracy_prestim = accuracy_prestim.(design_testset{1})';
for kk = 2:length(design_testset)
    mean_accuracy_prestim = cat(1, mean_accuracy_prestim, accuracy_prestim.(design_testset{kk})');
end

mean_accuracy_poststim_72ms = accuracy_poststim_72ms.(design_testset{1})';
for kk = 2:length(design_testset)
    mean_accuracy_poststim_72ms = cat(1, mean_accuracy_poststim_72ms, accuracy_poststim_72ms.(design_testset{kk})');
end

trlno = trlno_testset.(design_testset{1})';
for kk = 2:length(design_testset)
    trlno = cat(1, trlno, trlno_testset.(design_testset{kk})');
end


[rho,pval] = corr(trl_testset,acc, 'type', 'kendall')
[rho,pval] = corr(trl_testset,acc, 'type', 'spearman')

[rho,pval] = corr(like_ratio_training,acc, 'type', 'kendall')
[rho,pval] = corr(like_ratio_training,acc, 'type', 'spearman')

% Zusammenhang bildlich darstellen:
scatter(trlno, mean_accuracy_prestim)
scatter(like_ratio_training, acc)

[p,tbl,stats] = anova1(acc,design)

end