function [] = adi_mvpa_accuracy_significance_test(subjects_dir, path2file, balldesign)
% Description:  tests for significance of accuracy, dprime, and balanced accuracy values using permution
% 

% Created: 21.01.2020
% Latest update:


%--------------------------------------------------------------------------

%%











%%
fsample = 256;
sample_steps = 1/fsample;
time = -0.5:sample_steps:1;


for ii = 3:length(subjects_dir)
    load([subjects_dir(ii).folder filesep subjects_dir(ii).name filesep path2file ])
    balldesigns_enough_trls_testset = ones(1, length(balldesign));
    for pp = 1:length(balldesign)
        if numel(perf.ratings.(balldesign{pp})) <= 1 % if too few trials in testset (<=6), classification results are not used
            perf.lda.accuracy(pp,:) = NaN;
            balldesigns_enough_trls_testset(pp) = 0; 
        end
    end

    
    lda_acc_all(ii,:) = nanmean(perf.lda.accuracy);
    for kk = find(balldesigns_enough_trls_testset)
        num_of_trials.(subjects_dir(ii).name).(perf.CV(kk).design) = perf.number_of_trials{kk};    
    end
    perf_all(ii).accuracy = perf.lda.accuracy; 
    lda_num_of_trials(ii,:) = cell2mat(perf.number_of_trials);
    lda_num_of_trials(ii,find(balldesigns_enough_trls_testset==0)) = NaN;

% ....

%% using prestim-interval for significance test:
prestim = time(1:nearest(time, 0)-1);
% find 95. percentil:
perctl_95 = prctile(perf.lda.accuracy(1, 1:numel(prestim)),95);
perctl_99 = prctile(perf.lda.accuracy(1, 1:numel(prestim)),99);
perctl_999 = prctile(perf.lda.accuracy(1, 1:numel(prestim)),100);
indx_perctl_95 = find(perf.lda.accuracy(1, numel(prestim)+1:size(perf.lda.accuracy,2)) > perctl_95);
indx_perctl_99 = find(perf.lda.accuracy(1, numel(prestim)+1:size(perf.lda.accuracy,2)) > perctl_99);
indx_perctl_999 = find(perf.lda.accuracy(1, numel(prestim)+1:size(perf.lda.accuracy,2)) > perctl_999);

figure;
plot(time, perf.lda.accuracy(1,:))
hold on
plot(time(indx_perctl_95+128), ones(1,numel(indx_perctl_95)), '*')
hold on
plot(time, perctl_95*ones(1,numel(time)), 'g')

figure;
plot(time, perf.lda.accuracy(1,:))
hold on
plot(time(indx_perctl_99+128), ones(1,numel(indx_perctl_99)), '*')
hold on
plot(time, perctl_99*ones(1,numel(time)), 'g')

%% mean of prestim-interval: 
% problem: does not take standard error into account
mean_prestim = mean(perf.lda.accuracy(:, 1:numel(time(1:nearest(time, 0)-1))));
perctl_99_mean_prestim = prctile(mean_prestim,99);
mean_poststim = mean(perf.lda.accuracy(:, numel(mean_prestim)+1:end));
indx_perctl_99_poststim = find(mean_poststim > perctl_99_mean_prestim);

figure;
plot(time, [mean_prestim mean_poststim])
hold on
plot(time(indx_perctl_99_poststim+128), ones(1,numel(indx_perctl_99_poststim)), '*')
hold on
plot(time, perctl_99_mean_prestim*ones(1,numel(time)), 'Color', [0.5,0.5,0.5])

end











end