function [] = adi_leave_out_exemplar_mean_subj(subjects_dir, path2file)
%% computes mean accuracy across subjects for classification with channelselection (results from searchlight)
% includes calculation of MCC, Sensitivity, False Alarm rate (not done yet)
% includes calculation of weighted mean and weighted CI and figures
%% bei einigen nochmal berechnen, bei denen im Nachhinein runs entfernt worden ==> hier evtl. pro run einzeln berechnen
%% Nachberechnung: bei einigen Bällen sind zu wenige trials in Testset => accuracy-Werte müssen an der Stelle entfernt werden bzw. durch NaNs ersetzt werden
fsample = 256;
sample_steps = 1/fsample;
time = -0.5:sample_steps:1;
balldesign = {'gbf';'gbs';'gbv';'ggf';'ggs';'ggv';'rwf';'rws';'rwv'};
subjects_dir([1 17 19 26],:) = [];
%
for ii = 1:length(subjects_dir)
    load([subjects_dir(ii).folder filesep subjects_dir(ii).name filesep path2file ])
    balldesigns_enough_trls_testset = ones(1, length(balldesign));
    for pp = 1:length(balldesign)
        if numel(perf.ratings.(balldesign{pp})) <= 6 % if too few trials in testset (<=6), classification results are not used
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
%     for pp = 1:length(perf.lda.rwv)
%         
%         C = confusionmat(perf.lda.rwv(pp).testlabels, perf.lda.rwv(pp).predicted_labels);
%         try
%             hits_sensitivity_TPR(pp) = C(1,1)/(C(1,1)+C(1,2));
%         catch
%             hits_sensitivity_TPR
%         end
%         
%         corr_rejected_specificity(pp) = C(2,2)/(C(2,2)+C(2,1));
%         false_alarms_FPR(pp) = C(2,1)/(C(2,1)+C(2,2));
%         miss_FNR(pp) = C(1,2)/(C(1,2)+C(1,1));
%         
%     
%         %% Matthews correlation coefficient:
%         % MCC = (TP*TN-FP*FN)/((sqrt(TP+FP)*(FP+FN)*(TN+FP)*(TN+FN))
%         TP = C(1,1);
%         TN = C(2,2);
%         FP = C(2,1);
%         FN = C(1,2);
%         MCC(pp) = (TP*TN-FP*FN)/(sqrt((TP+FP)*(FP+FN)*(TN+FP)*(TN+FN)));
%         accuracy =  (TP + TN)/(TP + TN + FP + FN)
%         clear C
%     end
   
%     figure
%     plot(time, perf.lda.accuracy(1,:))
%     hold on
%     plot(time, MCC)
%     
%     figure
%     plot(time, perf.lda.accuracy(1,:))
%     hold on
%     plot(time, hits_sensitivity_TPR)
%     hold on
%     plot(time, corr_rejected_specificity)
%     hold on
%     plot(time, false_alarms_FPR)
%     ylim([-0.1 1.1])
    clear perf
    close all
end

%% evtl. anstelle der Accuracy auch Hits und false alarm rate bzw. AUC? (Yuval)
% hit = correctly predicted like/number of actual likes (#Hit/(#Hit+#Miss)
% correct_rejected = correctly predicted dislikes/ number of acutal dislikes
% false_alarm: incorrectly predicted likes/number of acutal dislikes
% (#FA/(#FA/)
% miss: incorrectly predicted dislikes/number of acutal likes
% F-Score:
% what could be the changing criterion in ROC?
% mit ROC measures nicht weitergemacht
figure
plot(time,mean(lda_acc_all))
trials_per_subj = nansum(lda_num_of_trials');
sum_all_trials = sum(trials_per_subj);
weights = trials_per_subj./sum_all_trials;

for kk = 1:size(lda_acc_all,2)
    weighted_mean_lda(kk) = sum(lda_acc_all(:,kk) .* weights');
end
% standard deviation: the standard error of the weighted mean is sqrt( (w1*s1)^2 + (w2*s2)^2 + (w3*s3)^2). 

for kk = 1:size(lda_acc_all,2)
    for pp = 1:length(lda_acc_all(:,pp))
        weighted_acc_lda(pp, kk) = lda_acc_all(pp,kk) .* weights(pp);
        weighted_acc_lda2(pp, kk) = weighted_acc_lda(pp, kk)^2;
    end
end
std_weighted_temp = sum(weighted_acc_lda2);
std_weighted = sqrt(std_weighted_temp);

% confidence interval

for kk = 1:length(weighted_mean_lda)
    CI_weighted(kk,:) = ci_weighted(weighted_mean_lda(kk), std_weighted(kk), size(lda_acc_all,1));
end

figure;
plot(time, weighted_mean_lda)
xlabel('time')
ylabel('accuracy')
box off
hold on
plot(time, CI_weighted(:,1), 'k:')
hold on
plot(time, CI_weighted(:,2), 'k:')
title('weighted mean with 95% CI')

figure;
plot(time, mean(weighted_mean_lda'), 'b')
hold on
plot(time, CI(:,1), 'k:')
hold on
plot(time, CI(:,2), 'k:')
box off; xlabel('time'); ylabel('accuracy lda')
title('weighted mean')



% 
% 
% for kk = 1:length(perf_all)
% 
%     for jj = 1:length(balldesign)
%         accuracy_all.(balldesign{jj})(kk,:) = perf_all(kk).accuracy(jj,:);
%     end
% end
% 
% %% Favoritenrangfolge
% for kk = 1:length(subjects_dir)
%     if 1 == strcmp(subjects_dir(kk).name, list_favorite_balls(kk,1))
%         for pp = 1:4
%             fav_balls(pp) = list_favorite_balls(kk, pp+1);
%         end
%         
%         for jj = 1:length(balldesign)
%             if ~isempty(find(strcmp(fav_balls, balldesign{jj})))
%                 rank_balldesign.(balldesign{jj})(kk,1) = find(strcmp(fav_balls, balldesign{jj}));
%             else 
%                 rank_balldesign.(balldesign{jj})(kk,1) = 5;
%             end
%         end
% 
%     end
%     clear fav_balls
% end
% 
% 
% 
% 
% 
% 
% 
% %% gbf:
% mean_gbf = mean(accuracy_all.gbf,1);
% for kk = 1:length(accuracy_all.gbf)
%     CI_gbf(kk,:) = ci(accuracy_all.gbf(:,kk));
% end
% figure; hold on;
% plot(time, mean_gbf)
% plot(time, CI_gbf(:,1), 'k:')
% hold on
% plot(time, CI_gbf(:,2), 'k:')
% title('mean accuracy gbf')
% 
% %% gbs:
% mean_gbs = mean(accuracy_all.gbs,1);
% for kk = 1:length(accuracy_all.gbs)
%     CI_gbs(kk,:) = ci(accuracy_all.gbs(:,kk));
% end
% figure; hold on;
% plot(time, mean_gbs)
% plot(time, CI_gbs(:,1), 'k:')
% hold on
% plot(time, CI_gbs(:,2), 'k:')
% title('mean accuracy gbs')
% 
% %% gbv
% mean_gbv = mean(accuracy_all.gbv,1);
% for kk = 1:length(accuracy_all.gbv)
%     CI_gbv(kk,:) = ci(accuracy_all.gbv(:,kk));
% end
% figure; hold on;
% plot(time, mean_gbv)
% plot(time, CI_gbv(:,1), 'k:')
% hold on
% plot(time, CI_gbv(:,2), 'k:')
% title('mean accuracy gbv')
% 
% %% ggf
% 
% mean_ggf = mean(accuracy_all.ggf,1);
% for kk = 1:length(accuracy_all.ggf)
%     CI_ggf(kk,:) = ci(accuracy_all.ggf(:,kk));
% end
% figure; hold on;
% plot(time, mean_ggf)
% plot(time, CI_ggf(:,1), 'k:')
% hold on
% plot(time, CI_ggf(:,2), 'k:')
% title('mean accuracy ggf')
% 
% 
% %% ggs:
% mean_ggs = mean(accuracy_all.ggs,1);
% for kk = 1:length(accuracy_all.ggs)
%     CI_ggs(kk,:) = ci(accuracy_all.ggs(:,kk));
% end
% figure; hold on;
% plot(time, mean_ggs)
% plot(time, CI_ggs(:,1), 'k:')
% hold on
% plot(time, CI_ggs(:,2), 'k:')
% title('mean accuracy ggs')
% 
% 
% %% ggv:
% 
% mean_ggv = mean(accuracy_all.ggv,1);
% for kk = 1:length(accuracy_all.ggv)
%     CI_ggv(kk,:) = ci(accuracy_all.ggv(:,kk));
% end
% figure; hold on;
% plot(time, mean_ggv)
% plot(time, CI_ggv(:,1), 'k:')
% hold on
% plot(time, CI_ggv(:,2), 'k:')
% title('mean accuracy ggv')
% 
% 
% %% rwf
% mean_rwf = mean(accuracy_all.rwf,1);
% for kk = 1:length(accuracy_all.rwf)
%     CI_rwf(kk,:) = ci(accuracy_all.rwf(:,kk));
% end
% figure; hold on;
% plot(time, mean_rwf)
% plot(time, CI_rwf(:,1), 'k:')
% hold on
% plot(time, CI_rwf(:,2), 'k:')
% title('mean accuracy rwf')
% 
% %% rws
% mean_rws = mean(accuracy_all.rws,1);
% for kk = 1:length(accuracy_all.rws)
%     CI_rws(kk,:) = ci(accuracy_all.rws(:,kk));
% end
% figure; hold on;
% plot(time, mean_rws)
% plot(time, CI_rws(:,1), 'k:')
% hold on
% plot(time, CI_rws(:,2), 'k:')
% title('mean accuracy rws')
% 
% 
% %% rwv
% 
% mean_rwv = mean(accuracy_all.rwv,1);
% for kk = 1:length(accuracy_all.rwv)
%     CI_rwv(kk,:) = ci(accuracy_all.rwv(:,kk));
% end
% figure; hold on;
% plot(time, mean_rwv)
% plot(time, CI_rwv(:,1), 'k:')
% hold on
% plot(time, CI_rwv(:,2), 'k:')
% title('mean accuracy rwv')
% 
% 
% 
% for kk = 1:length(lda_acc_all)
%     CI(kk,:) = ci(lda_acc_all(:,kk));
% end
% 
% 
% 
% fsample = 256;
% sample_steps = 1/fsample;
% time = -0.5:sample_steps:1;
% 
% figure;
% plot(time, mean(lda_acc_all), 'k')
% hold on
% plot(time, CI(:,1), 'k:')
% hold on
% plot(time, CI(:,2), 'k:')
% box off; xlabel('time'); ylabel('accuracy lda')
% 
% 
% 
% 
% upper = CI(:,1);
% lower = CI(:,2);
% x = time;
% if find(size(x)==(max(size(x))))<2
% x=x'; end
% if find(size(lower)==(max(size(lower))))<2
% lower=lower'; end
% if find(size(upper)==(max(size(upper))))<2
% upper=upper'; end
% 
% 
% colour='b';
% figure;
% fill([x fliplr(x)],[upper fliplr(lower)],colour)
% 
% if length(lower)~=length(upper)
%     error('lower and upper vectors must be same length')
% end
% 
% 
% 
% hold on; plot(time, mean(lda_acc_all), 'k')
% 
% 
% 
% 
% 
% 




end

%% subfunctions
function CI = ci_weighted(x_weighted, std_weighted, N)

SEM = std_weighted/sqrt(N);               % Standard Error
ts = tinv([0.025  0.975],N-1);      % T-Score
CI = x_weighted + ts*SEM;     
end