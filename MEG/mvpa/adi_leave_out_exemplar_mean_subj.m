function [] = adi_leave_out_exemplar_mean_subj(subjects_dir, path2file)

%% bei einigen nochmal berechnen, bei denen im Nachhinein runs entfernt worden ==> hier evtl. pro run einzeln berechnen
%% zur Berechnung einer Regression, um herauszufinden, welche Faktoren eine gute bzw. schlechte Klassifikation bedingen
list_favorite_balls = import_favouriteballlist('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\list_favourite_balls.txt');

fsample = 256;
sample_steps = 1/fsample;
time = -0.5:sample_steps:1;
samples_beg = nearest(time, 0.56);
samples_end = nearest(time, 0.740);

balldesign = {'gbf';'gbs';'gbv';'ggf';'ggs';'ggv';'rwf';'rws';'rwv'};

for ii = 1:length(subjects_dir)

    load([subjects_dir(ii).folder filesep subjects_dir(ii).name filesep path2file ])
    lda_acc_all(ii,:) = mean(perf.lda.accuracy);
    for kk = 1:length(perf.CV)
        num_of_trials.(subjects_dir(ii).name).(perf.CV(kk).design) = perf.number_of_trials{kk};    
    end
    perf_all(ii).accuracy = perf.lda.accuracy;
    
    
%     figure;
%     plot(time, perf.lda.accuracy)
%     legend({'gbf';'gbs';'gbv';'ggf';'ggs';'ggv';'rwf';'rws';'rwv'})
%     title([subjects_dir(ii).name])
%     savefig([subjects_dir(ii).folder filesep subjects_dir(ii).name filesep 'MEG_analysis\noisereduced\1_95Hz\mvpa\perf_without_pca_bootstrap_pseudotrials_trainfold_separate_balldesigns.fig' ])
   
end
close all


for ii = 1:size(struct2cell(num_of_trials),1)
    for kk = 1:length(balldesign)

%         trials_all.(balldesign{kk})(ii,:) = num_of_trials.(subjects_dir(ii).name).(balldesign{kk});
        mean_accuracy_comp1.(balldesign{kk})(ii,1) = mean(perf_all(ii).accuracy(kk,samples_beg:samples_end));    
    end
end


for ii = 1:size(struct2cell(num_of_trials),1)
    for kk = 1:length(balldesign)

        trials_all.(balldesign{kk})(ii,:) = num_of_trials(ii).(balldesign(kk));
        mean_accuracy_30_156ms.(balldesign{kk})(ii,1) = mean(perf_all(ii).accuracy(kk,samples_30ms:samples_150ms));    
    end
end

for kk = 1:length(perf_all)

    for jj = 1:length(balldesign)
        accuracy_all.(balldesign{jj})(kk,:) = perf_all(kk).accuracy(jj,:);
        mean_accuracy_30_156ms.(balldesign{jj})(kk,1) = mean(perf_all(kk).accuracy(jj,samples_30ms:samples_150ms));    
    end
end

%% Favoritenrangfolge
for kk = 1:length(subjects_dir)
    if 1 == strcmp(subjects_dir(kk).name, list_favorite_balls(kk,1))
        for pp = 1:4
            fav_balls(pp) = list_favorite_balls(kk, pp+1);
        end
        
        for jj = 1:length(balldesign)
            if ~isempty(find(strcmp(fav_balls, balldesign{jj})))
                rank_balldesign.(balldesign{jj})(kk,1) = find(strcmp(fav_balls, balldesign{jj}));
            else 
                rank_balldesign.(balldesign{jj})(kk,1) = 5;
            end
        end

    end
    clear fav_balls
end







%% gbf:
mean_gbf = mean(accuracy_all.gbf,1);
for kk = 1:length(accuracy_all.gbf)
    CI_gbf(kk,:) = ci(accuracy_all.gbf(:,kk));
end
figure; hold on;
plot(time, mean_gbf)
plot(time, CI_gbf(:,1), 'k:')
hold on
plot(time, CI_gbf(:,2), 'k:')
title('mean accuracy gbf')

%% gbs:
mean_gbs = mean(accuracy_all.gbs,1);
for kk = 1:length(accuracy_all.gbs)
    CI_gbs(kk,:) = ci(accuracy_all.gbs(:,kk));
end
figure; hold on;
plot(time, mean_gbs)
plot(time, CI_gbs(:,1), 'k:')
hold on
plot(time, CI_gbs(:,2), 'k:')
title('mean accuracy gbs')

%% gbv
mean_gbv = mean(accuracy_all.gbv,1);
for kk = 1:length(accuracy_all.gbv)
    CI_gbv(kk,:) = ci(accuracy_all.gbv(:,kk));
end
figure; hold on;
plot(time, mean_gbv)
plot(time, CI_gbv(:,1), 'k:')
hold on
plot(time, CI_gbv(:,2), 'k:')
title('mean accuracy gbv')

%% ggf

mean_ggf = mean(accuracy_all.ggf,1);
for kk = 1:length(accuracy_all.ggf)
    CI_ggf(kk,:) = ci(accuracy_all.ggf(:,kk));
end
figure; hold on;
plot(time, mean_ggf)
plot(time, CI_ggf(:,1), 'k:')
hold on
plot(time, CI_ggf(:,2), 'k:')
title('mean accuracy ggf')


%% ggs:
mean_ggs = mean(accuracy_all.ggs,1);
for kk = 1:length(accuracy_all.ggs)
    CI_ggs(kk,:) = ci(accuracy_all.ggs(:,kk));
end
figure; hold on;
plot(time, mean_ggs)
plot(time, CI_ggs(:,1), 'k:')
hold on
plot(time, CI_ggs(:,2), 'k:')
title('mean accuracy ggs')


%% ggv:

mean_ggv = mean(accuracy_all.ggv,1);
for kk = 1:length(accuracy_all.ggv)
    CI_ggv(kk,:) = ci(accuracy_all.ggv(:,kk));
end
figure; hold on;
plot(time, mean_ggv)
plot(time, CI_ggv(:,1), 'k:')
hold on
plot(time, CI_ggv(:,2), 'k:')
title('mean accuracy ggv')


%% rwf
mean_rwf = mean(accuracy_all.rwf,1);
for kk = 1:length(accuracy_all.rwf)
    CI_rwf(kk,:) = ci(accuracy_all.rwf(:,kk));
end
figure; hold on;
plot(time, mean_rwf)
plot(time, CI_rwf(:,1), 'k:')
hold on
plot(time, CI_rwf(:,2), 'k:')
title('mean accuracy rwf')

%% rws
mean_rws = mean(accuracy_all.rws,1);
for kk = 1:length(accuracy_all.rws)
    CI_rws(kk,:) = ci(accuracy_all.rws(:,kk));
end
figure; hold on;
plot(time, mean_rws)
plot(time, CI_rws(:,1), 'k:')
hold on
plot(time, CI_rws(:,2), 'k:')
title('mean accuracy rws')


%% rwv

mean_rwv = mean(accuracy_all.rwv,1);
for kk = 1:length(accuracy_all.rwv)
    CI_rwv(kk,:) = ci(accuracy_all.rwv(:,kk));
end
figure; hold on;
plot(time, mean_rwv)
plot(time, CI_rwv(:,1), 'k:')
hold on
plot(time, CI_rwv(:,2), 'k:')
title('mean accuracy rwv')



for kk = 1:length(lda_acc_all)
    CI(kk,:) = ci(lda_acc_all(:,kk));
end



fsample = 256;
sample_steps = 1/fsample;
time = -0.5:sample_steps:1;

figure;
plot(time, mean(lda_acc_all), 'k')
hold on
plot(time, CI(:,1), 'k:')
hold on
plot(time, CI(:,2), 'k:')
box off; xlabel('time'); ylabel('accuracy lda')




upper = CI(:,1);
lower = CI(:,2);
x = time;
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end


colour='b';
figure;
fill([x fliplr(x)],[upper fliplr(lower)],colour)

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end



hold on; plot(time, mean(lda_acc_all), 'k')










end