function [] = adi_leave_out_exemplar_mean_subj(subjects_dir, path2file)

%% bei einigen nochmal berechnen, bei denen im Nachhinein runs entfernt worden ==> hier pro run einzeln berechnen
list_favorite_balls = import_favouriteballlist('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\list_favourite_balls.txt');
% entferne Probanden mit nur einem run
subjects_dir([1 17 19 26],:) = [];
list_favorite_balls([1 17 19 26],:) = [];

for ii = 1:length(subjects_dir)

    load([subjects_dir(ii).folder filesep subjects_dir(ii).name filesep path2file ])
    perf_all(ii) = perf;
    lda_acc_all(ii,:) = mean(perf.lda.accuracy);
    lda_num_of_trials(ii) = perf.number_of_trials;

end

weights = lda_num_of_trials/sum(lda_num_of_trials);

for kk = 1:size(lda_acc_all,2)
    weighted_mean_lda(kk) = sum(lda_acc_all(:,kk) .* weights');
end



for kk = 1:length(lda_acc_all)
    CI(kk,:) = ci(lda_acc_all(:,kk));
end


fsample = 256;
sample_steps = 1/fsample;
time = -0.5:sample_steps:1;

figure;
plot(time,weighted_mean_lda, 'b')
hold on
plot(time, CI(:,1), 'k:')
hold on
plot(time, CI(:,2), 'k:')
box off; xlabel('time'); ylabel('accuracy lda')
title('weighted mean')

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


%% performance of favorite_ball no 1:

fav_ball = list_favorite_balls(:,2);
 

for kk = 1:length(fav_ball)
    
    for ll = 1:length(perf_all(kk).crossval_balldesign)
        order_CV{ll} = perf_all(kk).crossval_balldesign(ll).design;
    end
    ind_ = [];
    ind_ = find(strcmp(order_CV, fav_ball(kk)));

    perf_fav1(kk,:) = perf_all(kk).lda.accuracy(ind_, :);
    
end

%% Mittelung und Gewichtung je nach Anzahl der zugrundeliegenden Trials:

for kk = 1:size(perf_fav1,2)
    weighted_mean_Fav1(kk) = sum(perf_fav1(:,kk) .* weights');
end


for kk = 1:length(perf_fav1)
    CI_fav1(kk,:) = ci(perf_fav1(:,kk));
end

figure;
plot(time, weighted_mean_Fav1, 'b')
hold on
plot(time, CI_fav1(:,1), 'k:')
hold on
plot(time, CI_fav1(:,2), 'k:')
title('mean MVPA Favorite ball No 1 in testset');
box off
% hier weitermachen, evlt. noch CI gewichten

%% performance of favorite_ball no 2

fav_ball_no2 = list_favorite_balls(:,3);

for kk = 1:length(fav_ball_no2)
    ind_ = [];
    order_CV = [];
    for ll = 1:length(perf_all(kk).crossval_balldesign)
        order_CV{ll} = perf_all(kk).crossval_balldesign(ll).design;
    end
   
    ind_ = find(strcmp(order_CV, fav_ball_no2(kk)));

    perf_fav2(kk,:) = perf_all(kk).lda.accuracy(ind_, :);
    
end



for kk = 1:size(perf_fav2,2)
    weighted_mean_Fav2(kk) = sum(perf_fav2(:,kk) .* weights');
end


figure; hold on;
plot(time, weighted_mean_Fav1, 'k')
plot(time, weighted_mean_Fav2, 'b')
title('mean MVPA Favorite ball No 1 and 2');
legend({'Favorite 1'; 'Favorite 2'})

for kk = 1:length(perf_fav2)
    CI_Fav2(kk,:) = ci(perf_fav2(:,kk));
end




for kk = 1:size(perf_fav2,2)
    mean_Fav2(kk) = mean(perf_fav2(:,kk));
end

figure;
plot(time,weighted_mean_Fav2, 'b')
hold on
plot(time, mean_Fav2, 'k:')


figure;
plot(time,mean_Fav2, 'b')
hold on
plot(time, CI_Fav2(:,1), 'k:')
hold on
plot(time, CI_Fav2(:,2), 'k:')
box off; xlabel('time'); ylabel('accuracy lda')
title('mean Fav2')



%% performance of favorite_ball no 3 :

fav_ball_no3 = list_favorite_balls(:,4);

for kk = 1:length(fav_ball_no3)
    ind_ = [];
    order_CV = [];
    for ll = 1:length(perf_all(kk).crossval_balldesign)
        order_CV{ll} = perf_all(kk).crossval_balldesign(ll).design;
    end
   
    ind_ = find(strcmp(order_CV, fav_ball_no3(kk)));

    perf_fav3(kk,:) = perf_all(kk).lda.accuracy(ind_, :);
    
end


for kk = 1:size(perf_fav3,2)
    weighted_mean_Fav3(kk) = sum(perf_fav3(:,kk) .* weights');
end

figure; hold on;
plot(time, weighted_mean_Fav1, 'k')
plot(time, weighted_mean_Fav2, 'b')
plot(time, weighted_mean_Fav3, 'g')
title('mean MVPA Favorite ball No 1,2 and 3');
legend({'Favorite 1'; 'Favorite 2'; 'Favorite 3' })


%% performance of favorite_ball no 4 :

fav_ball_no4 = list_favorite_balls(:,5);

for kk = 1:length(fav_ball_no4)
    ind_ = [];
    order_CV = [];
    for ll = 1:length(perf_all(kk).crossval_balldesign)
        order_CV{ll} = perf_all(kk).crossval_balldesign(ll).design;
    end
   
    ind_ = find(strcmp(order_CV, fav_ball_no4(kk)));

    perf_fav4(kk,:) = perf_all(kk).lda.accuracy(ind_, :);
    
end


for kk = 1:size(perf_fav4,2)
    weighted_mean_Fav4(kk) = sum(perf_fav4(:,kk) .* weights');
end

figure; hold on;
plot(time, weighted_mean_Fav1, 'k')
plot(time, weighted_mean_Fav2, 'b')
plot(time, weighted_mean_Fav3, 'g')
plot(time, weighted_mean_Fav4, 'r')

title('mean MVPA Favorite ball No 1,2,3 and 4');
legend({'Favorite 1'; 'Favorite 2'; 'Favorite 3'; 'Favorite 4' })

%% mean performance of unpreferred balldesigns:

preferred_balls = list_favorite_balls(:,2:5);


for kk = 1:size(preferred_balls,1)
    ind_ = [];
    order_CV = [];    
    for ll = 1:length(perf_all(kk).crossval_balldesign)
        order_CV{ll} = perf_all(kk).crossval_balldesign(ll).design;
    end
    balls_fav1_4 = preferred_balls(kk,:);
    for pp = 1:length(balls_fav1_4)
        ind_(pp) = find(strcmpi(order_CV, balls_fav1_4(pp)));
    end
    ind_reverse = [1:9];
    ind_reverse(ind_) = [];
    perf_unpreferred(kk,:) = mean(perf_all(kk).lda.accuracy(ind_reverse, :));
    
end


for kk = 1:size(perf_unpreferred,2)
    weighted_mean_unpreferred(kk) = sum(perf_unpreferred(:,kk) .* weights');
end

figure; hold on;
plot(time, mean([weighted_mean_Fav1; weighted_mean_Fav2; weighted_mean_Fav3; weighted_mean_Fav4]), 'k')
plot(time, weighted_mean_unpreferred, 'b')
legend({'Favorite 1 to 4'; 'unpreferred' })



figure; hold on;
plot(time, weighted_mean_Fav1, 'k')
plot(time, weighted_mean_Fav2, 'b')
plot(time, weighted_mean_Fav3, 'g')
plot(time, weighted_mean_Fav4, 'r')
plot(time, weighted_mean_unpreferred, 'y')
legend({'Favorite 1'; 'Favorite 2'; 'Favorite 3'; 'Favorite4';'unpreferred balldesigns' })


Fav_2_4 = mean([weighted_mean_Fav2; weighted_mean_Fav3; weighted_mean_Fav4]);
figure; hold on;
plot(time, weighted_mean_Fav1, 'k')
plot(time, Fav_2_4, 'b')
legend({'Favorite 1'; 'Favorite 2 to 4' })

%%  separate performance of unpreferred balldesigns :

preferred_balls = list_favorite_balls(:,2:5);

perf_unpreferred = [];
for kk = 1:size(preferred_balls,1)
    ind_ = [];
    order_CV = [];    
    for ll = 1:length(perf_all(kk).crossval_balldesign)
        order_CV{ll} = perf_all(kk).crossval_balldesign(ll).design;
    end
    balls_fav1_4 = preferred_balls(kk,:);
    for pp = 1:length(balls_fav1_4)
        ind_(pp) = find(strcmpi(order_CV, balls_fav1_4(pp)));
    end
    ind_reverse = [1:9];
    ind_reverse(ind_) = [];
    
    for pp = ind_reverse
        perf_unpreferred(kk).(order_CV{pp}(1,:)) = perf_all(kk).lda.accuracy(pp, :);
    end
end

for pp = 1:length(perf_unpreferred)
    if ~isempty(perf_unpreferred(pp).gbv)
        perf_unpreferred_gbv(pp,:) = perf_unpreferred(pp).gbv;
    else
        perf_unpreferred_gbv(pp,1:385) = NaN; 
    end
    if ~isempty(perf_unpreferred(pp).gbs)
        perf_unpreferred_gbs(pp,:) = perf_unpreferred(pp).gbs;
    else
        perf_unpreferred_gbs(pp,1:385) = NaN;     
    end
    
    if ~isempty(perf_unpreferred(pp).ggs)
        perf_unpreferred_ggs(pp,:) = perf_unpreferred(pp).ggs;
    else
        perf_unpreferred_ggs(pp,1:385) = NaN;
    end
    
    if ~isempty(perf_unpreferred(pp).rws)
        perf_unpreferred_rws(pp,:) = perf_unpreferred(pp).rws;
    else
        perf_unpreferred_rws(pp,1:385) = NaN;        
    end
    
    if ~isempty(perf_unpreferred(pp).gbf)
        perf_unpreferred_gbf(pp,:) = perf_unpreferred(pp).gbf;
    else
        perf_unpreferred_gbf(pp,1:385) = NaN;        
    end
    
    if ~isempty(perf_unpreferred(pp).ggf)
        perf_unpreferred_ggf(pp,:) = perf_unpreferred(pp).ggf;
    else
        perf_unpreferred_ggf(pp,1:385) = NaN;        
    end
    
    if ~isempty(perf_unpreferred(pp).ggv)
        perf_unpreferred_ggv(pp,:) = perf_unpreferred(pp).ggv;
    else
        perf_unpreferred_ggv(pp,1:385) = NaN;        
    end
end



%gbv
ind_nan_gbv = find(isnan(perf_unpreferred_gbv(:,1)));
ind_true_gbv = 1:size(perf_unpreferred_gbv,1);
ind_true_gbv(ind_nan_gbv) = [];

for kk = 1:size(perf_unpreferred_gbv,2)
    weighted_mean_unpreferred_gbv(kk) = sum(perf_unpreferred_gbv(ind_true_gbv,kk) .* weights(ind_true_gbv)');
end

%gbs
ind_nan_gbs = find(isnan(perf_unpreferred_gbs(:,1)));
ind_true_gbs = 1:size(perf_unpreferred_gbs,1);
ind_true_gbs(ind_nan_gbs) = [];

for kk = 1:size(perf_unpreferred_gbs,2)
    weighted_mean_unpreferred_gbs(kk) = sum(perf_unpreferred_gbs(ind_true_gbs,kk) .* weights(ind_true_gbs)');
end

%ggs
ind_nan_ggs = find(isnan(perf_unpreferred_ggs(:,1)));
ind_true_ggs = 1:size(perf_unpreferred_ggs,1);
ind_true_ggs(ind_nan_ggs) = [];

for kk = 1:size(perf_unpreferred_ggs,2)
    weighted_mean_unpreferred_ggs(kk) = sum(perf_unpreferred_ggs(ind_true_ggs,kk) .* weights(ind_true_ggs)');
end

%rws
ind_nan_rws = find(isnan(perf_unpreferred_rws(:,1)));
ind_true_rws = 1:size(perf_unpreferred_rws,1);
ind_true_rws(ind_nan_rws) = [];

for kk = 1:size(perf_unpreferred_rws,2)
    weighted_mean_unpreferred_rws(kk) = sum(perf_unpreferred_rws(ind_true_rws,kk) .* weights(ind_true_rws)');
end

% gbf
ind_nan_gbf = find(isnan(perf_unpreferred_gbf(:,1)));
ind_true_gbf = 1:size(perf_unpreferred_gbf,1);
ind_true_gbf(ind_nan_gbf) = [];

for kk = 1:size(perf_unpreferred_gbf,2)
    weighted_mean_unpreferred_gbf(kk) = sum(perf_unpreferred_gbf(ind_true_gbf,kk) .* weights(ind_true_gbf)');
end

% ggf
ind_nan_ggf = find(isnan(perf_unpreferred_ggf(:,1)));
ind_true_ggf = 1:size(perf_unpreferred_ggf,1);
ind_true_ggf(ind_nan_ggf) = [];

for kk = 1:size(perf_unpreferred_ggf,2)
    weighted_mean_unpreferred_ggf(kk) = sum(perf_unpreferred_ggf(ind_true_ggf,kk) .* weights(ind_true_ggf)');
end




% ggv
ind_nan_ggv = find(isnan(perf_unpreferred_ggv(:,1)));
ind_true_ggv = 1:size(perf_unpreferred_ggv,1);
ind_true_ggv(ind_nan_ggv) = [];

for kk = 1:size(perf_unpreferred_ggv,2)
    weighted_mean_unpreferred_ggv(kk) = sum(perf_unpreferred_ggv(ind_true_ggv,kk) .* weights(ind_true_ggv)');
end

% perf_unpreferred_ggv = perf_unpreferred_ggv(ind_true_ggv,:);

figure; hold on;
plot(time, weighted_mean_unpreferred_ggv, 'b')
plot(time, weighted_mean_unpreferred_ggf, 'b')


figure; hold on;
plot(time, weighted_mean_Fav1, 'k')
plot(time, Fav_2_4, 'b')
plot(time, weighted_mean_unpreferred_gbs, 'g')
plot(time, weighted_mean_unpreferred_ggs, 'r')
plot(time, weighted_mean_unpreferred_gbv, 'y')
plot(time, weighted_mean_unpreferred_rws, '--')
legend({'Favorite 1'; 'Favorite 2 to 4';  'gbs (unpreferred)'; 'ggs (unpreferred)'; 'gbv (unpreferred)'; 'rws (unpreferred)'})

figure; hold on;
plot(time, weighted_mean_Fav1, 'k')
plot(time, Fav_2_4, 'b')
plot(time, weighted_mean_unpreferred_gbs, 'g')
plot(time, weighted_mean_unpreferred_ggs, 'r')
plot(time, weighted_mean_unpreferred_gbv, 'y')
plot(time, weighted_mean_unpreferred_rws, '--')
legend({'Favorite 1'; 'Favorite 2 to 4';  'gbs'; 'ggs'; 'gbv'; 'rws'})


figure; 
plot (time, perf_unpreferred_gbs(ind_true_gbs,:))

%% Konfidenzintervall unpreferred:
perf_unpreferred_rws = perf_unpreferred_rws(ind_true_rws,:);
perf_unpreferred_gbv = perf_unpreferred_gbv(ind_true_gbv,:);
perf_unpreferred_ggs = perf_unpreferred_ggs(ind_true_ggs,:);
perf_unpreferred_gbs = perf_unpreferred_gbs(ind_true_gbs,:);
perf_unpreferred_gbf = perf_unpreferred_gbf(ind_true_gbf,:);

for kk = 1:length(perf_unpreferred_rws)
    CI_rws(kk,:) = ci(perf_unpreferred_rws(:,kk));
end

for kk = 1:length(perf_unpreferred_gbv)
    CI_gbv(kk,:) = ci(perf_unpreferred_gbv(:,kk));
end

for kk = 1:length(perf_unpreferred_gbs)
    CI_gbs(kk,:) = ci(perf_unpreferred_gbs(:,kk));
end

for kk = 1:length(perf_unpreferred_ggs)
    CI_ggs(kk,:) = ci(perf_unpreferred_ggs(:,kk));
end
for kk = 1:length(perf_unpreferred_gbf)
    CI_gbf(kk,:) = ci(perf_unpreferred_gbf(:,kk));
end

figure; hold on;
plot(time, mean(perf_unpreferred_rws))
plot(time, CI_rws(:,1), 'k:')
hold on
plot(time, CI_rws(:,2), 'k:')
title('unpreferred rws')

figure; hold on;
plot(time, mean(perf_unpreferred_ggs))
plot(time, CI_ggs(:,1), 'k:')
hold on
plot(time, CI_ggs(:,2), 'k:')
title('unpreferred ggs')

figure; hold on;
plot(time, mean(perf_unpreferred_gbs))
plot(time, CI_gbs(:,1), 'k:')
hold on
plot(time, CI_gbs(:,2), 'k:')
title('unpreferred gbs')

figure; hold on;
plot(time, mean(perf_unpreferred_gbv))
plot(time, CI_gbv(:,1), 'k:')
hold on
plot(time, CI_gbv(:,2), 'k:')
title('unpreferred gbv')

figure; hold on;
plot(time, mean(perf_unpreferred_gbf))
plot(time, CI_gbf(:,1), 'k:')
hold on
plot(time, CI_gbf(:,2), 'k:')
title('unpreferred gbf')


end