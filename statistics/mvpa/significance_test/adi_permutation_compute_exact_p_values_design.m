function [] = adi_permutation_compute_p_values(subjectdir, balldesign, timerange, path2files)

if length(balldesign) == 1
    
    mkfigures(subjectdir, balldesign, timerange, path2files)
    
else

for ii = 1:length(subjectdir)

load([subjectdir(ii).folder filesep subjectdir(ii).name path2files 'performance_real_data.mat'], 'perf')

for kk = 1:length(balldesign)
    if exist([subjectdir(ii).folder filesep subjectdir(ii).name filesep path2files 'performance_permutation_null_distribution_' balldesign{kk} '.mat'], 'file')
        perf_perm.(balldesign{kk}) = load([subjectdir(ii).folder filesep subjectdir(ii).name filesep ...
            path2files 'performance_permutation_null_distribution_' balldesign{kk} '.mat'], ['perm_' balldesign{kk}])
    end
end
designs = fieldnames(perf_perm);
temp = zeros(1, numel(fieldnames(perf_perm)));
for kk = 1:numel(fieldnames(perf_perm))
    if 1 == strcmp(perf.CV(kk).design, designs{kk}) && 1 == isfield(perf.CV, 'trialnumer_likes_trainingsset') && (perf.CV(kk).trialnumer_likes_trainingsset <10 || perf.CV(kk).trialnumer_dislikes_trainingsset <10)
        temp(kk) = 1; 
    elseif 1 == strcmp(perf.CV(kk).design, designs{kk}) && isfield(perf.CV, 'trialnumber_likes_trainingsset') && (perf.CV(kk).trialnumber_likes_trainingsset <10 || perf.CV(kk).trialnumber_dislikes_trainingsset <10)
        temp(kk) = 1;
    end
end

indx = find(temp);
if ~isempty(indx) && isfield(perf, (perf.CV(indx).design))
    perf = rmfield(perf, (perf.CV(indx).design));
    perf_perm = rmfield(perf_perm, (perf.CV(indx).design));
end
clear indx
designs = fieldnames(perf_perm);
fsample = 256.001;
fsample_sec = 1/fsample;
time = -0.5:fsample_sec:1;

%% collect permutation statistic over subjects
acc_perf_balldesigns = zeros(numel(designs), length(time));
for kk = 1:numel(designs)
    acc_perf_balldesigns(kk,:) = perf.(designs{kk}).lda.accuracy;
end
acc_mean_balldesigns_orig_data = nanmean(acc_perf_balldesigns);

acc_perm_mean_balldesigns =  perf_perm.(designs{1}).(['perm_' designs{1}]).accuracy;
for kk = 2:numel(fieldnames(perf_perm))
   acc_perm_mean_balldesigns = (acc_perm_mean_balldesigns+perf_perm.(designs{kk}).(['perm_' designs{kk}]).accuracy)./2;
end

switch ii
    case 1
        acc_perm_mean_balldesigns_all_subjects(:,:) = acc_perm_mean_balldesigns;
    otherwise
        acc_perm_mean_balldesigns_all_subjects = (acc_perm_mean_balldesigns_all_subjects+acc_perm_mean_balldesigns)./2;
end

acc_mean_balldesigns_orig_data_all_subjects(ii,:) = acc_mean_balldesigns_orig_data;

%% compute statistic (exact pvalue) per subject and give out proportion of subjects with significant p-values: 

 perf.accuracy_sign = zeros(1, size(time,2));
 perm_timerange = nearest(time, timerange(1)):1:nearest(time, timerange(2));
for tt = perm_timerange
    prc_ = prctile(acc_perm_mean_balldesigns(:,tt), 95);
    if acc_mean_balldesigns_orig_data(1,tt) > prc_
       perf.accuracy_sign(tt) = 1;
    else
       perf.accuracy_sign(tt) = 0;
    end
    clear prc_
end

 % 2. calculation of p value with FDR-correction
p = nan(1, numel(perm_timerange));
for tt = perm_timerange
    [indx, acc_taller_orig] = find(acc_perm_mean_balldesigns(:,tt) >  acc_mean_balldesigns_orig_data(tt));
    p(tt-(perm_timerange(1)-1)) = (sum(acc_taller_orig)+1)./(1000+1);  % p=(b+1)/(m+1) => see Dima et al(2018), Human Brain Mapping
end
p_ = nan(1, size(acc_perm_mean_balldesigns,2));
p_(perm_timerange) = p;
sign_p = p_<0.05;

p_all_subjects(ii,:) = p(1,:);

% confidence intervals:
for kk = 1:length(acc_perf_balldesigns)
   CI_orig(kk).lda = ci(acc_perf_balldesigns(:,kk));
end
for kk = 1:length(CI_orig)
    ci_orig_lower(kk) = CI_orig(kk).lda(1);
    ci_orig_upper(kk) = CI_orig(kk).lda(2);
end

for kk = 1:size(acc_perm_mean_balldesigns,2)
   CI_perm(kk).lda = ci(acc_perm_mean_balldesigns(:,kk));
end

for kk = 1:length(CI_perm)
    ci_perm_lower(kk) = CI_perm(kk).lda(1);
    ci_perm_upper(kk) = CI_perm(kk).lda(2);
end
    
figure; hold on;
plot(time, acc_mean_balldesigns_orig_data)
plot(time, ci_orig_lower, 'color', [0.5 0.5 0.5])
plot(time, ci_orig_upper, 'color', [0.5 0.5 0.5])
plot(time, nanmean(acc_perm_mean_balldesigns))
plot(time, sign_p*0.4, 'r*'  )
ylim([0.01 1.1])
title('accuracy uncorrected')
xlabel('time')
ylabel('accuracy')


transparency_color = 0.1;
edge = 'b';
transparency_edge = 0.1;
filled = [ci_orig_upper,fliplr(ci_orig_lower)];
xpoints = [time, fliplr(time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

filled_perm = [ci_perm_lower,fliplr(ci_perm_upper)];
fillhandle = fill(xpoints,filled_perm,'r');%plot the data

transparency_color = 0.5;
transparency_edge = 0.5;
edge = 'r';
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

plot(time, nanmean(acc_perm_mean_balldesigns), 'r')
plot(time, ci_perm_lower, 'color', [0.5 0.5 0.5])
plot(time, ci_perm_upper, 'color', [0.5 0.5 0.5])
    

savefig([subjectdir(ii).folder filesep subjectdir(ii).name filesep path2files 'accuracy_orig_data_uncorrected.fig'])
close 
% FDR-correction:
p(2,:) = perm_timerange;
[sort_p, indx_p] = sort(p(1,:)); % 
rank_p = 1:sum(numel(sort_p));%182;
p_adjusted = nan(1,numel(rank_p));
p_adjusted(end) = sort_p(end);  % größter p-Wert kommt ans Ende 
%     
for tt = fliplr(rank_p(1:end-1))
    temp = sort_p(tt)*(numel(sort_p)/rank_p(tt));
    if temp < p_adjusted(rank_p(tt+1))
        p_adjusted(tt) = temp;
    elseif temp > p_adjusted(rank_p(tt+1))
        p_adjusted(tt) = p_adjusted(tt+1);
    else
        p_adjusted(tt) = temp;
    end
    clear temp
end
    
[~, old_order] = sort(indx_p);
p_adjusted_ = p_adjusted(old_order);

p_adjusted_ = [nan(1, perm_timerange(1)-1) p_adjusted_  nan(1,385-numel(nan(1, perm_timerange(1)-1))-numel(p_adjusted_) )];
indx_p_sign = find(p_adjusted_<0.05);
p_sign = nan(1, 385);
p_sign(indx_p_sign) = 1;

p_adjusted_all_subj(ii, :) = p_adjusted;


 % plot significance as stars:
figure; hold on;
plot(time, acc_mean_balldesigns_orig_data)
plot(time, ci_orig_lower, 'color', [0.5 0.5 0.5])
plot(time, ci_orig_upper, 'color', [0.5 0.5 0.5])
plot(time, nanmean(acc_perm_mean_balldesigns), 'r')
plot(time, ci_perm_lower, 'color', [0.5 0.5 0.5])
plot(time, ci_perm_upper, 'color', [0.5 0.5 0.5])
plot(time, p_sign*0.73, 'r*'  )
transparency_color = 0.1;
edge = 'b';
transparency_edge = 0.1;
filled = [ci_orig_upper,fliplr(ci_orig_lower)];
xpoints = [time, fliplr(time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

ylim([0.1 1])
title('accuracy FDR corrected')
xlabel('time')
ylabel('accuracy')
savefig([subjectdir(ii).folder filesep subjectdir(ii).name filesep path2files '\accuracy_orig_data_fdr_ci.fig'])
close  
clear acc_mean_balldesigns_orig_data acc_perf_balldesigns acc_perm_mean_balldesigns perf perf_perm p_adjusted CI_orig ci_perm_lower ci_perm_upper CI_perm ci_orig_lower ci_orig_upper p
end

figure;
plot(time, nanmean(acc_mean_balldesigns_orig_data_all_subjects(:,:))) 
hold on
plot(time, nanmean(acc_perm_mean_balldesigns_all_subjects))
for kk = 1:length(acc_mean_balldesigns_orig_data_all_subjects)
   CI_orig_all_subjects(kk).lda = ci(acc_mean_balldesigns_orig_data_all_subjects(:,kk));
end
for kk = 1:length(CI_orig_all_subjects)
    ci_orig_all_subj_lower(kk) = CI_orig_all_subjects(kk).lda(1);
    ci_orig_all_subj_upper(kk) = CI_orig_all_subjects(kk).lda(2);
end


transparency_color = 0.1;
edge = 'b';
transparency_edge = 0.1;
filled = [ci_orig_all_subj_upper,fliplr(ci_orig_all_subj_lower)];
xpoints = [time, fliplr(time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color



for kk = 1:size(acc_perm_mean_balldesigns_all_subjects, 2)
   CI_perm_all_subjects(kk).lda = ci(acc_perm_mean_balldesigns_all_subjects(:,kk));
end
for kk = 1:size(acc_perm_mean_balldesigns_all_subjects, 2)
    ci_perm_all_subj_lower(kk) = CI_perm_all_subjects(kk).lda(1);
    ci_perm_all_subj_upper(kk) = CI_perm_all_subjects(kk).lda(2);
end
transparency_color = 0.3;
transparency_edge = 0.3;
filled_perm = [ci_perm_all_subj_lower,fliplr(ci_perm_all_subj_upper)];
fillhandle = fill(xpoints,filled_perm,'r');%plot the data
edge = 'r';
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

% find corrected p_values per subject to figure out the proportion of subjects with significant samples:
p_sign_all_subj = zeros(size(p_adjusted_all_subj, 1), size(p_adjusted_all_subj, 2));
for kk = 1:size(p_adjusted_all_subj, 1)
    indx_p_sign_all_subj = find(p_adjusted_all_subj(kk,:)<0.05);
    p_sign_all_subj(kk,indx_p_sign_all_subj) = 1;
end

perf.accuracy_sign = zeros(1, size(time,2));
 perm_timerange = nearest(time, timerange(1)):1:nearest(time, timerange(2));
for tt = perm_timerange
    prc_ = prctile(acc_perm_mean_balldesigns_all_subjects(:,tt), 95);
    if acc_mean_balldesigns_orig_data_all_subjects(1,tt) > prc_
       perf.accuracy_sign(tt) = 1;
    else
       perf.accuracy_sign(tt) = 0;
    end
    clear prc_
end

% find uncorrected p_values per subject to figure out the proportion of subjects with significant samples:
p_sign_all_subj = zeros(size(p_all_subjects, 1), size(p_all_subjects, 2));
for kk = 1:size(p_all_subjects, 1)
    indx_p_sign_all_subj = find(p_all_subjects(kk,:)<0.05);
    p_sign_all_subj(kk,indx_p_sign_all_subj) = 1;
end
sum(p_sign_all_subj)
temp = nan(1, size(acc_perm_mean_balldesigns_all_subjects,2));
temp(perm_timerange)=sum(p_sign_all_subj);

 % 2. calculation of p value over subjects with FDR-correction
p_all_subj = nan(1, numel(perm_timerange));
for tt = perm_timerange
    [indx_all_subj, acc_taller_orig_all_subj] = find(acc_perm_mean_balldesigns_all_subjects(:,tt) >  nanmean(acc_mean_balldesigns_orig_data_all_subjects(:,tt)));
    p_all_subj(tt-(perm_timerange(1)-1)) = (sum(acc_taller_orig_all_subj)+1)./(1000+1);  % p=(b+1)/(m+1) => see Dima et al(2018), Human Brain Mapping
end
p_all_subj_ = nan(1, size(acc_perm_mean_balldesigns_all_subjects,2));
p_all_subj_(perm_timerange) = p_all_subj;
sign_p_all_subj = p_all_subj_<=0.05;

plot(time, sign_p_all_subj*0.45, 'r*')
ylim([0.38 0.68])
ylabel('accuracy')
xlabel('time')
box off

temp = nan(28, size(acc_perm_mean_balldesigns_all_subjects,2));
temp(:, perm_timerange) = p_all_subjects;
p_comp_1 = temp(:, nearest(time, 0.03): nearest(time, 0.15));

test = zeros(size(p_comp_1,1), size(p_comp_1,2));
for kk = 1:size(p_comp_1,1)
    ind = find(p_comp_1(kk,:)<=0.05);
    test(kk,ind) = 1;
end
sum(test,2)
% ci:
% saved in 'E:\adidas\fieldtrip_Auswertung\group_analysis\sensor_space\significance_testing_mvpa'

end
end


function [] = mkfigures(subjectdir, balldesign, timerange, path2files)

data_balldesign_all_subj = nan(length(subjectdir), 385);
data_perm_balldesign_all_subj = nan(1000, 385);

for ii = 1:length(subjectdir)
   
    load([subjectdir(ii).folder filesep subjectdir(ii).name filesep path2files 'performance_real_data.mat'], 'perf')
    order_balls = [];
    for kk = 1:length(perf.CV)
        order_balls{kk, 1} = perf.CV(kk).design;
    end
    indx = find(strcmp(order_balls, balldesign));
    
    if 1 == ~isempty(indx) && 1 == isfield(perf.CV, 'trialnumer_likes_trainingsset') && (perf.CV(indx).trialnumer_likes_trainingsset <10 || perf.CV(indx).trialnumer_dislikes_trainingsset <10);
        data_balldesign_all_subj(ii,:) = NaN;
    elseif 1 == ~isempty(indx) && isfield(perf.CV, 'trialnumber_likes_trainingsset') && (perf.CV(indx).trialnumber_likes_trainingsset <10 || perf.CV(indx).trialnumber_dislikes_trainingsset <10);
         data_balldesign_all_subj(ii,:) = NaN;
    elseif 1 == isempty(indx)
         data_balldesign_all_subj(ii,:) = NaN;
    else
        data_balldesign_all_subj(ii,:) = perf.(balldesign{1}).lda.accuracy;
    end
    if 2 == exist([subjectdir(ii).folder filesep subjectdir(ii).name path2files 'performance_permutation_null_distribution_' balldesign{1} '.mat'], 'file')
        perf_perm.(balldesign{1}) = load([subjectdir(ii).folder filesep subjectdir(ii).name filesep ...
            path2files 'performance_permutation_null_distribution_' balldesign{1} '.mat'], ['perm_' balldesign{1}]);
    end
    if 1 == exist('perf_perm', 'var')
        switch ii
            case 1
                data_perm_balldesign_all_subj =  perf_perm.(balldesign{1}).(['perm_' balldesign{1}]).accuracy;   
            otherwise
                data_perm_balldesign_all_subj = (data_perm_balldesign_all_subj+perf_perm.(balldesign{1}).(['perm_' balldesign{1}]).accuracy)./2;
        end  
    end
end
fsample = 256.001;
fsample_sec = 1/fsample;
time = -0.5:fsample_sec:1;

data_balldesign_all_subj(find(isnan(data_balldesign_all_subj(:,1))),:) = [];

% confidence intervals:
for kk = 1:length(data_balldesign_all_subj)
   CI_orig(kk).lda = ci(data_balldesign_all_subj(:,kk));
end
for kk = 1:length(data_balldesign_all_subj)
    ci_orig_lower(kk) = CI_orig(kk).lda(1);
    ci_orig_upper(kk) = CI_orig(kk).lda(2);
end

for kk = 1:size(data_perm_balldesign_all_subj,2)
   CI_perm(kk).lda = ci(data_perm_balldesign_all_subj(:,kk));
end

for kk = 1:length(CI_perm)
    ci_perm_lower(kk) = CI_perm(kk).lda(1);
    ci_perm_upper(kk) = CI_perm(kk).lda(2);
end

figure; hold on;
plot(time, nanmean(data_balldesign_all_subj))
plot(time, nanmean(data_perm_balldesign_all_subj))

transparency_color = 0.1;
edge = 'b';
transparency_edge = 0.1;
filled = [ci_orig_upper,fliplr(ci_orig_lower)];
xpoints = [time, fliplr(time)];
fillhandle = fill(xpoints,filled,'b');%plot the data
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

filled_perm = [ci_perm_lower,fliplr(ci_perm_upper)];
fillhandle = fill(xpoints,filled_perm,'r');%plot the data

transparency_color = 0.5;
transparency_edge = 0.5;
edge = 'r';
set(fillhandle,'EdgeColor',edge,'FaceAlpha',transparency_edge,'EdgeAlpha',transparency_color);%set edge color

perm_timerange = nearest(time, timerange(1)):1:nearest(time, timerange(2));

p_all_subj = nan(1, numel(perm_timerange));
for tt = perm_timerange
    [indx_all_subj, acc_taller_orig_all_subj] = find(data_perm_balldesign_all_subj(:,tt) >  nanmean(data_balldesign_all_subj(:,tt)));
    p_all_subj(tt-(perm_timerange(1)-1)) = (sum(acc_taller_orig_all_subj)+1)./(1000+1);  % p=(b+1)/(m+1) => see Dima et al(2018), Human Brain Mapping
end
p_all_subj_ = nan(1, size(data_perm_balldesign_all_subj,2));
p_all_subj_(perm_timerange) = p_all_subj;
sign_p_all_subj = p_all_subj_<=0.05;

plot(time, sign_p_all_subj*0.45, 'r*')
ylim([0.3 0.85])
ylabel('accuracy')
xlabel('time')
box off
set(gcf, 'Position', [434 664 762 303])
gcf

% save
 
 end