function stats(time, grouppath, balldesigns, balldesigns_names, roi, condition)

% channels raussuchen (ohne Kleinhirn, Thalamus, Basalganglien)
load ('U:\My Documents\MATLAB\atlas_clusterstatistic_ohne_cerebellum_thalamus_basalganglien.mat')
% atlas_cluster = sortrows(atlas_cluster, 3);
% for k=1:length(atlas_cluster)
%     atlas_cluster(k,4)={k}; 
% end

cond.(condition)=struct();
ind_cluster = cell2mat(atlas_cluster(:,3));   
    
    
if ~isempty(roi)
    
    ind=zeros(1,1);
    for p=1:length(roi)
        [indx.(roi{p})] = find(strcmp(atlas_cluster(:,2), roi{p}));
    end

    vox_index = indx.(roi{1});

    for p=2:length(roi)
        vox_index =cat(1,vox_index, indx.(roi{p}));
    end
    for k=1:length(balldesigns)
        data.(balldesigns{k}) = load ([grouppath balldesigns{k} '_ball\session_' balldesigns{k} '.mat']); 
        data.(balldesigns{k}).session.data=data.(balldesigns{k}).session.data(:,vox_index,:);
        [data.(balldesigns{k}).session, cond] = regroup_data(data.(balldesigns{k}).session, cond, k, condition); 
    end
    for k=1:length(vox_index)
        label{k} = [num2str(atlas_cluster{vox_index(k),3}) '_' atlas_cluster{vox_index(k),2}]; 
    end

else
    ind_cluster = cell2mat(atlas_cluster(:,3));

    for k=1:length(balldesigns)
  
        data.(balldesigns{k}) = load ([grouppath balldesigns{k} '_ball\session_' balldesigns{k} '.mat']); 
           if 1==isequal(size(data.(balldesigns{k}).session.data,2), length(ind_cluster))
                data.(balldesigns{k}).session.data=data.(balldesigns{k}).session.data(:,:,:);
           else
                data.(balldesigns{k}).session.data=data.(balldesigns{k}).session.data(:,ind_cluster,:);
           end
        [data.(balldesigns{k}).session, cond] = regroup_data(data.(balldesigns{k}).session, cond, k, condition); 

    end   
    
    for k=1:length(atlas_cluster)
        label{k} = [num2str(atlas_cluster{k,3}) '_' atlas_cluster{k,2}]; 
    end

end

   
trials = [];

for p = 1:size(cond.(condition).trial, 1)   
    trials{p} = squeeze(cond.(condition).trial(p,:,:));
end

data_cond.trial = trials;

k=1;
data_cond.label = label';
% data_like.dimord = 'rpt_chan_time';
data_cond.time = repmat(data.(balldesigns{k}).session.time(1), 1, size(data_cond.trial,2));


for i=1:length(cond.(condition).balldesign)
    switch char(cond.(condition).balldesign{i})
        
        case balldesigns_names{1}
            data_cond.balldesign_num(i) = 1;
        case balldesigns_names{2}
            data_cond.balldesign_num(i) = 2;
        case balldesigns_names{3}
            data_cond.balldesign_num(i) = 3;
        case balldesigns_names{4}
            data_cond.balldesign_num(i) = 4;
    end
   
end

counts=histcounts(data_cond.balldesign_num);
[indx_min_counts ] = find(counts==min(counts));
trls_min = min(counts);

like=[data_cond.trial(randperm(counts(1), trls_min)) data_cond.trial(counts(1)+randperm(counts(2), trls_min)) data_cond.trial(counts(1)+counts(2)+randperm(counts(3), trls_min))];
balldesign_like=[cond.(condition).balldesign(randperm(counts(1), trls_min)) cond.(condition).balldesign(counts(1)+randperm(counts(2), trls_min)) cond.(condition).balldesign(counts(1)+counts(2)+randperm(counts(3), trls_min)) ];
time_red = [data_cond.time(1:trls_min) data_cond.time(counts(1)+1:counts(1)+trls_min) data_cond.time(counts(1)+counts(2)+1:counts(1)+counts(2)+trls_min)];
response_like=[cond.(condition).response_label(1:trls_min) cond.(condition).response_label(counts(1)+1:counts(1)+trls_min) cond.(condition).response_label(counts(1)+counts(2)+1:counts(1)+counts(2)+trls_min)];

data_like.trial=like;
data_like.balldesign=balldesign_like;
data_like.time = time_red;
data_like.response_label = response_like;
data_like.label=label';

data_dislike.trial(661:670)=[];
data_dislike.trial(326:335)=[];

data_dislike.balldesign(661:670)=[];
data_dislike.balldesign(326:335)=[];

data_dislike.time(661:670)=[];
data_dislike.time(326:335)=[];

%%
cfg = [];
cfg.latency = time;
cfg.method           = 'montecarlo';    % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_indepsamplesT'; % independent t-test statt dependent t-test
cfg.correctm         = 'cluster';
cfg.correcttail = 'prob';
cfg.clusterthreshold = 'nonparametric_common';
cfg.clusteralpha     = 0.05;     % alpha level of the sample-specific test statistic that will be used for thresholding                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum';
cfg.clustercritval = 0.05;                               % in the clustering algorithm (default=0).
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;   % alpha level of the permutation test
cfg.numrandomization = 500;      % number of draws from the permutation distribution
cfg.minnbchan = 2;  % minimum number of neighborhood channels that is required for a selected sample to be included
% cfg.parameter = 'individual'; % oder cfg.parameter = 'trial';
% ntrials_cond = histc(sessions.labels, unique(sessions.labels));
ntrials_cond = [length(data_like.trial) length(data_dislike.trial)];
design  = zeros(2,min(ntrials_cond)+min(ntrials_cond));
design(1,1:min(ntrials_cond)) = 1;
design(1,min(ntrials_cond)+1:min(ntrials_cond)+min(ntrials_cond)) = 2;

% design(2,1:ntrials_cond(1)) = [1:ntrials_cond(1)];
% design(2,ntrials_cond(1)+1:ntrials_cond(1)+ntrials_cond(2)) = [1:ntrials_cond(2)];

% design(2,:)=1:min(ntrials_cond)+min(ntrials_cond);
cfg.design = design;
cfg.ivar     = 1;
% cfg.uvar     = 2;
load ('E:\adidas\fieldtrip_Auswertung\single_subjects\nl_adi_28\MEG\sourcespace\source_avg\run1\source_avg_appended_conditions_bp1_45Hz.mat', 'source_avg')
chanpos = source_avg.pos(source_avg.inside==1,:);
chanpos = chanpos(ind_cluster,:);

% [neighbours] = prepare_neighbours (chanpos, label, 'roi'); % Method = roi or distance
load('E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\neighbours_rois.mat')
if 1==exist('vox_index', 'var')
    cfg.neighbours = neighbours(vox_index);
else
    cfg.neighbours = neighbours;
end

[stat]=ft_timelockstatistics(cfg, data_like, data_dislike);

save([grouppath 'stats'], 'stat')


%%

cfg=[];
cfg.keeptrials = 'yes';
data_like_avg = ft_timelockanalysis(cfg, data_like);

cfg=[];
cfg.keeptrials = 'yes';
data_dislike_avg = ft_timelockanalysis(cfg, data_dislike);

dislike_avg.avg = squeeze(mean(data_dislike_avg.trial));
like_avg.avg = squeeze(mean(data_like_avg.trial));



time1=nearest(data_like_avg.time, time(1));
time2=nearest(data_like_avg.time, time(2));

like_avg_selected = like_avg;
like_avg_selected.avg=like_avg.avg(:,time1:time2);
like_avg_selected.time=data_like_avg.time(1,time1:time2);

dislike_avg_selected = dislike_avg;
dislike_avg_selected.avg=dislike_avg.avg(:,time1:time2);
dislike_avg_selected.time=data_dislike_avg.time(1,time1:time2);

dislike_avg_selected.dimord = 'chan_time';
like_avg_selected.dimord = 'chan_time';
dislike_avg_selected.label = label';
like_avg_selected.label = label';

cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
raweffectLike_vs_Dislike = ft_math(cfg,like_avg_selected,dislike_avg_selected); % geht das überhaupt im source space?


% Then, find which clusters are significant, outputting their indices as held in stat.posclusters
% In case you have downloaded and loaded the data, ensure stat.cfg.alpha exist
if ~isfield(stat.cfg,'alpha') 
    stat.cfg.alpha = 0.05;

end % stat.cfg.alpha was moved as the downloaded data was processed by an additional FieldTrip function to anonymize the data.

for i=1:size(stat.posclusters,2)
    poscluster_pval(i)= stat.posclusters(i).prob;
end

pos_signif_clust = find(poscluster_pval < stat.cfg.alpha);
% (stat.cfg.alpha is the alpha level we specified earlier for cluster comparisons; In this case, 0.025)
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust);

% and now for the negative clusters...
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
neg = ismember(stat.negclusterslabelmat, neg_signif_clust);


timestep = 0.01;	% timestep between time windows for each subplot (in seconds) 
sampling_rate = 256;	% Data has a temporal resolution of 300 Hz
sample_count = length(stat.time); % number of temporal samples in the statistics object
% j = [time(1):timestep:time(2)]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
j = [stat.time(1):timestep:stat.time(end)]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot

m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in MEEG samples
m=ceil(m);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order  
[i1,i2] = match_str(raweffectLike_vs_Dislike.label, stat.label);

pos_neg = pos+neg;

ind_sum_pos_neg=find(sum(pos_neg));
sign_time = stat.time(ind_sum_pos_neg);


for k=1:length(sign_time)

    h(k)=figure
    scatter3(chanpos(:,1),chanpos(:,2),chanpos(:,3), 'k')
    hold on
    ind = find(pos_neg(:,ind_sum_pos_neg(k)));

    if ~isempty(roi)
        roi_num = cell2mat(atlas_cluster(vox_index(ind,1)));           
    else
        roi_num = cell2mat(atlas_cluster(ind,1));
    end

    if 1 == length(unique(roi_num))
        temp = find(cell2mat(atlas_cluster(:,1))== unique(roi_num));
        roi_name = atlas_cluster(temp(1),2);
        scatter3(chanpos(ind,1),chanpos(ind,2),chanpos(ind,3), [], raweffectLike_vs_Dislike.avg(ind, ind_sum_pos_neg(k)), 'filled')
        h_legend = legend({'brain', char(roi_name)}, 'Location','Northeast')  ;  
        clearvars roi_name  
    else 
        for i=1:length(unique(roi_num))
            unique_roi = unique(roi_num);
            temp = find(cell2mat(atlas_cluster(:,1))== unique_roi(i));
            roi_name(i) = atlas_cluster(temp(1),2);
            clearvars temp
        end

        scatter3(chanpos(ind,1),chanpos(ind,2),chanpos(ind,3), [], raweffectLike_vs_Dislike.avg(ind, ind_sum_pos_neg(k)), 'filled')
        h_legend = legend({'brain', char(roi_name)}, 'Location','Northeast') ;   
        clearvars roi_name  
    end

    title([balldesign ' ' num2str(sign_time(k)) ' s' ])
%     set(h_legend,'Position',[0 0 0.4429 0.3349])
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    ax.SortMethod ='ChildOrder'; 
    position = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), position(3)+ti(1)+ti(3), position(4)+ti(2)+ti(4)];
    M(k)=getframe(ax,rect);

    if k==1 
        export_fig ([ grouppath filesep balldesign '_sign_clusters.pdf']);
    else
        export_fig ([grouppath filesep balldesign '_sign_clusters.pdf'], '-append')
    end
end


savefig(h, [grouppath balldesign '_sign_clusters.fig'])
close(h)



myVideo = VideoWriter([balldesign '_movie' '.avi']);
myVideo.FrameRate = 2;  % Default 30
open(myVideo);
cd(grouppath)
for i=1:length(M)
    writeVideo(myVideo, M(i));
end
close(myVideo);
save([grouppath 'fig_cdata'], 'M')
fig = figure;
movie(fig,M,5,3)

% 
% import mlreportgen.ppt.*
% slides = Presentation([balldesign '_sign_cluster']);
% add(slides,'Title Slide');
% contents = find(slides,'Title');
% replace(contents(1),'My First Presentation');
% close(slides);

figure
plot(data_like.time{1,1}, rms(like_avg.avg))
hold on
plot(data_like.time{1,1}, rms(dislike_avg.avg))
title('like_vs_dislike');
legend('like', 'dislike')

sign_clusters = {'Angular_R'; 'Temporal_Mid_R'; 'Temporal_Inf_R'; 'Temporal_Sup_R';} % nr. 66, 86, 90, 82
ind=zeros(1,1);
for p=1:length(sign_clusters)
    
    [indx.(sign_clusters{p})] = find(strcmp(atlas_cluster(:,2), sign_clusters{p}));
    
end

vox_index = indx.(sign_clusters{1});

for p=2:length(sign_clusters)
    vox_index =cat(1,vox_index, indx.(sign_clusters{p}));
end

figure
plot(data_like.time{1,1}, rms(like_avg.avg(vox_index,:)))
hold on
plot(data_like.time{1,1}, rms(dislike_avg.avg(vox_index,:)))
title('like_vs_dislike');
legend('like', 'dislike')

figure
plot(data_like.time{1,1}, rms(like_avg.avg(indx.(sign_clusters{4}),:)))
hold on
plot(data_like.time{1,1}, rms(dislike_avg.avg(indx.(sign_clusters{4}),:)))
title(['like_vs_dislike_' sign_clusters{4}]);
legend('like', 'dislike')

figure
plot(data_like.time{1,1}, mean(like_avg.avg))
hold on
plot(data_like.time{1,1}, mean(dislike_avg.avg))
title('like_vs_dislike');
legend('like', 'dislike')


end

function x()

cfg=[];
cfg.latency = [0.57 0.7]; 
for i=1:length(dislike_avg)
    data_selected = ft_selectdata(cfg, dislike_avg{i});
    dislike_avg_selected{i} = data_selected;
    clear data_selected
end

for i=1:length(like_avg)
    data_selected = ft_selectdata(cfg, like_avg{i})
    like_avg_selected{i} = data_selected;
    clear data_selected
end

cfg=[];
cfg.latency = [0.57 0.7]; 
session_selected = ft_selectdata(cfg, session);
session_sel = session;
for i=1:length(session.trial)
    data_selected = session.trial{i}(:, nearest(session.time{1}, 0.57): nearest(session.time{1}, 0.7));
    session_sel.trial{i} = data_selected;
    clear data_selected
end
for i=1:length(session.time)
    data_selected = session.time{i}(1, nearest(session.time{1}, 0.57): nearest(session.time{1}, 0.7));
    session_sel.time{i} = data_selected;
    clear data_selected
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [neighbours] = prepare_neighbours(chanpos, label, method)

switch method
    
    case 'distance'

        % use a smart default for the distance
        % neighbourdist = 40 * ft_scalingfactor('mm', 'cm');
        neighbourdist = 40 * 0.1; % da chanpos in cm
        fprintf('using a distance threshold of %g\n', neighbourdist);

        nsensors = length(label);

        % compute the distance between all sensors
        dist = zeros(nsensors,nsensors);
        for i=1:nsensors
          dist(i,:) = sqrt(sum((chanpos(1:nsensors,:) - repmat(chanpos(i,:), nsensors, 1)).^2,2))';
        end

        % find the neighbouring electrodes based on distance
        % later we have to restrict the neighbouring electrodes to those actually selected in the dataset
        channeighbstructmat = (dist<neighbourdist);

        % electrode istelf is not a neighbour
        channeighbstructmat = (channeighbstructmat .* ~eye(nsensors));

        % construct a structured cell array with all neighbours
        neighbours=struct;
        for i=1:nsensors
          neighbours(i).label       = label{i};
          neighbours(i).neighblabel = label(find(channeighbstructmat(i,:)));
        end

        for k=1:60:length(label)

            figure
            scatter3(chanpos(:,1), chanpos(:, 2),chanpos(:,3), 'filled', 'black')
            hold on
            scatter3(chanpos(str2double(neighbours(k).label),1), chanpos(str2double(neighbours(k).label),2),chanpos(str2double(neighbours(k).label),3), 'filled')
            hold on
            scatter3(chanpos(str2double(neighbours(k).neighblabel),1), chanpos(str2double(neighbours(k).neighblabel),2),chanpos(str2double(neighbours(k).neighblabel),3), 'filled')
            hold on
            grid off

        end

    case 'roi'
        
        nsensors = length(label);
        mat_cluster =cell2mat(atlas_cluster(:,1));
        index_rois = unique(mat_cluster);
        
        
        % construct a structured cell array with all neighbours
        neighbours=struct;
        for i=1:nsensors
              neighbours(i).label       = label{i};
              index_roi = cell2mat(atlas_cluster(i,1));
              all_ind_roi = find(ismember(cell2mat(atlas_cluster(:,1)), index_roi ));
              neighbours(i).neighblabel = label(all_ind_roi);
              neighbours(i).label_num = i;
              neighbours(i).neighblabel_num = all_ind_roi;
              clearvars all_ind_roi       
        end
        
        % electrode istelf is not a neighbour:
        for i=1:nsensors
            for k=1:length(neighbours(i).neighblabel)
                temp(k) = strcmp(neighbours(i).neighblabel{k}, neighbours(i).label) ;
            end
            neighbours(i).neighblabel(find(temp)) = [];
            neighbours(i).neighblabel_num(find(temp),:) = [];
            clearvars temp
        end
        
         for k=1:60:nsensors

            figure
            scatter3(chanpos(:,1), chanpos(:, 2),chanpos(:,3))
            hold on
            scatter3(chanpos(neighbours(k).label_num,1), chanpos(neighbours(k).label_num, 2), chanpos(neighbours(k).label_num,3), 'filled')
            hold on
            scatter3(chanpos(neighbours(k).neighblabel_num,1), chanpos(neighbours(k).neighblabel_num,2),chanpos(neighbours(k).neighblabel_num,3), 'filled')
            hold on
            grid off
            title(neighbours(k).label)

        end
        
    end
end

function [data, cond]=regroup_data(data, cond, k, condition)

switch condition
    case 'like'
        ind_cond=find(data.labels==1);
        trials = data.data(ind_cond,:,:);
        
    case 'dislike'
        ind_cond=find(data.labels==2);
        trials = data.data(ind_cond,:,:);
end

switch k
    
    case 1
        cond.(condition).trial = trials;  
        cond.(condition).response_label = data.response_label(ind_cond);
        cond.(condition).balldesign = data.balldesign(ind_cond);
    otherwise
        cond.(condition).trial = cat(1, cond.(condition).trial, trials);
        cond.(condition).balldesign = [cond.(condition).balldesign, data.balldesign(ind_cond)];
        cond.(condition).response_label = [cond.(condition).response_label data.response_label(ind_cond)];
end

end


function test()

 for i=1:size(virtsens_all_subj_zscore(1).trial{1,1},1)
      label{i}=num2str(i);
 end
 
cfg=[];
for k=1:length(virtsens_all_subj_zscore)
    
    virtsens_all_subj_zscore(k).label = label';
    rws_virtsens_all_avg{k} = ft_timelockanalysis(cfg, virtsens_all_subj_zscore(k))

end

cfg = [];
cfg.channel   = 'all';
cfg.keepindividual = 'yes';
cfg.latency   = 'all';
cfg.parameter = 'avg';

rws_dislike = ft_timelockgrandaverage(cfg, rws_virtsens_all_avg{:});


%------------------------------------
dislike=gbs_dislike;
dislike.individual = cat(1, gbs_dislike.individual, ggs_dislike.individual, rws_dislike.individual)

like=rwf_like;
like.individual = cat(1, rwf_like.individual, ggv_like.individual, rwv_like.individual)
like.label =label';
dislike.label =label';
dislike.individual = dislike.individual(:, ind_cluster, :);
like.individual = like.individual(:, ind_cluster, :);



 for i=1:1213
      label{i}=num2str(i);
 end

 for i=1:size(dislike.individual,1)
     time_s{i}=rws_dislike.time;
 end
 
dislike.time=time_s
like.time=time_s;

trials = [];

for p = 1:size(like.individual, 1)   
    trials_like{p} = squeeze(like.individual(p,:,:));
end
like.trial = trials_like;


for p = 1:size(dislike.individual, 1)   
    trials_dislike{p} = squeeze(dislike.individual(p,:,:));
end
dislike.trial = trials_dislike;

%-----------------------------------
rws_dislike.individual = rws_dislike.individual(:,ind_cluster,:);
rwv_like.individual = rwv_like.individual(:,ind_cluster,:);




%--------------------------------
cfg = [];
cfg.latency = time;
cfg.method           = 'montecarlo';    % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_indepsamplesT'; % independent t-test statt dependent t-test
cfg.correctm         = 'cluster';
cfg.clusterthreshold = 'nonparametric_common';
cfg.clusteralpha     = 0.05;     % alpha level of the sample-specific test statistic that will be used for thresholding                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum';
cfg.clustercritval = 0.05;                               % in the clustering algorithm (default=0).
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;   % alpha level of the permutation test
cfg.numrandomization = 500;      % number of draws from the permutation distribution
cfg.minnbchan = 2;  % minimum number of neighborhood channels that is required for a selected sample to be included
% cfg.parameter = 'individual'; % oder cfg.parameter = 'trial';
% ntrials_cond = histc(sessions.labels, unique(sessions.labels));
ntrials_cond = [size(data_like.trial,2) size(data_dislike.trial,2)];
design  = zeros(2,min(ntrials_cond)+min(ntrials_cond));
design(1,1:min(ntrials_cond)) = 1;
design(1,min(ntrials_cond)+1:min(ntrials_cond)+min(ntrials_cond)) = 2;

% design(2,1:ntrials_cond(1)) = [1:ntrials_cond(1)];
% design(2,ntrials_cond(1)+1:ntrials_cond(1)+ntrials_cond(2)) = [1:ntrials_cond(2)];

% design(2,:)=1:min(ntrials_cond)+min(ntrials_cond);
cfg.design = design;
cfg.ivar     = 1;
% cfg.uvar     = 2;
load ('E:\adidas\fieldtrip_Auswertung\single_subjects\nl_adi_28\MEG\sourcespace\source_avg\run1\source_avg_appended_conditions_bp1_45Hz.mat', 'source_avg')
chanpos = source_avg.pos(source_avg.inside==1,:);
chanpos = chanpos(ind_cluster,:);

% [neighbours] = prepare_neighbours (chanpos, label, 'roi'); % Method = roi or distance
load('E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\neighbours_rois.mat')
if 1==exist('vox_index', 'var')
    cfg.neighbours = neighbours(vox_index);
else
    cfg.neighbours = neighbours;
end

data_like.individual = data_like.individual(:,ind_cluster, :);
data_dislike.individual = data_dislike.individual(:,ind_cluster, :);

rwv_like.individual = rwv_like.individual(:,ind_cluster, :);
rws_dislike.individual = rws_dislike.individual(:,ind_cluster, :);

rwv_like.time=repmat(rwv_like.time, 21)

rws_dislike.time
trials = [];

for p = 1:size(rwv_like.individual, 1)   
    trials_like{p} = squeeze(rwv_like.individual(p,:,:));
end
rwv_like.trial = trials_like;


trials = [];

for p = 1:size(rws_dislike.individual, 1)   
    trials_dislike{p} = squeeze(rws_dislike.individual(p,:,:));
end
rws_dislike.trial = trials_dislike;


 for i=1:1213
      label{i}=num2str(i);
 end

 for i=1:length(rws_dislike.trial)
     time_s{i}=rws_dislike.time;
 end
 
 cat(3, rwv_like.time, rwv_like.time)
 repmat(rwv_like.time,2,1)
 
 rwv_like.label = label';
 rws_dislike.label = label';
 
[stat]=ft_timelockstatistics(cfg, data_like, data_dislike);


avg_like=squeeze(mean(data_like.individual,1));
avg_dislike=squeeze(mean(data_dislike.individual,1));

figure
plot(time_s{1}, mean(avg_like)) 
hold on
plot(time_s{1}, mean(avg_dislike)) 

end