function stats(time, grouppath, balldesigns)

% channels raussuchen (ohne Kleinhirn, Thalamus, Basalganglien)
load ('U:\My Documents\MATLAB\atlas_clusterstatistic_ohne_cerebellum_thalamus_basalganglien.mat')
% atlas_cluster = sortrows(atlas_cluster, 3);
% for k=1:length(atlas_cluster)
%     atlas_cluster(k,4)={k}; 
% end

ind_cluster = cell2mat(atlas_cluster(:,3));
data_like = struct();
data_dislike = struct();
for k=1:length(balldesigns)
  
   data.(balldesigns{k}) = load ([grouppath balldesigns{k} '_ball\session_' balldesigns{k} '.mat']);   
   data.(balldesigns{k}).session.data=data.(balldesigns{k}).session.data(:,ind_cluster,:);
     
   [data.(balldesigns{k}).session, data_like, data_dislike] = regroup_data(data.(balldesigns{k}).session, data_like, data_dislike); 

end   
   
trials_like = [];
trials_dislike = [];
for p = 1:size(data_like.trial, 1)   
    trials_like{p} = squeeze(data_like.trial(p,:,:));
    trials_dislike{p} = squeeze(data_dislike.trial(p,:,:));   
end

data_like.trial = trials_like;
data_dislike.trial = trials_dislike;
clearvars trials_like trials_dislike


for k=1:length(atlas_cluster)
    label{k} = [num2str(atlas_cluster{k,3}) '_' atlas_cluster{k,2}]; 
end
data_like.label = label';
% data_like.dimord = 'rpt_chan_time';
data_like.time = repmat(data.(balldesigns{k}).time(1), 1, size(data_like.trial,1));

data_dislike.label = label';
% data_dislike.dimord = 'rpt_chan_time';
data_dislike.time = repmat(data.(balldesigns{k}).time(1), 1, size(data_dislike.trial,1));


%%
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
cfg.numrandomization = 100;      % number of draws from the permutation distribution
cfg.minnbchan = 2;  % minimum number of neighborhood channels that is required for a selected sample to be included
% cfg.parameter = 'individual'; % oder cfg.parameter = 'trial';
% ntrials_cond = histc(sessions.labels, unique(sessions.labels));
ntrials_cond = [length(data_like.trial) length(data_dislike.trial)];
design  = zeros(2,min(ntrials_cond)+min(ntrials_cond));
design(1,1:min(ntrials_cond)) = 1;
design(1,min(ntrials_cond)+1:min(ntrials_cond)+min(ntrials_cond)) = 2;

% design(2,1:ntrials_cond(1)) = [1:ntrials_cond(1)];
% design(2,ntrials_cond(1)+1:ntrials_cond(1)+ntrials_cond(2)) = [1:ntrials_cond(2)];

design(2,:)=1:min(ntrials_cond)+min(ntrials_cond);
cfg.design = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
load ('E:\adidas\fieldtrip_Auswertung\single_subjects\nl_adi_28\MEG\sourcespace\source_avg\run1\source_avg_appended_conditions_bp1_45Hz.mat', 'source_avg')
chanpos = source_avg.pos(source_avg.inside==1,:);
chanpos = chanpos(ind_cluster,:);

% [neighbours] = prepare_neighbours (chanpos, label, 'roi'); % Method = roi or distance
load('E:\adidas\fieldtrip_Auswertung\group_analysis\source_space\MEG\new_data\neighbours_rois.mat')
cfg.neighbours = neighbours;

[stat]=ft_timelockstatistics(cfg, data_like, data_dislike);

save([path 'stats'], 'stat')


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
    stat.cfg.alpha = 0.025;

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

    roi_num = cell2mat(atlas_cluster(ind,1));
    if 1 == length(unique(roi_num))
        temp = find(cell2mat(atlas_cluster(:,1))== unique(roi_num));
        roi_name = atlas_cluster(temp(1),2);
        scatter3(chanpos(ind,1),chanpos(ind,2),chanpos(ind,3), [], raweffectLike_vs_Dislike.avg(ind, ind_sum_pos_neg(k)), 'filled')
        h_legend = legend({'brain', char(roi_name)}, 'Location','Northeast')  ;  
        
    else 
        for i=1:length(unique(roi_num))
            unique_roi = unique(roi_num);
            temp = find(cell2mat(atlas_cluster(:,1))== unique_roi(i));
            roi_name(i) = atlas_cluster(temp(1),2);
            clearvars temp
        end

        scatter3(chanpos(ind,1),chanpos(ind,2),chanpos(ind,3), [], raweffectLike_vs_Dislike.avg(ind, ind_sum_pos_neg(k)), 'filled')
        h_legend = legend({'brain', char(roi_name)}, 'Location','Northeast') ;   
      
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
        export_fig ([ grouppath balldesign '_sign_clusters.pdf']);
    else
        export_fig ([grouppath balldesign '_sign_clusters.pdf'], '-append')
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

function [data, data_like, data_dislike]=adi_regroup_data(data, data_like, data_dislike)

ind_like=find(data.(balldesigns{k}).labels==1);
ind_dislike=find(data.(balldesigns{k}).labels==2);

like = data.(balldesigns{k}).data(ind_like,:,:);
dislike = data.(balldesigns{k}).data(ind_dislike,:,:);


ntrials_cond = histc(data.(balldesigns{k}).labels, unique(data.(balldesigns{k}).labels));
smaller_cond.label = find(ntrials_cond==min(ntrials_cond));
switch smaller_cond.label
    case 1
       smaller_cond.name = 'like';
    case 2
       smaller_cond.name = 'dislike';
end

cond_tall.label = find(ntrials_cond==max(ntrials_cond));
ind_trial_perm = randperm(max(ntrials_cond), min(ntrials_cond));


switch cond_tall.label
    case 1
       cond_tall.name = 'like';
       like = like(ind_trial_perm,:,:);
    case 2
       cond_tall.name = 'dislike';
       dislike = dislike(ind_trial_perm,:,:);
end   

% data.(balldesigns{k}).data = cat(1, like, dislike);
% data.(balldesigns{k}).labels = [ones(1, size(like,1)), ones(1, size(dislike,1))*2];

switch k
    
    case 1
        data_like.trial = like;
        data_dislike.trial = dislike;  
    otherwise
        data_like.trial = cat(1, data_like.trial, like);
        data_dislike.trial = cat(1, data_dislike.trial, dislike);
end

clearvars like dislike




end