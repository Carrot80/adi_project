function stats(session, time, path, balldesign)

% channels raussuchen (ohne Kleinhirn, Thalamus, Basalganglien)
load ('U:\My Documents\MATLAB\atlas_clusterstatistic_ohne_cerebellum_thalamus_basalganglien.mat')
% atlas_cluster = sortrows(atlas_cluster, 3);
% for k=1:length(atlas_cluster)
%     atlas_cluster(k,4)={k}; 
% end

ind_cluster = cell2mat(atlas_cluster(:,3));

ind_like = find(session.labels==1);
like = [];
like.trial=session.data(ind_like, ind_cluster,:);
% like.time=session.time(ind_like);
like.time=session.time{1,1};

for k=1:length(atlas_cluster)
    label{k} = [num2str(atlas_cluster{k,3}) '_' atlas_cluster{k,2}]; 
end
like.label = label';
like.dimord = 'rpt_chan_time';
like.fsample = median(1./diff(session.time{1,1}));

cfg=[];
cfg.keeptrials = 'yes';
like_avg = ft_timelockanalysis(cfg, like);


% zeitintervall anpassen, evtl. an rms anpassen?

ind_dislike = find(session.labels==2);
dislike = [];
dislike.trial = session.data(ind_dislike, ind_cluster,:);
% dislike.time=session.time(ind_dislike);
dislike.time=session.time{1,1};
dislike.label = label';
dislike.dimord = 'rpt_chan_time';
dislike.fsample = median(1./diff(session.time{1,1}));

cfg=[];
cfg.keeptrials = 'yes';
dislike_avg = ft_timelockanalysis(cfg, dislike);

ntrials_cond = histc(session.labels, unique(session.labels));
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
       like.trial = like.trial(ind_trial_perm,:,:);
       like.time = like.time(ind_trial_perm);
    case 2
       cond_tall.name = 'dislike';
       dislike.trial = dislike.trial(ind_trial_perm,:,:);
       dislike.time = dislike.time;
end       

%%
cfg = [];
cfg.latency = time;
cfg.method           = 'montecarlo';    % use the Monte Carlo Method to calculate the significance probability
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterthreshold = 'nonparametric_common';
cfg.clusteralpha     = 0.05;     % alpha level of the sample-specific test statistic that will be used for thresholding                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum';
cfg.clustercritval = 0.05;                               % in the clustering algorithm (default=0).
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;   % alpha level of the permutation test
cfg.numrandomization = 100;      % number of draws from the permutation distribution, 800 minimum
cfg.minnbchan = 2;  % minimum number of neighborhood channels that is required for a selected sample to be included
% cfg.parameter = 'individual'; % oder cfg.parameter = 'trial';
% ntrials_cond = histc(sessions.labels, unique(sessions.labels));
ntrials_cond = [size(like.trial,1) size(dislike.trial,1)];

%% hier weitermachen
design = zeros(2,2*ntrials_cond(1));
design(1,:)=session.subject_num(1:2*ntrials_cond(1));

for i = 1:subj
design(1,i) = i;
end
for i = 1:subj
design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;



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

[stat]=ft_timelockstatistics(cfg, like, dislike);

save([path 'stats'], 'stat')
% 
% design = zeros(1,size(like.trial,1)) + size(dislike.trial,1);
% design(1,1:size(like.trial,1)) = 1;
% design(1,size(like.trial,1)+1:(size(like.trial,1)+size(dislike.trial,1))) = 2;
% cfg.design = session.labels;
% stat = ft_sourcestatistics(cfg,sourcepstS1,sourcepreS1);
% stat.pos=template_grid.pos;% keep positions for plotting later

%%

% Then take the difference of the averages using ft_math
% avg_dislike = grandavg_dislike;
% avg_dislike = rmfield(avg_dislike, 'individual');
% avg_dislike.label = grandavg_dislike.label;
% avg_like = grandavg_like;
% avg_like = rmfield(avg_like, 'individual');
% 

dislike_avg.avg = squeeze(mean(dislike_avg.trial));
like_avg.avg = squeeze(mean(like_avg.trial));

dislike_avg.dimord = 'chan_time';
like_avg.dimord = 'chan_time';

time1=nearest(like_avg.time, time(1));
time2=nearest(like_avg.time, time(2));

like_avg_selected = like_avg;
like_avg_selected = rmfield(like_avg_selected, 'trial');
like_avg_selected.avg=like_avg.avg(:,time1:time2);
like_avg_selected.time=like_avg.time(1,time1:time2);

dislike_avg_selected = dislike_avg;
dislike_avg_selected = rmfield(dislike_avg_selected, 'trial');
dislike_avg_selected.avg=dislike_avg.avg(:,time1:time2);
dislike_avg_selected.time=dislike_avg.time(1,time1:time2);

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
        export_fig ([ path balldesign '_sign_clusters.pdf']);
    else
        export_fig ([path balldesign '_sign_clusters.pdf'], '-append')
    end
end

savefig(h, [path balldesign '_sign_clusters.fig'])
close(h)



myVideo = VideoWriter([balldesign '_movie' '.avi']);
myVideo.FrameRate = 2;  % Default 30
open(myVideo);
cd(path)
for i=1:length(M)
    writeVideo(myVideo, M(i));
end
close(myVideo);
save([path 'fig_cdata'], 'M')
fig = figure;
movie(fig,M,5,3)

% 
% import mlreportgen.ppt.*
% slides = Presentation([balldesign '_sign_cluster']);
% add(slides,'Title Slide');
% contents = find(slides,'Title');
% replace(contents(1),'My First Presentation');
% close(slides);


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