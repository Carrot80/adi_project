
function adi_source_ROI (path2vol, path2data, mriPath, outPath_extdisc, freqbandname, like, dislike, dontcare)

% calculations spatial filter per Run and Subject, 

%% run 1:
if ~exist([outPath_extdisc 'MEG\sourcespace\run1\vs_mean_all_rois_' dislike '_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath,  num2str(1), freqbandname, like, dislike, dontcare, outPath_extdisc)                                                                   
end

%% run2:

if ~exist([outPath_extdisc 'MEG\sourcespace\run2\vs_mean_all_rois_' dislike '_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(2), freqbandname, like, dislike, dontcare, outPath_extdisc)
end

%% run3:

if ~exist([outPath_extdisc 'MEG\sourcespace\run3\vs_mean_all_rois_' dislike '_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(3), freqbandname, like, dislike, dontcare, outPath_extdisc)

end

end

function [virtsens] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, run, freq, like, dislike, dontcare, outPath_extdisc)
    
% load data:
if exist([path2data, like '500_' run '.mat'], 'file')
    data_like = load ([path2data, like '500_' run '.mat'], 'cleanMEG_interp');
    [data_like_bpfreq] = adi_bpfilter(data_like.cleanMEG_interp, freq);
    clear data_like
end

if exist([path2data, dislike '500_' run '.mat'], 'file')
    data_dislike = load ([path2data, dislike '500_' run '.mat'], 'cleanMEG_interp');
    [data_dislike_bpfreq] = adi_bpfilter(data_dislike.cleanMEG_interp, freq);
    clear data_dislike
end

if exist([path2data, dontcare '500_' run '.mat'], 'file')
    data_dontcare = load ([path2data, dontcare '500_' run '.mat'], 'cleanMEG_interp');
    [data_dontcare_bpfreq] = adi_bpfilter(data_dontcare.cleanMEG_interp, freq);
    clear data_dontcare
end

if exist('data_dontcare_bpfreq', 'var') && exist('data_like_bpfreq', 'var') && exist('data_dislike_bpfreq', 'var') 
    cfg = [];
    data_appended = ft_appenddata(cfg, data_like_bpfreq, data_dislike_bpfreq, data_dontcare_bpfreq);
elseif exist('data_dontcare_bpfreq', 'var') && exist('data_dislike_bpfreq', 'var')     
    cfg = [];
    data_appended = ft_appenddata(cfg, data_dontcare_bpfreq, data_dislike_bpfreq);
elseif exist('data_like_bpfreq', 'var') && exist('data_dislike_bpfreq', 'var') 
    cfg = [];
    data_appended = ft_appenddata(cfg, data_like_bpfreq, data_dislike_bpfreq);
elseif exist('data_dontcare_bpfreq', 'var')  && exist('data_like_bpfreq', 'var')
    cfg = [];
    data_appended = ft_appenddata(cfg, data_dontcare_bpfreq, data_like_bpfreq);
end
        
   
%% Computing source modell: 
% Following the construction of the volume conduction model (singleshell), we need to discretize the brain into a source model or grid.
% For each grid point in the brain, the lead field matrix is calculated later When constructing the source model, you might want to keep in mind 
% that averaging and statistics over subjects can only be done if the individual subjects source reconstructed results are mapped onto a common space. 
% Now, we will construct a regular grid in MNI template space and spatially deform this grid to the individual subjects brain. The following code loads the template grid that is included in the FieldTrip release:

% load templates:

% load template grid:
template_vol = load('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template/headmodel/standard_singleshell');

%    template_sourcemodel = load('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\sourcemodel\standard_sourcemodel3d8mm');
load ([mriPath, 'mri_realigned.mat']);
mri_realigned_cm = ft_convert_units(mri_realigned, 'cm');
% inverse-warp the template grid to subject specific coordinates 
%oder:
cfg = [];
cfg.grid.xgrid  = -20:1:20;
cfg.grid.ygrid  = -20:1:20;
cfg.grid.zgrid  = -20:1:20;
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.headmodel  = template_vol.vol;
template_grid  = ft_prepare_sourcemodel(cfg); % = sourcemodel
     
atlas = ft_read_atlas('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\atlas\aal\ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas,'cm'); % assure that atlas and template_grid are expressed in the %same unitscfg = [];
cfg = [];
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel;
cfg.inputcoord = 'mni';
mask = ft_volumelookup(cfg,template_grid);
    
% create temporary mask according to the atlas entries
tmp                  = repmat(template_grid.inside,1,1);
tmp(tmp==1)          = 0;
tmp(mask)            = 1;

% define inside locations according to the atlas based mask
template_grid.inside = tmp;

% plot the atlas based grid
figure;ft_plot_mesh(template_grid.pos(template_grid.inside,:));

cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri_realigned_cm;
cfg.grid.unit      ='cm';
sourcemodel_inversewarped = ft_prepare_sourcemodel(cfg); % sourcemodel = grid
figure;ft_plot_mesh(sourcemodel_inversewarped.pos(sourcemodel_inversewarped.inside,:));
    
% Finally, it is wise to check whether all computed objects align well with one another, i.e. whether the grid is correctly 
% placed within the volume conduction model, which both have to be aligned with the MEG sensors. 
% Note that all objects that we plot need to be expressed in the same units and the same coordinate space. 
% Here, we need to transform the head model from 'mm' into 'cm'.
hdm_ind = load([path2vol 'vol.mat']);
hdm_ind_cm = ft_convert_units(hdm_ind.vol, 'cm');

figure; hold on     % plot all objects in one figure
ft_plot_vol(hdm_ind_cm, 'edgecolor', 'none', 'facealpha', 0.4)
alpha 0.4           % make the surface transparent
ft_plot_mesh(sourcemodel_inversewarped.pos(sourcemodel_inversewarped.inside,:));
%      ft_plot_mesh(sourcemodel_inversewarped.pos(:,:));
hold on
data_bpfreq.grad = ft_convert_units(data_appended.grad, 'cm');
hold on
ft_plot_sens(data_bpfreq.grad);
pathFig = [outPath_extdisc 'MEG\sourcespace\run' run]; 
if ~exist(pathFig, 'dir')
    mkdir(pathFig)
end
savefig([pathFig '\sensor_position_check.fig'])
close all

cfg             = [];
cfg.grid        = sourcemodel_inversewarped;
cfg.headmodel   = hdm_ind_cm;
cfg.vol   = hdm_ind_cm;
cfg.channel     = data_appended.label;
cfg.grad        = data_appended.grad;
cfg.normalize = 'yes'; % If you are not contrasting the activity of interest against another condition or baseline time-window, then you may choose to normalize the lead field in this step (cfg.normalize='yes'), which will help control against the power bias towards the center of the head.
cfg.reducerank  = '2'; % or number (default = 3 for EEG, 2 for MEG
sourcemodel_lf  = ft_prepare_leadfield(cfg, data_appended); % sourcemodel = grid
    
cfg = [];
cfg.covariance='yes';
cfg.covariancewindow = [0 1];
cfg.keeptrials = 'yes';
cfg.vartrllength = 2;
avg_data_appended = ft_timelockanalysis(cfg, data_appended);
% avg_data_like = ft_timelockanalysis(cfg, data_like_bpfreq);
% avg_data_dislike = ft_timelockanalysis(cfg, data_dislike_bpfreq);
% if exist ('data_dontcare_bpfreq', 'var')
%     avg_data_dontcare = ft_timelockanalysis(cfg, data_dontcare_bpfreq);
% end

cfg=[];
cfg.method = 'lcmv';
cfg.grid = sourcemodel_lf;
cfg.vol = hdm_ind_cm;
cfg.lcmv.keepfilter = 'yes';
cfg.channel = data_appended.label;
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.fixedori = 'no';
cfg.lcmv.lamda = '5%';
source_avg = ft_sourceanalysis(cfg, avg_data_appended);
source_avg.pos = template_grid.pos;
% hier weitermachen: common filter multiplizeren mit sensor data
% condition like, dislike, dontcare ;
% ROI-Analyse: vorher filter aus ROI extrahieren und mitteln; dann
% euklidsche norm berechnen
% auch source statistics durchführen, einmal korrected und einmal
% uncorrected

% exctract spatialfilter from ROIs:
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.interpmethod = 'nearest';
cfg.parameter    = 'tissue';
atlas_downsampled2source_avg  = ft_sourceinterpolate(cfg, atlas, source_avg);
atlas_downsampled2source_avg.coordsys = 'mni';
atlas_downsampled2source_avg.tissuelabel = atlas_downsampled2source_avg.tissuelabel(1,1:90);
cfg=[];
cfg.inputcoord = 'mni';
cfg.atlas = atlas_downsampled2source_avg;
for k = 1:length(atlas_downsampled2source_avg.tissuelabel)
    cfg.roi = atlas_downsampled2source_avg.tissuelabel{1,k};
%         mask_roi_pos_source = ft_volumelookup(cfg, source_avg); % template_grid oder source_avg? scheint egal zu sein
    mask_roi_pos_template_grid = ft_volumelookup(cfg, template_grid);
%         mast_roi_ind = find(mask_roi_pos_source);
    mast_roi_ind2 = find(mask_roi_pos_template_grid); % unterschiedliche positionen
%         filter_ROI = source_avg.avg.filter(mast_roi_ind); % pos überprüfen!
    filter_ROI2 = source_avg.avg.filter(mast_roi_ind2); %
%         isequal(filter_ROI, filter_ROI2);
%         isequal(mast_roi_ind, mast_roi_ind2);
    
%     figure
%     plot3(template_grid.pos(find(template_grid.inside),1),template_grid.pos(find(template_grid.inside),2),template_grid.pos(find(template_grid.inside),3), '*')
%     hold on
%     plot3(template_grid.pos(mast_roi_ind2,1),template_grid.pos(mast_roi_ind2,2),template_grid.pos(mast_roi_ind2,3), '*')
%     savefig([pathFig filesep 'ROI_' atlas_downsampled2source_avg.tissuelabel{1,k} '.fig'])
%     close

    % condition like
    if exist('data_like_bpfreq', 'var')
        virtsens_roi_like = {};

        for i = 1:length(data_like_bpfreq.trial)
            virtsens_roi_like.trial{i} = filter_ROI2*data_like_bpfreq.trial{i};
        end
%         % noise abziehen: Mittelwert von kompletten Filter abziehen?
%         ind_filter = find(source_avg.inside);
%         spatialfilter_inside = source_avg.avg.filter(ind_filter,1);
% %         spatialfilter = cat(1,spatialfilter_inside);
%         spatialfilter_inside = cat(1,source_avg.avg.filter{:})
%        
%         
%         filter_ROI2_row1 = zeros(length(filter_ROI2), 248);
%         filter_ROI2_row2 = zeros(length(filter_ROI2), 248);
%         filter_ROI2_row2 = zeros(length(filter_ROI2), 248);
%         
%         for i = 1:length(filter_ROI2)
%          filter_ROI2_row1(i,:) = filter_ROI2{i,1}(1,:);
%          filter_ROI2_row2(i,:) = filter_ROI2{i,1}(2,:);
%          filter_ROI2_row3(i,:) = filter_ROI2{i,1}(3,:);
%         end
%         
%         weights_noise(1,:) = mean(abs(filter_ROI2_row1),2);
%         weights_noise(2,:) = mean(abs(filter_ROI2_row2),2);
%         weights_noise(3,:) = mean(abs(filter_ROI2_row3),2);
%        
%         for i = 1:length(virtsens_roi_like.trial)
%             virtsens_ns_roi_like.trial{i} = virtsens_roi_like.trial{1,i}./repmat(weights_noise,1,size(virtsens_roi_like.trial{1,i},2)); % => all will have a similar noise level
%         end
% 
%         virtsens_ns_roi_like.time = data_like_bpfreq.time;
%         virtsens_ns_roi_like.fsample = data_like_bpfreq.fsample;
%         virtsens_ns_roi_like.tissue = atlas_downsampled2source_avg.tissuelabel{1,k};
%         virtsens_ns_roi_like_all.(atlas_downsampled2source_avg.tissuelabel{1,k}) = virtsens_ns_roi_like;     

        for p = 1:length(virtsens_roi_like.trial)
%             euclid_roi = zeros(length(virtsens_roi_like.trial{1,p}), length(virtsens_roi_like.trial{1,p}{1,1}));
            for i = 1:length(virtsens_roi_like.trial{1,p})
                Euklidsche_dist = sqrt(virtsens_roi_like.trial{1,p}{i,:}(1,:).^2+virtsens_roi_like.trial{1,p}{i,:}(2,:).^2+virtsens_roi_like.trial{1,p}{i,:}(3,:).^2);
                euclid_roi(i,:) = Euklidsche_dist;
                clear Euklidsche_dist
            end
            vs_euclid_roi_like.trial{1,p} = euclid_roi;
            clear euclid_roi
        end

            vs_euclid_roi_like_all.(atlas_downsampled2source_avg.tissuelabel{1,k}) = vs_euclid_roi_like;
            vs_euclid_roi_like_all.(atlas_downsampled2source_avg.tissuelabel{1,k}).time = data_appended.time{1,1};
            clear vs_euclid_roi_like
    end
        
 %% condition dislike
        
    if exist('data_dislike_bpfreq', 'var')
        virtsens_roi_dislike = {};

        for i = 1:length(data_dislike_bpfreq.trial)
            virtsens_roi_dislike.trial{i} = filter_ROI2*data_dislike_bpfreq.trial{i};
        end

        for p = 1:length(virtsens_roi_dislike.trial)
%             euclid_roi = zeros(length(virtsens_roi_dislike.trial{1,p}), length(virtsens_roi_dislike.trial{1,p}{1,1}));
            for i = 1:length(virtsens_roi_dislike.trial{1,p})
                Euklidsche_dist = sqrt(virtsens_roi_dislike.trial{1,p}{i,:}(1,:).^2+virtsens_roi_dislike.trial{1,p}{i,:}(2,:).^2+virtsens_roi_dislike.trial{1,p}{i,:}(3,:).^2);
                euclid_roi(i,:) = Euklidsche_dist;
                clear Euklidsche_dist
            end
            vs_euclid_roi_dislike.trial{1,p} = euclid_roi;
            clear euclid_roi
        end

        vs_euclid_roi_dislike_all.(atlas_downsampled2source_avg.tissuelabel{1,k}) = vs_euclid_roi_dislike;
        vs_euclid_roi_dislike_all.(atlas_downsampled2source_avg.tissuelabel{1,k}).time = data_appended.time{1,1};
        clear vs_euclid_roi_dislike    
    end

 %% condition dontcare
    if exist('data_dontcare_bpfreq', 'var')
       virtsens_roi_dontcare = {};

        for i = 1:length(data_dontcare_bpfreq.trial)
            virtsens_roi_dontcare.trial{i} = filter_ROI2*data_dontcare_bpfreq.trial{i};
        end

%             vs_euclid_roi_dontcare = zeros(length(virtsens_roi_dontcare.trial), 248);
        for p = 1:length(virtsens_roi_dontcare.trial)
            euclid_roi = zeros(length(virtsens_roi_dontcare.trial{1,p}), length(virtsens_roi_dontcare.trial{1,p}{1,1}));
            for i = 1:length(virtsens_roi_dontcare.trial{1,p})
                Euklidsche_dist = sqrt(virtsens_roi_dontcare.trial{1,p}{i,:}(1,:).^2+virtsens_roi_dontcare.trial{1,p}{i,:}(2,:).^2+virtsens_roi_dontcare.trial{1,p}{i,:}(3,:).^2);
                euclid_roi(i,:) = Euklidsche_dist;
                clear Euklidsche_dist
            end
            vs_euclid_roi_dontcare.trial{1,p} = euclid_roi;
            clear euclid_roi
        end

        vs_euclid_roi_dontcare_all.(atlas_downsampled2source_avg.tissuelabel{1,k}) = vs_euclid_roi_dontcare;
        vs_euclid_roi_dontcare_all.(atlas_downsampled2source_avg.tissuelabel{1,k}).time = data_appended.time{1,1};
        clear vs_euclid_roi_dontcare
    end
end

if exist ('vs_euclid_roi_like_all', 'var')
    [max_all_rois_like] =  extract_mean_max(vs_euclid_roi_like_all, 'max', atlas_downsampled2source_avg.tissuelabel); 
    save ([outPath_extdisc 'MEG\sourcespace\run' run '\vs_max_all_rois_like_' freq], 'max_all_rois_like')
    [mean_all_rois_like] =  extract_mean_max(vs_euclid_roi_like_all, 'mean', atlas_downsampled2source_avg.tissuelabel); 
    save ([outPath_extdisc 'MEG\sourcespace\run' run '\vs_mean_all_rois_like_' freq], 'mean_all_rois_like')
end
if exist ('vs_euclid_roi_dislike_all', 'var')
    [max_all_rois_dislike] =  extract_mean_max(vs_euclid_roi_dislike_all, 'max', atlas_downsampled2source_avg.tissuelabel); 
    save ([outPath_extdisc 'MEG\sourcespace\run' run '\vs_max_all_rois_dislike_' freq], 'max_all_rois_dislike')
    [mean_all_rois_dislike] =  extract_mean_max(vs_euclid_roi_dislike_all, 'mean', atlas_downsampled2source_avg.tissuelabel); 
    save ([outPath_extdisc 'MEG\sourcespace\run' run '\vs_mean_all_rois_dislike_' freq], 'mean_all_rois_dislike')
end
if exist ('vs_euclid_roi_dontcare_all', 'var')
    [mean_all_rois_dontcare] =  extract_mean_max(vs_euclid_roi_dontcare_all, 'mean', atlas_downsampled2source_avg.tissuelabel); 
    save ([outPath_extdisc 'MEG\sourcespace\run' run '\vs_mean_all_rois_dontcare_' freq], 'mean_all_rois_dontcare')
    [max_all_rois_dontcare] =  extract_mean_max(vs_euclid_roi_dontcare_all, 'max', atlas_downsampled2source_avg.tissuelabel); 
    save ([outPath_extdisc 'MEG\sourcespace\run' run '\vs_max_all_rois_dontcare_' freq], 'max_all_rois_dontcare')
end


   %%
    
%     adi_sourcestatistics (source_avg, avg_data_like, avg_data_dislike, sourcemodel_lf, hdm_ind_cm, data_bpfreq, template_grid)
    

    
end


function [data_bpfreq_res_sel] = adi_bpfilter(filename, bpname)


switch bpname
    case 'bp1-45Hz'
        bpfreq = 45;
    case 'delta'
        bpfreq = 4;
    case 'theta'
        bpfreq = [4 8];
    case 'alpha'
        bpfreq = [8 13];
    case 'beta'
        bpfreq = [13 25];
    case 'low_gamma'
        bpfreq = [25 45];
    case 'high_gamma'
        bpfreq = [55 90];
end

cfg=[];
cfg.trials        = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') || 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

[data_bpfreq] = ft_preprocessing(cfg, filename);
data_bpfreq.ChannelFlag_Bst = filename.ChannelFlag_Bst;
if isfield(filename, 'trialinfo')
    data_bpfreq.trialinfo = filename.trialinfo;
end

for k = 1:length(data_bpfreq.trial)
    data_bpfreq.trial{1,k} = data_bpfreq.trial{1,k}(1:248,:);
end
    
for k = 1:length(data_bpfreq.trial)
    data_bpfreq.grad.label(249:end) = [];
    data_bpfreq.grad.chanori(249:end, :) = [];
    data_bpfreq.grad.chanpos(249:end, :) = [];
    data_bpfreq.grad.tra(249:end, :) = [];
    data_bpfreq.label(249:end) = [];
end
if 3053 ==length(data_bpfreq.time{1,1}) 
    for k = 1:length(data_bpfreq.trial)
        data_bpfreq.trial{1,k}(:,1) = [];
        data_bpfreq.time{1,k}(:,1) = [];
    end
end

cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);

end
     
function [vs_rois_condition_red] = extract_mean_max(vs_rois_condition, measure, tissuelabel) 
 
for i = 1:length(tissuelabel)
    for p = 1:length(vs_rois_condition.(tissuelabel{i}).trial)
        retval = vs_rois_condition.(tissuelabel{i}).trial{1,p};
        switch measure
            case 'max'
                vs_rois_condition_red.trial{i,p} =  max(retval);        
            case 'mean'
                 vs_rois_condition_red.trial{i,p} =  mean(retval);
        end
    end
 end
vs_rois_condition_red.tissuelabel = tissuelabel;
vs_rois_condition_red.time = vs_rois_condition.(tissuelabel{i}).time;
end

function [] = adi_sourcestatistics(source_avg, avg_data_like, avg_data_dislike, sourcemodel_lf, hdm_ind_cm, data_bpfreq, template_grid)

%%
    
    % Now we make a first call to ft_sourceanalysis in order to compute the spatial filters on the basis of the entire data and keep them in the output for a later use.
    cfg = [];
    cfg.method = 'lcmv';
    cfg.grid = sourcemodel_lf;  % Stefan: cfg.grid = sourcemodel
    cfg.vol = hdm_ind_cm; % Stefan: cfg.headmodel = vol;
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectnoise = 'yes';
    cfg.grid.filter = source_avg.avg.filter;
    cfg.channel = data_bpfreq.grad.label;
    cfg.lcmv.fixedori='no';
    cfg.lcmv.lamda='5%';
%     cfg.lcmv.weightnorm = 'nai'; % nicht klar ob redundant, da cfg.normalize = 'yes' bei ft_prepare_leadfield
%     cfg.lcmv.realfilter   = 'yes'; % keep imaginary part, macht kein Unterschied 
        source_avg_like = ft_sourceanalysis(cfg, avg_data_like);
        source_avg_dislike = ft_sourceanalysis(cfg, avg_data_dislike);

    
    % comput percent change:
    cfg = [];
    cfg.parameter = 'avg.pow';
    cfg.operation = '((x1-x2)./x2)*100';
    percent_change=ft_math(cfg,source_avg_like,source_avg_dislike);
    template_mri = ft_read_mri('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\anatomy\single_subj_T1.nii');
    template_mri_cm = ft_convert_units(template_mri,'cm');

    percent_change.pos=template_grid.pos;
    cfg              = [];
    cfg.voxelcoord   = 'no';
    cfg.parameter    = 'pow';
    cfg.interpmethod = 'nearest';
    percent_change_int  = ft_sourceinterpolate(cfg, percent_change, template_mri_cm);
    
    cfg               = [];
    cfg.method        = 'ortho';
    cfg.funparameter  = 'pow';
    cfg.location = [64 -32 8];
    cfg.funcolormap = 'jet';
    ft_sourceplot(cfg,percent_change_int); % grid überprüfen

% 
%     % statistical analysis: 
%     cfg = [];
%     cfg.parameter    = 'pow';
%     cfg.dim          = sourcemodel_lf.dim;
%     cfg.method           = 'montecarlo';
%     cfg.statistic        = 'ft_statfun_depsamplesT';
%     cfg.correctm         = 'cluster';
% %     cfg.clusterthreshold = 'nonparametric_individual';
%     cfg.clusteralpha     = 0.05;
%     cfg.clusterstatistic = 'maxsum';
%     cfg.tail             = 0;
%     cfg.clustertail      = 0;
%     cfg.alpha            = 0.025;
%     cfg.numrandomization = 1000;
%  
%     trials_like = numel(source_data_like.trial);
%     trials_dislike = numel(source_data_dislike.trial);
%     design  = zeros(2,trials_like+trials_dislike);
%     design(1,1:trials_like) = 1;
%     design(1,trials_like+1:trials_like+trials_dislike) = 2;
%     design(2,1:trials_like) = [1:trials_like];
%     design(2,trials_like+1:trials_like+trials_dislike) = [1:trials_dislike];
%  
%     cfg.design   = design;
%     cfg.ivar     = 1;
%     cfg.uvar     = 2;
%     stat = ft_sourcestatistics(cfg,source_data_like,source_data_dislike); % Anzahl der Trials in beiden Bedingungen muss gleich sein
%     stat.inside=template_grid.inside;
%  
%     cfg            = [];
%     cfg.voxelcoord   = 'no';
%     cfg.parameter    = 'stat';
%     cfg.interpmethod = 'nearest';
%     cfg.downsample = 5; % für plotten nicht downsamplen, sieht nicht gut aus
%     statint  = ft_sourceinterpolate(cfg, stat, template_mri_cm);
%     cfg.parameter    = 'mask';
%     maskint  = ft_sourceinterpolate(cfg, stat, template_mri_cm);
%     statint.mask = maskint.mask;
%     
%     statint.coordsys = 'mni';
%     cfg               = [];
%     cfg.method        = 'ortho';
%     cfg.funparameter  = 'stat';
%     cfg.maskparameter = 'mask';
%     cfg.atlas         = atlas;
%     cfg.location = 'max';
% %     cfg.funcolorlim   = [-5 5];
%     cfg.funcolormap = 'jet';
%     try ft_sourceplot(cfg,statint);
%         
%     catch
%         warning('sourceplot funktioniert nicht')
%     end
%     


end

