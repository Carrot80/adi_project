function   [virtsens_like, virtsens_dislike] = adi_source_reconstruction_allRuns(path2vol, path2data, mriPath, SubjOutPath, outPath_extdisc, freqbandname, name_condition1, name_condition2, latency, cfg_virtsens)

% calculations spatial filter, multiplies it with avgdata and runs Support
% Vector machine

time =  latency(1,:);

% if ~exist([SubjOutPath '\run1\SVM_result\virtSens\virtsens_like_vs_dislike_1_' freqbandname '.mat'])
    
    [virtsens_like] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, 'allRuns', freqbandname, name_condition1, outPath_extdisc);                                                                   
    [virtsens_dislike] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, 'allRuns', freqbandname, name_condition2, outPath_extdisc);

% end

end

function [virtsens] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, run, freq, condition, outPath_extdisc)

cond = [condition '_' run];   
% load data:
if exist([path2data, '03_appended_data\' cond '_' freq '.mat'], 'file')
   data = load ([path2data, '03_appended_data\' condition '_' run '_' freq '.mat']);
else
    return
end
if exist([path2data '02b_bpfreq\' condition '500_2_' freq '.mat'], 'file')
    load([path2data '02b_bpfreq\' condition '500_2_' freq '.mat'])
else
    load([path2data '02b_bpfreq\' condition '500_1_' freq '.mat'])
end
data.(cond).grad = data_bpfreq.grad;
clear data_bpfreq

for k = 1:length(data.(cond).trial)
    data.(cond).grad.label(249:end) = [];
    data.(cond).grad.chanori(249:end, :) = [];
    data.(cond).grad.chanpos(249:end, :) = [];
    data.(cond).grad.tra(249:end, :) = [];
    data.(cond).label(249:end) = [];
end
    
for k = 1:length(data.(cond).trial)
    data.(cond).trial{1,k} = data.(cond).trial{1,k}(1:248,:);
end
   
% load template vol:
if ~exist('template_grid', 'var')
    load('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\headmodel\standard_singleshell');
    cfg = [];
    cfg.grid.xgrid  = -20:1:20;
    cfg.grid.ygrid  = -20:1:20;
    cfg.grid.zgrid  = -20:1:20;
    cfg.grid.unit   = 'cm';
    cfg.grid.tight  = 'yes';
    cfg.inwardshift = -1.5;
    cfg.headmodel  = vol;
    template_grid  = ft_prepare_sourcemodel(cfg);
end
    
figure;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));
hold on
ft_plot_vol(vol,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; 
camlight;

%% Load atlas and create a binary mask
if ~exist('mask', 'var')
    atlas = ft_read_atlas('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\atlas\aal\ROI_MNI_V4.nii');
    atlas = ft_convert_units(atlas,'cm'); % assure that atlas and template_grid are expressed in the %same units
    cfg = [];
    cfg.atlas = atlas;
    cfg.roi = atlas.tissuelabel;
    cfg.inputcoord = 'mni';
    mask = ft_volumelookup(cfg, template_grid);
end
    
    % create temporary mask according to the atlas entries
tmp                  = repmat(template_grid.inside,1,1);
tmp(tmp==1)          = 0;
tmp(mask)            = 1;

% define inside locations according to the atlas based mask
template_grid.inside = tmp;

% plot the atlas based grid
figure;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));
view ([-90 -90 -90])

%% Inverse-warp the subject specific grid to the atlas based template grid

% For this step the individual volume is required:
if ~exist('sourcmodel', 'var')
    load ([mriPath, 'mri_realigned.mat'])
    cfg            = [];
    cfg.resolution = 1;
    cfg.dim        = [256 256 256];
    mri_resliced = ft_volumereslice(cfg, mri_realigned)
    mri_resliced = ft_convert_units(mri_resliced, 'cm')
    grad = ft_convert_units(data.(cond).grad, 'cm');
    data.(cond).grad = grad;

    cfg                = [];
    cfg.grid.warpmni   = 'yes';
    cfg.grid.template  = template_grid;
    cfg.grid.nonlinear = 'yes';
    cfg.mri            = mri_resliced; % nicht sicher
    sourcemodel        = ft_prepare_sourcemodel(cfg);
end
if ~exist('hdm_ind', 'var')
    hdm_ind = load ([path2vol, 'vol.mat'])
    hdm_ind = ft_convert_units(hdm_ind.vol, 'cm'); 
end
% Plot the final source model together with the individual head model and
% the sensor array:
close all
% the 4D/bti system is expressed in units of 'm', therefore we force all geometrical objects to have the same unit
% vol_indiv = ft_convert_units(hdm, 'm');
% sourcemodel = ft_convert_units(sourcemodel, 'm');

figure; hold on     % plot all objects in one figure

ft_plot_vol(hdm_ind,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; %camlight;
alpha 0.4           % make the surface transparent

ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:),'vertexcolor','b');% plot only locations inside the volume

ft_plot_sens(data.(cond).grad,'style','g*', 'coil', 1);% plot the sensor array

view ([0 -90 0])


savefig([SubjOutPath 'check_vol_sourcemodel_grad.fig']);
fig = ([SubjOutPath 'check_vol_sourcemodel_grad']);
print('-dpng', fig); 
close all


%We first create the leadfield using ft_prepare_leadfield using the individual head model from the previous step, the sensor array and the sourcemodel.
cfg                 = [];
cfg.channel         = data.(cond).label;% ensure that rejected sensors are not present
cfg.grad            = data.(cond).grad;
cfg.vol             = hdm_ind;
cfg.lcmv.reducerank = 2; % default for MEG is 2, for EEG is 3
cfg.grid = sourcemodel;
% cfg.grid.pos=rois; 
[grid] = ft_prepare_leadfield(cfg);

    %%
% 
%     cfg = [];
%     cfg.toilim = prestim;
%     datapre = ft_redefinetrial(cfg, data.(cond));
%     cfg.toilim = poststim;
%     datapost = ft_redefinetrial(cfg, data.(cond));

cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-.5 1]; % oder von 0 bis 1?
cfg.keeptrials = 'yes';
avg_data.(cond) = ft_timelockanalysis(cfg, data.(cond));
%     
%     cfg = [];
%     cfg.covariance='yes';
%     cfg.keeptrials='yes';
%     avgpre = ft_timelockanalysis(cfg, datapre);
%     avgpst = ft_timelockanalysis(cfg, datapost);

    % Now we make a first call to ft_sourceanalysis in order to compute the spatial filters on the basis of the entire data and keep them in the output for a later use.
cfg = [];
cfg.method = 'lcmv';
cfg.grid = grid;  % Stefan: cfg.grid = sourcemodel
cfg.vol = hdm_ind; % Stefan: cfg.headmodel = vol;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.channel = data.(cond).label;
cfg.lcmv.fixedori='yes';
cfg.lcmv.lamda='5%';
source_avg = ft_sourceanalysis(cfg, avg_data.(cond));
    
%     cfg = [];
%     cfg.method = 'lcmv';
%     cfg.grid = grid;
%     cfg.grid.filter = source_avg.avg.filter;
%     cfg.headmodel = hdm_ind;
%     cfg.lcmv.projectnoise = 'yes';
%     source_pre = ft_sourceanalysis(cfg, avgpre);
%     source_post = ft_sourceanalysis(cfg, avgpst);
  
    
%     outPath_extdisc_subj_time = [ outPath_extdisc 'sourcespace\_run' run filesep prestim '_' poststim ];
%     if ~exist('outPath_extdisc_subj_time', 'dir')
%         mkdir (outPath_extdisc_subj_time)
%     end
%     save([outPath_extdisc_subj_time 'source_pre' ], 'source_pre')
%     save([outPath_extdisc_subj_time 'source_post' ], 'source_post')
    
    spatialfilter_orig = source_avg.avg.filter;
    save([outPath_extdisc 'spatialfilter_fixed_orientation_singletrials_' condition '_' freq], 'spatialfilter_orig')
    
%     cfg = [];
%     cfg.parameter = 'avg.pow';
%     cfg.operation = '((x2-x1)./x1)*100';
%     source_ratio = ft_math(cfg, source_pre, source_post);

    
    % On the basis of the computed filters, kept in the output, it is now possible to multiply them with the data. This operation will yield time series commonly known as virtual sensors.
    
    spatialfilter = cat(1,source_avg.avg.filter{:});
    virtsens = [];
    for i = 1:length(data.(cond).trial)
        virtsens.trial{i} = spatialfilter*data.(cond).trial{i};
    end;
    virtsens.time=data.(cond).time;

    for i=1:length(virtsens.trial{1}(:,1))
        virtsens.label{i}=num2str(i);
    end;

    % siehe Skript von Yuval: => all will have a similar noise level
%     ns = mean(abs(spatialfilter),2);
% 
%     for k = 1:length(virtsens.trial)
%         virtsens_ns.trial{k} = virtsens.trial{1,k}./repmat(ns,1,size(virtsens.trial{1,k},2)); % => all will have a similar noise level
%     end
% 
%     virtsens_ns.time = virtsens.time;
%     virtsens_ns.fsample = virtsens.fsample;
%     
%     for i = 1:length(virtsens_ns.trial{1}(:,1))
%         virtsens_ns.label{i} = num2str(i);
%     end

    
    
%     cfg=[];
%     avg = ft_timelockanalysis(cfg, virtsens)
%     figure
%     plot(avg_data.(cond).time, avg_data.(cond).avg(1:248, :))
%     plot(virtsens_ns.time{1,1}, avg.avg(3,:))
%     save ([outPath 'virtsens_ns_' condition '_' freq '.mat'], 'virtsens_ns')
%     save ([outPath 'virtsens_' condition '_' freq '.mat'], 'virtsens')

%     %% reduce the source reconstructed data to the dominant orientation
%     cfg = [];
%     cfg.projectmom = 'yes';
%     cfg.keepnoisemom = 'yes';
%     source_proj = ft_sourcedescriptives(cfg, source_post); % mir ist unklar, worauf dies hinausläuft
% 
%     load('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180305\template\sourcemodel/standard_sourcemodel3d10mm');
%     template_grid = sourcemodel;
%     clear sourcemodel;
%     template_grid = ft_convert_units(template_grid, 'cm');
%     
%     templatefile = fullfile('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180305\template\anatomy\single_subj_T1.nii');
%     template_mri = ft_read_mri(templatefile);
%     source_proj.pos = template_grid.pos;
%     cfg              = [];
%     cfg.voxelcoord   = 'no';
%     cfg.parameter    = parameter; %?
%     cfg.interpmethod = 'linear';
%     source_int  = ft_sourceinterpolate(cfg, source_proj, template_mri);
%     source_int.coordsys = 'mni';
%     
%     cfg               = [];
%     cfg.method        = 'ortho';
%     cfg.funparameter  = parameter;
%     mask = source_int.(cfg.funparameter);
%     mask(isnan(mask)) = 0;
%     mask(abs(mask)>0)      = 1;
%         
%     if ~isnan(thresholdMin) 
%         fun = source_int.(cfg.funparameter);
%         mask(fun > thresholdMin & fun <= 0) = 0;
%     end
%     
%     if ~isnan(thresholdMax) 
%         fun = source_int.(cfg.funparameter);
%         mask(fun < thresholdMax & fun >= 0) = 0;
%     end
%     
%     source_int.mask = mask;
%     cfg.maskparameter = 'mask';
%     cfg.opacitylim = [0.9 1];
%     cfg.funcolorlim   = op;
%     cfg.funcolormap = 'jet';
%     cfg.atlas = fullfile('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180305\template\atlas\aal\ROI_MNI_V4.nii');
%     
%     if ~isempty(location)
%         cfg.location = location;
%     end
%     
%     ft_sourceplot(cfg, source_int);
% 
%     cfg              = [];
%     cfg.parameter    = parameter;
%     cfg.interpmethod = 'linear';
%     source_int  = ft_sourceinterpolate(cfg, source_proj, mri);
%     
%     %% plot the neural activity index (power/noise)
%     cfg               = [];
%     cfg.method        = 'ortho';
%     cfg.funparameter  = parameter;
%     mask = source_int.(cfg.funparameter);
%     mask(isnan(mask)) = 0;
%     mask(mask>0)      = 1;
%     
%     if ~isnan(thresholdMin) 
%         fun = source_int.(cfg.funparameter);
%         mask(fun > thresholdMin & fun <= 0) = 0;
%     end
%     
%     if ~isnan(thresholdMax) 
%         fun = source_int.(cfg.funparameter);
%         mask(fun < thresholdMax & fun >= 0) = 0;
%     end
%     
%     source_int.mask = mask;
%     cfg.maskparameter = 'mask';
%     cfg.opacitylim = [0.9 1];
%     cfg.funcolorlim   = op;
%     cfg.funcolormap   = 'jet';
%     
%     ft_sourceplot(cfg, source_int);
    
end

