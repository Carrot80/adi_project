
function adi_source_SVM (path2vol, path2data, mriPath, outPath_extdisc, freqbandname, like, dislike, latency)

% calculations spatial filter per Run and Subject, 
time =  latency(1,:);

%% run 1:
if ~exist([outPath_extdisc 'MEG\sourcespace\ROIs\run1\spatialfilter_loose_orientation_singletrials_like_' freqbandname '.mat'], 'file')
  
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath,  num2str(1), freqbandname, 'like', outPath_extdisc)                                                                   
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(1), freqbandname, 'dislike', outPath_extdisc)

    clear virtsens_like virtsens_ns_like virtsens_dislike virtsens_ns_dislike

end

%% run2:

if ~exist([outPath_extdisc 'MEG\sourcespace\run2\spatialfilter_loose_orientation_singletrials_like_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(2), freqbandname, 'like', outPath_extdisc)
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(2), freqbandname, 'dislike', outPath_extdisc)

    clear virtsens_like virtsens_ns_like virtsens_dislike virtsens_ns_dislike
end

%% run3:

if ~exist([outPath_extdisc 'MEG\sourcespace\run3\spatialfilter_loose_orientation_singletrials_like_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(3), freqbandname, 'like', outPath_extdisc)
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(3), freqbandname, 'dislike', outPath_extdisc)

    clear virtsens_like virtsens_ns_like virtsens_dislike virtsens_ns_dislike

end

end

function [virtsens_ns] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, run, freq, condition, outPath_extdisc)
    
% load data:
    if exist([path2data, condition '500_' run '.mat'], 'file')
        load ([path2data, condition '500_' run '.mat'], 'cleanMEG_interp')
    else
       return
    end
    
    [data_bpfreq] = adi_bpfilter(cleanMEG_interp, freq);
    
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
   
  %% Reconstruct the time course of activity at a particular brain location
   template_mri = ft_read_mri('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\anatomy\single_subj_T1.nii');
   template_mri_cm = ft_convert_units(template_mri,'cm')
   
load('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template/headmodel/standard_singleshell');
 
cfg = [];
cfg.grid.xgrid  = -20:1:20;
cfg.grid.ygrid  = -20:1:20;
cfg.grid.zgrid  = -20:1:20;
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.headmodel        = vol;
template_grid  = ft_prepare_sourcemodel(cfg);
 
figure;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));
hold on
ft_plot_vol(vol,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;

atlas = ft_read_atlas('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\atlas\aal\ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas,'cm'); % assure that atlas and template_grid are expressed in the %same unitscfg = [];

cfg = []
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel;
cfg.inputcoord = 'mni';
mask = ft_volumelookup(cfg, template_grid);
 
% create temporary mask according to the atlas entries
tmp                  = repmat(template_grid.inside,1,1);
tmp(tmp==1)          = 0;
tmp(mask)            = 1;
 
% define inside locations according to the atlas based mask
template_grid.inside = tmp;
 
% plot the atlas based grid
figure;ft_plot_mesh(template_grid.pos(template_grid.inside,:));
% 
% inverse-warp the subject specific grid to the atlas based template grid
% For this step the individual volume is required.
load ([mriPath, 'mri_realigned.mat']);
mri_realigned_cm = ft_convert_units(mri_realigned, 'cm');
cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes';
cfg.mri            = mri_realigned_cm;
sourcemodel        = ft_prepare_sourcemodel(cfg);

% Plot the final source model together with the individual head model and the sensor array
close all

% the 4D/bti system is expressed in units of 'm', therefore we force all geometrical objects to have the same unit


% Finally, it is wise to check whether all computed objects align well with one another, i.e. whether the grid is correctly 
% placed within the volume conduction model, which both have to be aligned with the MEG sensors. 
% Note that all objects that we plot need to be expressed in the same units and the same coordinate space. 
% Here, we need to transform the head model from 'mm' into 'cm'.

    hdm = load([path2vol 'vol.mat']);
    hdm_cm = ft_convert_units(hdm.vol, 'cm');
    sourcemodel = ft_convert_units(sourcemodel, 'cm');
    
    figure; hold on     % plot all objects in one figure
    ft_plot_vol(hdm_cm, 'facecolor', 'cortex', 'edgecolor', 'none')
    alpha 0.4           % make the surface transparent
    ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); alpha 0.5; %camlight;% plot only locations inside the volume
    hold on
    data_bpfreq.grad = ft_convert_units(data_bpfreq.grad, 'cm');
    ft_plot_sens(data_bpfreq.grad, 'style','*r');
    view ([0 -90 0])
    savefig([outPath_extdisc 'MEG\sourcespace\run' run '\sensor_position_check.fig'])

% %     Compute the leadfield
% We first create the leadfield using ft_prepare_leadfield using the individual head model from the previous step, the sensor array and the sourcemodel.

    cfg                 = [];
    cfg.channel         = data_bpfreq.label;% ensure that rejected sensors are not present
    cfg.grad            = data_bpfreq.grad;
    cfg.vol             = hdm.vol;
    cfg.lcmv.reducerank = 2; % default for MEG is 2, for EEG is 3
    cfg.grid = sourcemodel;
    [grid] = ft_prepare_leadfield(cfg);

    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = [-.5 1]; % oder von 0 bis 1?
    cfg.keeptrials = 'yes';
    avg_data_bpfreq = ft_timelockanalysis(cfg, data_bpfreq); 
    
    cfg=[];
    cfg.method='lcmv';
    cfg.grid=grid;
    cfg.vol=hdm_cm;
    cfg.lcmv.keepfilter='yes';
    cfg.lcmv.projectnoise = 'yes';
    cfg.channel = data_bpfreq.label;
    cfg.lcmv.lamda='5%';
    cfg.rawtrial = 'yes';
    cfg.lcmv.fixedori='no';
    cfg.lcmv.weightnorm = 'nai';
    cfg.lcmv.realfilter   = 'yes';
    source_avg=ft_sourceanalysis(cfg, avg_data_bpfreq);

source_avg.pos=template_grid.pos;
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, source_avg, template_mri);
cfg=[];
source_int_parcelled = ft_sourceparcellate(cfg, source_int, atlas);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel{1,roi_num};
mask_roi_num = ft_volumelookup(cfg, template_mri_cm);
x_roi = find(ismember(atlas.tissuelabel,atlas.tissuelabel{roi_num}));
indx_roi = find(atlas.tissue==x_roi);

%  Finally, we can plot the result using ft_sourceplot.

% cfg               = [];
% cfg.method        = 'ortho';
% cfg.funparameter  = 'pow';
% % cfg.location = [64 -32 8];
% cfg.funcolormap = 'jet';
% ft_sourceplot(cfg,source_int);

cfg=[];



cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [-30 30];
cfg.funcolormap    = 'jet'; 
cfg.projmethod     = 'nearest'; 
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projthresh     = 0.8;
cfg.camlight       = 'no';
ft_sourceplot(cfg, source_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull

   for roi_num = 1:90
       cfg = [];
       cfg.inputcoord = 'mni';
       cfg.atlas = atlas;
       cfg.roi = atlas.tissuelabel{1,roi_num};
       mask_roi_num = ft_volumelookup(cfg, template_mri_cm);
       x_roi = find(ismember(atlas.tissuelabel,atlas.tissuelabel{roi_num}));
       indx_roi = find(atlas.tissue==x_roi);
       
   end    
     %% Computing source modell: 
     % Following the construction of the volume conduction model (singleshell), we need to discretize the brain into a source model or grid.
     % For each grid point in the brain, the lead field matrix is calculated later When constructing the source model, you might want to keep in mind 
     % that averaging and statistics over subjects can only be done if the individual subjects source reconstructed results are mapped onto a common space. 
     % Now, we will construct a regular grid in MNI template space and spatially deform this grid to the individual subjects brain. The following code loads the template grid that is included in the FieldTrip release:

     % load template grid:

     template_sourcemodel = load('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\sourcemodel\standard_sourcemodel3d8mm');
     plot3(template_sourcemodel.sourcemodel.pos((template_sourcemodel.sourcemodel.inside==1),1), template_sourcemodel.sourcemodel.pos((template_sourcemodel.sourcemodel.inside==1),2), template_sourcemodel.sourcemodel.pos((template_sourcemodel.sourcemodel.inside==1),3))
     ft_plot_mesh(template_sourcemodel.sourcemodel, 'edgecolor', 'none')
     template_grid=ft_convert_units(template_sourcemodel.sourcemodel,'cm');% ensure no unit mismatch
     % Next, we normalise the individual MRI to derive parameters allowing to convert the mni- coordinates of the desired parcels into individual coordinates. 
     % For this we use ft_warp_apply:
     load ([mriPath, 'mri_realigned.mat']);
     mri_realigned_cm = ft_convert_units(mri_realigned, 'cm');
     norm_mri_realigned=ft_volumenormalise([],mri_realigned_cm);
%      cfg = []; 
%     cfg.interpmethod = 'nearest'; 
%     cfg.parameter = 'all'; 
%     stat_atlas = ft_sourceinterpolate(cfg, template_grid, atlas);
%      
    pos_roi=template_grid.pos(indx_roi,:);% xyz positions in mni coordinates
    posHGL=template_grid.pos(indxHGL,:);% xyz positions in mni coordinates
    posHGR=template_grid.pos(indxHGR,:);% xyz positions in mni coordinates

    posback=ft_warp_apply(norm.params,posCML,'sn2individual');
    btiposCML= ft_warp_apply(pinv(norm.initial),posback);% xyz positions in individual coordinates

    posback=ft_warp_apply(norm.params,posHGL,'sn2individual');
    btiposHGL= ft_warp_apply(pinv(norm.initial),posback);% xyz positions in individual coordinates

    posback=ft_warp_apply(norm.params,posHGR,'sn2individual');
    btiposHGR= ft_warp_apply(pinv(norm.initial),posback);% xyz positions in individual coordinates
    Now we create a source model for these particular locations only.

    cfg=[];
    cfg.vol=hdm;
    cfg.channel=dataica.label;  
    cfg.grid.pos=[btiposCML;btiposHGL;btiposHGR]./1000;% units of m
    cfg.grad=dataica.grad;
    sourcemodel_virt=ft_prepare_leadfield(cfg);

     % inverse-warp the template grid to subject specific coordinates 
    cfg                = [];
    cfg.grid.warpmni   = 'yes';
    cfg.grid.template  = template_sourcemodel.sourcemodel;
    cfg.grid.nonlinear = 'yes'; % use non-linear normalization
    cfg.mri            = mri_realigned;
    sourcemodel_inversewarped = ft_prepare_sourcemodel(cfg); % sourcemodel = grid
    plot3(sourcemodel_inversewarped.pos((sourcemodel_inversewarped.inside==1),1), sourcemodel_inversewarped.pos((sourcemodel_inversewarped.inside==1),2), sourcemodel_inversewarped.pos((sourcemodel_inversewarped.inside==1),3))
    


    cfg             = [];
    cfg.grid        = sourcemodel_inversewarped;
    cfg.headmodel   = hdm_cm;
    cfg.channel     = data_bpfreq.label;
    cfg.grad        = data_bpfreq.grad;
    cfg.normalize = 'yes'; % Tiefennormalisierung: If you are not contrasting the activity of interest against another condition or baseline time-window, then you may choose to normalize the lead field in this step (cfg.normalize='yes'), which will help control against the power bias towards the center of the head.
    cfg.reducerank  = '2'; % or number (default = 3 for EEG, 2 for MEG
    sourcemodel_lf  = ft_prepare_leadfield(cfg, data_bpfreq); % sourcemodel = grid; entpricht leadfield der Leitfähigkeit?
    plot3(sourcemodel_lf.pos((sourcemodel_lf.inside==1),1), sourcemodel_lf.pos((sourcemodel_lf.inside==1),2), sourcemodel_lf.pos((sourcemodel_lf.inside==1),3))

    %%
    
    
   
    
    % Determine the index of the label of interest:
    x = find(ismember(atlas.tissuelabel,'Heschl_L'));
    % and determine the index points of locations within the desired
    % parcel:
    indxHGL = find(atlas.tissue==x);
    % These steps can be repeated for all desired parcels. In the present case the ramining two:

    x=find(ismember(atlas.tissuelabel,'Heschl_R'));
    indxHGR = find(stat_atlas.tissue==x); 

    x=find(ismember(atlas.tissuelabel,'Cingulum_Mid_L'));
    indxCML = find(stat_atlas.tissue==x); 


%%
    
    
%     The grid of estimated power values can be plotted superimposed on the anatomical MRI.
%     This requires the output of ft_sourceanalysis to match position of the template MRI to
%     which we warped the sourcemodel. Because of this warping that we already did, we can simply overwrite the grid position information without any further mathematical operation:
%  
    outPath_extdisc_subj = [ outPath_extdisc 'MEG\sourcespace\ROI\run' run filesep];
    if ~exist(outPath_extdisc_subj, 'dir')
        mkdir (outPath_extdisc_subj)
    end
    
%     spatialfilter_orig = source_avg.avg.filter; 
%     save([outPath_extdisc_subj 'spatialfilter_loose_orientation_singletrials_' condition '_' freq], 'spatialfilter_orig')
   
    %%
    spatialfilter = cat(1,source_avg.avg.filter{:});
    
    virtsens = [];
    for i = 1:length(data_bpfreq.trial)
        virtsens.trial{i} = spatialfilter*data_bpfreq.trial{i};
    end
    virtsens.time=data_bpfreq.time;
    virtsens.fsample=data_bpfreq.fsample;

    for i=1:length(virtsens.trial{1}(:,1))
        virtsens.label{i}=num2str(i);
    end
    
    ns = mean(abs(spatialfilter),2);

    for k = 1:length(virtsens.trial)
        virtsens_ns.trial{k} = virtsens.trial{1,k}./repmat(ns,1,size(virtsens.trial{1,k},2)); % => all will have a similar noise level
    end

    virtsens_ns.time = virtsens.time;
    virtsens_ns.fsample = virtsens.fsample;

    for k = 1:length(virtsens_ns.trial{1}(:,1))
        virtsens_ns.label{k} = num2str(k);
    end
    virtsens = virtsens_ns;
    clear virtsens_ns
    
    
    %% 
%     https://mailman.science.ru.nl/pipermail/fieldtrip/2015-December/009880.html
    % you can simply interpolate your virtual channels to the standard MRI and use source plot with an atlas to find the corresponding label.
%     https://mailman.science.ru.nl/pipermail/fieldtrip/2014-July/008211.html
    
    
    %1. Compute the virtual channels for all grid points in the cortex:
    % done
%     2. If you have computed the virtual channels on an individual MRI, make sure that the positions correspond to the standard MRI (I think that's described in the source analysis-tutorial)
    % sollte ich getan haben, evtl. nochmal plotten
    
   % To find your ROI, you use ft_volumelookup to get a binary matrix which tells you, which voxels of the standard MRI are in your ROI.
   
     
   % 3. Build a high-resolution 3D-Grid on the standard MRI (It's important to use one grid point per voxel):
   
  vox1 = template_mri_cm.transform*[template_mri_cm.dim,1]'; % start
%   vox1 = vox1(1:3)*10; % convert to mm
    vox2 = template_mri_cm.transform*[1,1,1,1]';   %end
%     vox2 = vox2(1:3)*10; % convert to mm
% 
% % Choose the Minima
vox_min = min([vox1,vox2]')';
vox_max = max([vox1,vox2]')';
% 
% % Set the spacing
vox_delta = abs(template_mri_cm.transform*[1,1,1,1]'-template_mri_cm.transform*[2,2,2,1]');
% vox_delta = vox_delta(1:3)*10; % convert to mm
RES = vox_delta(1);

cfg = [];
%cfg.inwardshift = -0.5;
cfg.grid.xgrid = vox_min(1):RES:vox_max(1)+1;
cfg.grid.ygrid = vox_min(2):RES:vox_max(2)+1;
cfg.grid.zgrid = vox_min(3):RES:vox_max(3)+1;
cfg.grid.tight = 'no';
template_grid_cm = ft_prepare_sourcemodel(cfg, hdm_cm);
  
   
   % 4. Select all grid points from the grid of step 3 which belong to an atlas (see ft_volumelookup)
  
   
   % 5. Use pythagoras to find the virtual channels from step 2 closest to
   % the grid points from step 4. Using pythagoras, you can now find the closest high-resolution grid points for each of your virtual electrodes. As soon as you have these, you can see which of them have a value of 1 in your binary matrix of ft_volumelookup. Now you know which virtual electrodes are inside your ROI.
    
%    So for example, you have 2000 virtual electrodes on a -20:1:20 regular grid. You then compute a binary matrix and a high-resolution grid based on the standard MRI 
%    resulting in ~900000 grid points on a 91x109x91 grid. Using pythagoras, you find that grid point 6000
%    is closest to virtual electrode #1. Then you can see if element 6000 from your binary matrix is a 1. 
%    If so, your virtual electrode #1 is in your ROI, if not, your virtual electrode is not in your ROI.
   
cfg=[];
parcel = ft_sourceparcellate(cfg, virtsens.trial{1,1}, atlas);
parcelmask = ft_sourceparcellate(cfg, maskint, atlas);
%% create dummy struct
dummy=atlas;
dummymask = atlas;
for i=1:length(parcel.stat)
      dummy.tissue(find(dummy.tissue==i))=parcel.stat(i);
      dummymask.tissue(find(dummymask.tissue==i))=parcelmask.mask(i);
end;
%% plot the result
statint.parcel=dummy.tissue;
statint.coordsys = 'mni';
statint.mask  = dummymask.tissue;
cfg=[];
cfg.method = 'slice';
cfg.funparameter = 'parcel';
cfg.funcolormap    = 'jet';
cfg.maskparameter = 'mask';
cfg.renderer = 'zbuffer';
cfg.funcolorlim   = [-5 5];
cfg.atlas = atlas;
ft_sourceplot(cfg,statint);


   
     
%% alt:
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
    
    
    
%     cfg = [];
%     cfg.parameter = 'avg.pow';
%     cfg.operation = '((x2-x1)./x1)*100';
%     source_ratio = ft_math(cfg, source_pre, source_post);

    
    % On the basis of the computed filters, kept in the output, it is now possible to multiply them with the data. This operation will yield time series commonly known as virtual sensors.
%     
%     spatialfilter = cat(1,source_avg.avg.filter{:});
% %     spatialfilter = spatialfilter(:, 1:248);
% %     
% %     for k = 1:length(data_bpfreq.trial)
% %         data_bpfreq.trial{1,k} = data_bpfreq.trial{1,k}(1:248,:)
% %     end
%     
%     virtsens = [];
%     for i = 1:length(data_bpfreq.trial)
%         virtsens.trial{i} = spatialfilter*data_bpfreq.trial{i};
%     end
%     virtsens.time=data_bpfreq.time;
%     virtsens.fsample=data_bpfreq.fsample;
% 
%     for i=1:length(virtsens.trial{1}(:,1))
%         virtsens.label{i}=num2str(i);
%     end
% 
%     % siehe Skript von Yuval: => all will have a similar noise level
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
% 
%     for i=1:length(virtsens_ns.trial)
%         mean_VS_denoised.trial{i} = mean(virtsens_ns.trial{i})
%     end
    
    
%     cfg=[];
%     avg = ft_timelockanalysis(cfg, virtsens)
%     figure
%     plot(avg_data_bpfreq.time, avg_data_bpfreq.avg(1:248, :))
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


function [data_bpfreq] = adi_bpfilter(filename, bpname)


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

       
end
     


