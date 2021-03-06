
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
    cfg.headmodel        = template_vol.vol;
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
    data_bpfreq.grad = ft_convert_units(data_bpfreq.grad, 'cm');
    hold on
    ft_plot_sens(data_bpfreq.grad);
    savefig([outPath_extdisc 'MEG\sourcespace\run' run '\sensor_position_check.fig'])

    cfg             = [];
    cfg.grid        = sourcemodel_inversewarped;
    cfg.headmodel   = hdm_ind_cm;
    cfg.vol   = hdm_ind_cm;
    cfg.channel     = data_bpfreq.label;
    cfg.grad        = data_bpfreq.grad;
    cfg.normalize = 'yes'; % If you are not contrasting the activity of interest against another condition or baseline time-window, then you may choose to normalize the lead field in this step (cfg.normalize='yes'), which will help control against the power bias towards the center of the head.
    cfg.reducerank  = '2'; % or number (default = 3 for EEG, 2 for MEG
    sourcemodel_lf  = ft_prepare_leadfield(cfg, data_bpfreq); % sourcemodel = grid

    
      figure;ft_plot_mesh(sourcemodel_inversewarped.pos(sourcemodel_inversewarped.inside,:));
 

    %%
    
%     
%     % Determine the index of the label of interest:
%     x = find(ismember(atlas.tissuelabel,'Heschl_L'));
%     % and determine the index points of locations within the desired
%     % parcel:
%     indxHGL = find(atlas.tissue==x);
%     % These steps can be repeated for all desired parcels. In the present case the ramining two:
% 
%     x=find(ismember(atlas.tissuelabel,'Heschl_R'));
%     indxHGR = find(stat_atlas.tissue==x); 
% 
%     x=find(ismember(atlas.tissuelabel,'Cingulum_Mid_L'));
%     indxCML = find(stat_atlas.tissue==x); 

%%
% cfg = [];            
% cfg.toilim = [-.4 -.05];
% datapre = ft_redefinetrial(cfg, data_bpfreq); 
% cfg.toilim = [.05 .4];
% datapost = ft_redefinetrial(cfg, data_bpfreq); 

cfg = [];
cfg.covariance='yes';
cfg.covariancewindow = [0 1];
cfg.keeptrials = 'yes';
avg = ft_timelockanalysis(cfg,data_bpfreq);
%  
% cfg = [];
% cfg.covariance='yes';
% avgpre = ft_timelockanalysis(cfg,datapre);
% avgpst = ft_timelockanalysis(cfg,datapost);

cfg=[];
cfg.method='lcmv';
cfg.grid=sourcemodel_lf;
cfg.vol=hdm_ind_cm;
cfg.lcmv.keepfilter='yes';
cfg.channel = data_bpfreq.label;
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.fixedori='no';
cfg.lcmv.lamda='5%';
source_avg=ft_sourceanalysis(cfg, avg);
    
% cfg=[];
% cfg.method='lcmv';
% cfg.grid=sourcemodel_lf;
% cfg.grid.filter=sourceavg.avg.filter;
% cfg.vol=hdm_ind_cm;
% sourcepreS1=ft_sourceanalysis(cfg, avgpre);
% sourcepstS1=ft_sourceanalysis(cfg, avgpst);

% cfg = [];
% cfg.parameter = 'avg.pow';
% cfg.operation = '((x1-x2)./x2)*100';
% S1bl=ft_math(cfg,sourcepstS1,sourcepreS1);

% template_mri = ft_read_mri('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\anatomy\single_subj_T1.nii');
% template_mri_cm = ft_convert_units(template_mri,'cm')
% S1bl.pos=template_grid.pos;
% cfg              = [];
% cfg.voxelcoord   = 'no';
% cfg.parameter    = 'pow';
% cfg.interpmethod = 'nearest';
% source_int  = ft_sourceinterpolate(cfg, S1bl, template_mri);

% cfg               = [];
% cfg.method        = 'ortho';
% cfg.funparameter  = 'pow';
% cfg.location = [64 -32 8];
% cfg.funcolormap = 'jet';
% ft_sourceplot(cfg,source_int);

templatefile = '\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607/external/spm8/templates/T1.nii';
template_mri_spm = ft_read_mri(templatefile);
template_mri_spm_cm = ft_convert_units(template_mri_spm, 'cm');
% 
% cfg=[];
% cfg.downsample = 5 %integer number (default = 1, i.e. no downsampling)
% [template_mri_down_cm] = ft_volumedownsample(cfg, template_mri_spm_cm)
% atlas_down = ft_volumedownsample(cfg, atlas)
source_avg.pos = template_grid.pos;
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, source_avg, template_mri_spm_cm);
%%
cfg=[];
parcel = ft_sourceparcellate(cfg, source_int, atlas);

%%
dummy=atlas;
for i=1:length(parcel.pow)
      dummy.tissue(find(dummy.tissue==i))=parcel.pow(i);
end;
%%
source_int.parcel=dummy.tissue;
source_int.coordsys = 'mni';

cfg=[];
cfg.method = 'ortho';
cfg.funparameter = 'parcel';
cfg.funcolormap    = 'jet';
cfg.renderer = 'zbuffer';
% cfg.location = [-42 -20 6];
cfg.atlas = atlas;
% cfg.funcolorlim = [-30 30];
ft_sourceplot(cfg,source_int);

%Alternatively, the maximal activity in the left Heschl gyrus can be plotted on the brain surface as follows.
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
% view ([-70 20 50])
% light ('Position',[-70 20 50])
material dull
% bis hierhin scheint alles zu stimmen

%%


cfg = []; 
cfg.interpmethod = 'nearest'; 
cfg.parameter = 'tissue'; 
source_atlas = ft_sourceinterpolate(cfg, atlas, source_avg);
cfg.parameter = 'all'; 
atlas_source_atlas = ft_sourceinterpolate(cfg, source_avg, atlas);
% Determine the index of the label of interest,
x = find(ismember(atlas.tissuelabel,'Heschl_L'));
%and determine the index points of locations within the desired parcel.
indxHGL = find(source_atlas.tissue==x)


% These steps can be repeated for all desired parcels. In the present case the ramining two:
x = find(ismember(atlas.tissuelabel,'Frontal_Sup_L'));
indxFrontal_sup_L = find(source_atlas.tissue==x); 

x = find(ismember(atlas.tissuelabel,'Frontal_Sup_R'));
indxFrontal_sup_R = find(source_atlas.tissue==x); 

x = find(ismember(atlas.tissuelabel,'Frontal_Sup_Orb_L'));
indxFrontal_Sup_Orb_L = find(source_atlas.tissue==x);

x = find(ismember(atlas.tissuelabel,'Frontal_Sup_Orb_R'));
indxFrontal_Sup_Orb_L = find(source_atlas.tissue==x);

x=find(ismember(atlas.tissuelabel,'Heschl_R'));
indxHGR = find(source_atlas.tissue==x); 
 
x=find(ismember(atlas.tissuelabel,'Cingulum_Mid_L'));
indxCML = find(source_atlas.tissue==x); 

% Next, we normalise the individual MRI to derive parameters allowing to
% convert the mni- coordinates of the desired parcels into individual coordinates. For this we use ft_warp_apply.

template_grid=ft_convert_units(template_grid,'cm');% ensure no unit mismatch

norm_mri=ft_volumenormalise([],mri_realigned);
norm_mri_cm=ft_convert_units(norm_mri, 'cm');

% cfg = []; 
% cfg.interpmethod = 'nearest'; 
% cfg.parameter = 'anatomy'; 
% norm_mri_cm_interp = ft_sourceinterpolate(cfg, norm_mri_cm, source_avg);
% alternativ
cfg=[];
cfg.downsample = 10; %integer number (default = 1, i.e. no downsampling)
[norm_mri_cm_down] = ft_volumedownsample(cfg, norm_mri_cm)
cfg = []; 
cfg.interpmethod = 'nearest'; 
cfg.parameter = 'all'; 
norm_mri_cm_down_interp = ft_sourceinterpolate(cfg, norm_mri_cm_down, source_avg); % field initial fehlt
norm_mri_cm_down_interp.initial = norm_mri_cm_down.initial;
norm_mri_cm_down_interp.params = norm_mri_cm_down.params; % unklar ob das korrekt ist

posCML=template_grid.pos(indxCML,:);% xyz positions in mni coordinates
posHGL=template_grid.pos(indxHGL,:);% xyz positions in mni coordinates
posHGR=template_grid.pos(indxHGR,:);% xyz positions in mni coordinates

posFrontal_Sup_R=template_grid.pos(indxFrontal_sup_R,:);% xyz positions in mni coordinates
posFrontal_Sup_L=template_grid.pos(indxFrontal_sup_L,:);% xyz positions in mni coordinates

posback=ft_warp_apply(norm_mri_cm_down_interp.params,posFrontal_Sup_R,'sn2individual');
btiposindxFrontal_Sup_R= ft_warp_apply(pinv(norm_mri_cm_down_interp.initial),posback);% xyz positions in individual coordinates
posback=ft_warp_apply(norm_mri_cm_down_interp.params,posFrontal_Sup_L,'sn2individual');
btiposindxFrontal_Sup_L= ft_warp_apply(pinv(norm_mri_cm_down_interp.initial),posback);% xyz positions in individual coordinates



posback=ft_warp_apply(norm_mri_cm_down_interp.params,posCML,'sn2individual');
btiposCML= ft_warp_apply(pinv(norm_mri_cm_down.initial),posback);% xyz positions in individual coordinates

posback=ft_warp_apply(norm_mri_cm_down.params,posHGR,'sn2individual');
btiposHGR= ft_warp_apply(pinv(norm_mri_cm_down.initial),posback);% xyz positions in individual coordinates
posback=ft_warp_apply(norm_mri_cm_down.params,posHGL,'sn2individual');
btiposHGL= ft_warp_apply(pinv(norm_mri_cm_down.initial),posback);% xyz positions in individual coordinates
 
% Now we create a source model for these particular locations only.

cfg=[];
cfg.vol=hdm_ind_cm;
cfg.channel=data_bpfreq.label;  
cfg.grid.pos=[btiposCML;btiposHGL;btiposHGR;btiposindxFrontal_Sup_L;btiposindxFrontal_Sup_R];% units of cm (man kann alle Rois mit semicolon trennen)
cfg.grad=data_bpfreq.grad;
sourcemodel_virt=ft_prepare_leadfield(cfg);
figure
plot3(sourcemodel_virt.pos(:,1),sourcemodel_virt.pos(:,2),sourcemodel_virt.pos(:,3))
% And repeat the source analysis steps for above but now for 3 parcels represented in a total of 21 locations.


% cfg.grid.pos=[btiposCML;btiposHGL;btiposHGR;btiposindxFrontal_Sup_L;btiposindxFrontal_Sup_R];%
% plotten

cfg = [];
cfg.channel=data_bpfreq.label;
cfg.covariance='yes';
cfg.covariancewindow=[-0.4 1]; %?
cfg.keeptrials = 'yes';
avg = ft_timelockanalysis(cfg,data_bpfreq);
 
%% perform source analysis
cfg=[];
cfg.method='lcmv';
cfg.grid = sourcemodel_virt;
cfg.vol=hdm_ind_cm;
cfg.lcmv.keepfilter='yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.fixedori='no';
cfg.lcmv.lamda='5%';
source=ft_sourceanalysis(cfg, avg);

% On the basis of the computed filters, kept in the output, it is now possible to multiply them 
% with the data. This operation will yield time series commonly known as virtual sensors.

spatialfilter=cat(1,source.avg.filter{:});
virtsens=[];
for i=1:length(data_bpfreq.trial)
    virtsens.trial{i}=spatialfilter*data_bpfreq.trial{i};
 
end;
virtsens.time=data_bpfreq.time;
virtsens.fsample=data_bpfreq.fsample;
indx=[indxCML;indxHGL;indxHGR];
for i=1:length(virtsens.trial{1}(:,1))
    virtsens.label{i}=[num2str(i)];
end;

% Since our main interest is the time courses common to a given parcel we can average over within parcel locations.

cfg = [];
cfg.channel = virtsens.label(1:16);% cingulum is prepresented by 16 locations
cfg.avgoverchan = 'yes';
virtsensCML = ft_selectdata(cfg,virtsens);
virtsensCML.label = {'CML'};
 
cfg.channel = virtsens.label(17:19); % left heschl by 3
virtsensHGL = ft_selectdata(cfg,virtsens);
virtsensHGL.label = {'HGL'};
 
cfg.channel = virtsens.label(20:21); % right heschl by 2
virtsensHGR = ft_selectdata(cfg,virtsens);
virtsensHGR.label = {'HGR'};
 
%% append the data
virtsensparcel=ft_appenddata([],virtsensCML,virtsensHGL,virtsensHGR);








%%


    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = [-.5 1]; % oder von 0 bis 1?
    cfg.keeptrials = 'yes';
    avg_data_bpfreq = ft_timelockanalysis(cfg, data_bpfreq);
    
    % Now we make a first call to ft_sourceanalysis in order to compute the spatial filters on the basis of the entire data and keep them in the output for a later use.
    cfg = [];
    cfg.method = 'lcmv';
    cfg.grid = sourcemodel_lf;  % Stefan: cfg.grid = sourcemodel
    cfg.vol = hdm_cm; % Stefan: cfg.headmodel = vol;
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectnoise = 'yes';
    cfg.channel = data_bpfreq.label;
    cfg.lcmv.fixedori='no';
    cfg.lcmv.lamda='5%';
%     cfg.lcmv.weightnorm = 'nai'; % nicht klar ob redundant, da cfg.normalize = 'yes' bei ft_prepare_leadfield
%     cfg.lcmv.realfilter   = 'yes'; % keep imaginary part, macht kein Unterschied 
    try
        source_avg = ft_sourceanalysis(cfg, avg_data_bpfreq);
    catch
    end
%     The grid of estimated power values can be plotted superimposed on the anatomical MRI.
%     This requires the output of ft_sourceanalysis to match position of the template MRI to
%     which we warped the sourcemodel. Because of this warping that we already did, we can simply overwrite the grid position information without any further mathematical operation:
%  
%     outPath_extdisc_subj = [ outPath_extdisc 'MEG\sourcespace\ROI\run' run filesep];
%     if ~exist(outPath_extdisc_subj, 'dir')
%         mkdir (outPath_extdisc_subj)
%     end
    
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
  template_mri = ft_read_mri('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\anatomy\single_subj_T1.nii');
   template_mri_cm = ft_convert_units(template_mri,'cm')
 mri_realigned_cm = ft_convert_units(mri_realigned,'cm')
 
 cfg=[];
cfg.template ='\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\anatomy\single_subj_T1.nii';
mri_realigned_norm=ft_volumenormalise(cfg, mri_realigned)
mri_realigned_norm_cm =     ft_convert_units(mri_realigned_norm,'cm')



    cfg              = [];
    cfg.parameter    = 'avg.pow';
    cfg.interpmethod = 'nearest';
    source_avg_int  = ft_sourceinterpolate(cfg, source_avg, template_mri_cm);
    %passt nicht:
cfg               = [];
cfg.method        = 'slice';
cfg.funparameter  = 'avg.pow';
% cfg.maskparameter = cfg.funparameter;
% cfg.funcolorlim   = [0.0 1.2];
% cfg.opacitylim    = [0.0 1.2]; 
cfg.opacitymap    = 'rampup';  
ft_sourceplot(cfg,source_avg_int);

cfg=[];
cfg.template ='template_mri_cm';
source_avg_norm=ft_volumenormalise(cfg, source_avg) % wie kann ich source normalisieren?

    %passt nicht:
    figure
    plot3(source_avg.pos(source_avg.inside,1),source_avg.pos(source_avg.inside,2),source_avg.pos(source_avg.inside,3))
    hold on
    plot3(template_sourcemodel.sourcemodel.pos(template_sourcemodel.sourcemodel.inside,1),template_sourcemodel.sourcemodel.pos(template_sourcemodel.sourcemodel.inside,2),template_sourcemodel.sourcemodel.pos(template_sourcemodel.sourcemodel.inside,3))
    %     hold on
%     plot3(source_avg.avg.pow)
%     slice(source_avg.pos(source_avg.inside,1), source_avg.pos(source_avg.inside,2), source_avg.pos(source_avg.inside,3), source_avg.avg.pow(source_avg.inside))
   % To find your ROI, you use ft_volumelookup to get a binary matrix which tells you, which voxels of the standard MRI are in your ROI.
   
   template_mri = ft_read_mri('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\anatomy\single_subj_T1.nii');
   template_mri_cm = ft_convert_units(template_mri,'cm')
   atlas = ft_read_atlas('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\atlas\aal\ROI_MNI_V4.nii');
   atlas = ft_convert_units(atlas,'cm'); % assure that atlas and template_grid are expressed in the %same unitscfg = [];
   cfg = [];
   cfg.inputcoord = 'mni';
   cfg.atlas = atlas;
   cfg.roi = atlas.tissuelabel{1,1};
   mask = ft_volumelookup(cfg, template_mri_cm);
   

% create temporary mask according to the atlas entries
    tmp                  = repmat(template_sourcemodel.sourcemodel.inside,1,1);
    tmp(tmp==1)          = 0;
    tmp(mask)            = 1;
    
    % define inside locations according to the atlas based mask
    template_sourcemodel.sourcemodel.inside = tmp;
    
    % plot the atlas based grid
    figure;
    ft_plot_mesh(template_sourcemodel.sourcemodel.pos);
     figure;
     
    ft_plot_mesh(template_sourcemodel.sourcemodel.pos(template_sourcemodel.sourcemodel.inside,:));
    view ([-90 -90 -90])
   
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
     


