function [] = adi_source_reconstruction_stats(path2vol, path2data, mriPath, outPath, run, freq, condition1, condition2, condition3)

load ([path2data, condition1 '500_' run '_' freq '.mat'])

% load template vol:
load('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180305\template\headmodel\standard_singleshell');
cfg = [];
cfg.grid.xgrid  = -20:1:20;
cfg.grid.ygrid  = -20:1:20;
cfg.grid.zgrid  = -20:1:20;
cfg.grid.unit   = 'cm';
cfg.grid.tight  = 'yes';
cfg.inwardshift = -1.5;
cfg.headmodel  = vol;
template_grid  = ft_prepare_sourcemodel(cfg);

figure;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));
hold on
ft_plot_vol(vol,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;


%%
% Load atlas and create a binary mask
atlas = ft_read_atlas('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180305\template\atlas\aal\ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas,'cm');% assure that atlas and template_grid are expressed in the %same units
cfg = [];
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
view ([-90 -90 -90])
%%
% Inverse-warp the subject specific grid to the atlas based template grid
% For this step the individual volume is required:

load ([mriPath, 'mri_realigned.mat'])
cfg            = [];
cfg.resolution = 1;
cfg.dim        = [256 256 256];
mri_resliced = ft_volumereslice(cfg, mri_realigned)
mri_resliced = ft_convert_units(mri_resliced, 'cm')
grad = ft_convert_units(data_bpfreq.grad, 'cm');
data_bpfreq.grad = grad;

cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes';
cfg.mri            = mri_resliced; % nicht sicher
sourcemodel        = ft_prepare_sourcemodel(cfg);

hdm_ind = load ([path2vol, 'vol.mat'])
hdm_ind = ft_convert_units(hdm_ind.vol, 'cm'); 
% Plot the final source model together with the individual head model and
% the sensor array:
close all
% the 4D/bti system is expressed in units of 'm', therefore we force all geometrical objects to have the same unit
% vol_indiv = ft_convert_units(hdm, 'm');
% sourcemodel = ft_convert_units(sourcemodel, 'm');
 
figure; hold on     % plot all objects in one figure
 
ft_plot_vol(hdm_ind,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; %camlight;
alpha 0.4           % make the surface transparent
 
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));% plot only locations inside the volume
 
ft_plot_sens(data_bpfreq.grad,'style','*r');% plot the sensor array
view ([0 -90 0])

%We first create the leadfield using ft_prepare_leadfield using the individual head model from the previous step, the sensor array and the sourcemodel.
cfg                 = [];
cfg.channel         = data_bpfreq.label;% ensure that rejected sensors are not present
cfg.grad            = data_bpfreq.grad;
cfg.vol             = hdm_ind;
cfg.lcmv.reducerank = 2; % default for MEG is 2, for EEG is 3
cfg.grid = sourcemodel;
% cfg.grid.pos=rois; 
[grid] = ft_prepare_leadfield(cfg);

%% Select data
% We want to reconstruct the locations of the M100 component. In the present case this component is maximal around 50 to 180ms post stimulus onset. Therefore we segment the data around this interval using ft_redefinetrial. Furthermore we would contrast the source solution against a prestimulus baseline of equal length.
cfg = [];
cfg.toilim = [-0.18 0.02];
data_pre = ft_redefinetrial(cfg, data_bpfreq);
cfg.toilim = [0.02 0.18];
data_post = ft_redefinetrial(cfg, data_bpfreq);
% cfg = [];   
% dataAll = ft_appenddata(cfg, data_pre, data_post);
% dataAll.fsample = data_bpfreq.fsample;

cfg = [];
cfg.covariance = 'yes';
cfg.covariancewindow = [-.4 1]; % oder von 0 bis 1?
cfg.keeptrials = 'yes';
avg_data_bpfreq = ft_timelockanalysis(cfg, data_bpfreq);
 
cfg = [];
cfg.covariance='yes';
cfg.keeptrials = 'yes';
avg_data_pre = ft_timelockanalysis(cfg, data_pre);
avg_data_post = ft_timelockanalysis(cfg, data_post);

% Now we make a first call to ft_sourceanalysis in order to compute the spatial filters on the basis of the entire data and keep them in the output for a later use.
cfg = [];
cfg.method = 'lcmv';
cfg.grid = grid;
cfg.vol = hdm_ind;
cfg.lcmv.keepfilter = 'yes';
cfg.channel = data_bpfreq.label;
cfg.lcmv.fixedori='yes';
cfg.lcmv.lamda='5%';
source_avg_data_bpfreq = ft_sourceanalysis(cfg, avg_data_bpfreq);

% Subsequently we reconstruct the activity in the pre and post stimulus intervals using the precomputed filters.
cfg = [];
cfg.method = 'lcmv';
cfg.grid = grid;
cfg.grid.filter = source_avg_data_bpfreq.avg.filter;
cfg.vol = hdm_ind;
cfg.rawtrial = 'yes';
cfg.channel = data_bpfreq.label;
source_avg_datapre = ft_sourceanalysis(cfg, avg_data_pre);
source_avg_datapost = ft_sourceanalysis(cfg, avg_data_post);



cfg = [];
cfg.parameter = 'inside';
cfg.operation = '((x1-x2)./x2)*100';
sourceDiff = ft_math(cfg, source_avg_datapost, source_avg_datapre); % für SVM nehmen (mit vs multiplizieren?)


%% keep covariance in the output
% cfg = [];
% cfg.channel=data_bpfreq.label;
% cfg.covariance='yes';
% cfg.covariancewindow = [-0.5 0.5]; 
% avg = ft_timelockanalysis(cfg, data_bpfreq);
 
%% perform source analysis
% cfg=[];
% cfg.method='lcmv';
% cfg.grid = grid;
% cfg.vol=hdm_ind;
% cfg.lcmv.keepfilter='yes';
% cfg.lcmv.fixedori='yes';
% cfg.lcmv.lamda='5%';
% source=ft_sourceanalysis(cfg, avg);

% On the basis of the computed filters, kept in the output, it is now possible to multiply them with the data. This operation will yield time series commonly known as virtual sensors.

spatialfilter = cat(1,source_avg_data_bpfreq.avg.filter{:});
virtsens = [];
for i = 1:length(data_bpfreq.trial)
    virtsens.trial{i} = spatialfilter*data_bpfreq.trial{i};
 
end;
virtsens.time=data_bpfreq.time;
virtsens.fsample=data_bpfreq.fsample;
% indx=[indxFML;indxHGL;indxHGR];
for i=1:length(virtsens.trial{1}(:,1))
    virtsens.label{i}=[num2str(i)];
end;

figure
plot (virtsens.time{1, 1}, virtsens.trial{1,1}(1,:))

% Now statistical analysis can be performed.

cfg = [];
cfg.parameter    = 'pow';
cfg.dim          = grid.dim;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'max';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

ntrials = numel(source_avg_datapre.trial);
design  = zeros(2,2*ntrials);
design(1,1:ntrials) = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials) = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];
 
cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
stat = ft_sourcestatistics(cfg,source_avg_datapre,source_avg_datapost);
stat.pos = template_grid.pos;% keep positions for plotting later

templatefile = '\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180305\template\anatomy\single_subj_T1.nii';
template_mri = ft_read_mri(templatefile);

stat.inside = template_grid.inside;
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'stat';
cfg.interpmethod = 'nearest';
statint  = ft_sourceinterpolate(cfg, stat, template_mri); 
cfg.parameter    = 'mask'; 
maskint  = ft_sourceinterpolate(cfg, stat, template_mri); 
statint.mask = maskint.mask;

% Abbildung funktioniert nicht, da keine signifikanten Werte vorhanden sind:
statint.coordsys = 'mni';
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.maskparameter = 'mask';
atlas = ft_convert_units(atlas, 'mm')

cfg.atlas         = atlas;
cfg.location = 'max';
cfg.funcolorlim   = [-5 5];
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, statint); 

% We repeat the steps from above and plot the result in parceled brain space
cfg = [];
cfg.parameter    = 'mask';
statint2 = ft_sourceinterpolate(cfg, statint, atlas);
parcel = ft_sourceparcellate(cfg, statint2, atlas);
cfg.parameter    = 'mask';
maskint2 = ft_sourceinterpolate(cfg, maskint, atlas);
parcelmask = ft_sourceparcellate(cfg, maskint2, atlas);

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

% First, we interpolate the statistical result and the atlas.
cfg = []; 
cfg.interpmethod = 'nearest'; 
cfg.parameter = 'tissue'; 
stat_atlas = ft_sourceinterpolate(cfg, atlas, stat);

% Determine the index of the label of interest,

x = find(ismember(atlas.tissuelabel,'Heschl_L'));
% and determine the index points of locations within the desired parcel.

indxHGL = find(stat_atlas.tissue==x);
% These steps can be repeated for all desired parcels. In the present case the remaining two:

x=find(ismember(atlas.tissuelabel,'Heschl_R'));
indxHGR = find(stat_atlas.tissue==x); 
 
x=find(ismember(atlas.tissuelabel,'Cingulum_Mid_L'));
indxCML = find(stat_atlas.tissue==x); 

% Next, we normalise the individual MRI to derive parameters allowing to convert the mni- coordinates of the desired parcels into individual coordinates. For this we use ft_warp_apply.
template_grid=ft_convert_units(template_grid,'mm');% ensure no unit mismatch
norm=ft_volumenormalise([],mri_resliced);
 
 
posCML=template_grid.pos(indxCML,:);% xyz positions in mni coordinates
posHGL=template_grid.pos(indxHGL,:);% xyz positions in mni coordinates
posHGR=template_grid.pos(indxHGR,:);% xyz positions in mni coordinates
 
posback=ft_warp_apply(norm.params,posCML,'sn2individual');
btiposCML= ft_warp_apply(pinv(norm.initial),posback);% xyz positions in individual coordinates
 
posback=ft_warp_apply(norm.params,posHGL,'sn2individual');
btiposHGL= ft_warp_apply(pinv(norm.initial),posback);% xyz positions in individual coordinates
 
posback=ft_warp_apply(norm.params,posHGR,'sn2individual');
btiposHGR= ft_warp_apply(pinv(norm.initial),posback);% xyz positions in individual coordinates

% Now we create a source model for these particular locations only.

cfg=[];
cfg.vol=hdm;
cfg.channel=data_bpfreq.label;  
cfg.grid.pos=[btiposCML;btiposHGL;btiposHGR]./1000;% units of m
cfg.grad=data_bpfreq.grad;
sourcemodel_virt=ft_prepare_leadfield(cfg);

% And repeat the source analysis steps for above but now for 3 parcels represented in a total of 21 locations.

%% keep covariance in the output
cfg = [];
cfg.channel=data_bpfreq.label;
cfg.covariance='yes';
cfg.covariancewindow = [-0.5 0.5]; 
avg = ft_timelockanalysis(cfg, data_bpfreq);
 
%% perform source analysis
cfg=[];
cfg.method='lcmv';
cfg.grid = grid;
cfg.vol=hdm_ind;
cfg.lcmv.keepfilter='yes';
cfg.lcmv.fixedori='yes';
cfg.lcmv.lamda='5%';
source=ft_sourceanalysis(cfg, avg);

% On the basis of the computed filters, kept in the output, it is now possible to multiply them with the data. This operation will yield time series commonly known as virtual sensors.

spatialfilter = cat(1,source.avg.filter{:});
virtsens = [];
for i = 1:length(data_bpfreq.trial)
    virtsens.trial{i} = spatialfilter*data_bpfreq.trial{i};
 
end;
virtsens.time=data_bpfreq.time;
virtsens.fsample=data_bpfreq.fsample;
indx=[indxFML;indxHGL;indxHGR];
for i=1:length(virtsens.trial{1}(:,1))
    virtsens.label{i}=[num2str(i)];
end;

%Since our main interest is the time courses common to a given parcel we can average over within parcel locations.
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


%  now we compute the source wave forms, plot and evaluate the result.

cfg=[];
tlkvc=ft_timelockanalysis(cfg, virtsensparcel);
figure;
for i=1:length(tlkvc.label)
    cfg=[];
    cfg.channel = tlkvc.label{i};
    cfg.parameter = 'avg';
    cfg.xlim    = [-.1 1];
 
    subplot(2,2,i);ft_singleplotER(cfg,tlkvc);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Now we can subtract the two conditions, normalize by the power in the pre stimulus interval and multiply by 100. Thereby the data is expressed in percentage change from pre stimulus baseline.
cfg = [];
cfg.parameter = 'avg.pow';
cfg.operation = '((x1-x2)./x2)*100';
sourceDiff = ft_math(cfg, source_tdatapost, source_tdatapre);

% The result is then interpolated on the template mri after the individual dipole locations are set back to be equal with the template_grid locations. This is done with ft_sourceinterpolate.
templatefile = '\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180305\template\anatomy\single_subj_T1.nii';
template_mri = ft_read_mri(templatefile);
sourceDiff.pos=template_grid.pos;
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
sourceDiff_Int  = ft_sourceinterpolate(cfg, sourceDiff, template_mri);

% Finally, we can plot the result using ft_sourceplot:
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
% cfg.location = [64 -32 8]; % ?
cfg.funcolormap = 'jet';
ft_sourceplot(cfg, sourceDiff_Int);

templatefile = '\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180305\external\spm8\templates\T1.nii';
template_mri = ft_read_mri(templatefile);
template_mri = ft_convert_units(template_mri, 'cm');
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, sourceDiff, template_mri);
%%
cfg=[];
parcel = ft_sourceparcellate(cfg, source_int, atlas);

%% We create a dummy structure where we identify the power values per voxel and use this for subsequent plotting.

dummy = atlas;
for i=1:length(parcel.pow)
      dummy.tissue(find(dummy.tissue == i)) = parcel.pow(i);
end;
%% hier weitermachen:
source_int.parcel = dummy.tissue;
source_int.coordsys = 'mni';
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'parcel';
cfg.funcolormap    = 'jet';
cfg.renderer = 'zbuffer';
cfg.location = [-42 -20 6];
cfg.atlas = atlas;
cfg.funcolorlim = [-30 30];
ft_sourceplot(cfg, source_int); 

%% Alternatively, the maximal activity in the left Heschl gyrus can be plotted on the brain surface as follows.
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [-50 50];
cfg.funcolormap    = 'jet'; 
cfg.projmethod     = 'nearest'; 
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projthresh     = 0.6; % percent of maximum to present the hill of acitivty
cfg.camlight       = 'no';
ft_sourceplot(cfg, source_int);
view ([70 20 50])
light ('Position',[70 20 50])
material dull




%%

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




%%

% cfg              = [];
% cfg.method       = 'sam';
% cfg.channel       = 'MEG';
% % cfg.frequency    = [1 4];
% cfg.grid         = grid;
% cfg.headmodel    = vol;
% cfg.sam.projectnoise = 'yes';
% cfg.sam.lambda       = '5%';
% cfg.sam.keepfilter   = 'yes';
% cfg.sam.realfilter   = 'yes';
% cfg.keeptrials = 'yes';
% cfg.cov = tdataAll.cov;
% cfg.covariance =  tdataAll.cov;
% tsourceAll = ft_sourceanalysis(cfg, tdataAll, [])
% test:
% tsource_post = ft_sourceanalysis(cfg, tdata_post);
% tsource_pre = ft_sourceanalysis(cfg, tdata_pre);
% By placing this pre-computed filter inside cfg.grid.filter, it can now be applied to each condition separately.
cfg = [];
cfg.projectmom = 'yes';
cfg.keepnoisemom = 'yes';
[source] = ft_sourcedescriptives(cfg, sourceAll)

cfg.grid.filter = tsourceAll.avg.filter;
tsourcePre  = ft_sourceanalysis(cfg, tdata_pre);
tsourcePost = ft_sourceanalysis(cfg, tdata_post);
% save sourcePre_con sourcePre_con 
% save sourcePost_con sourcePost_con


sourceDiff = tsourcePost;
sourceDiff.avg.pow = (tsourcePost.avg.pow - tsourcePre.avg.pow) ./ tsourcePre.avg.pow;

cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'avg.pow';
sourceDiffInt  = ft_sourceinterpolate(cfg, sourceDiff , mri_resliced);
% Now plot the power ratios:
cfg = [];
cfg.method        = 'slice';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [-2 2];
cfg.opacitylim    = [-2 2]; 
cfg.opacitymap    = 'rampup';  
ft_sourceplot(cfg, sourceDiffInt);







cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [-2 2];
cfg.opacitylim    = [-1 6]; 
cfg.opacitymap    = 'rampup';  
ft_sourceplot(cfg, sourceDiffInt);


cfg = [];
cfg.nonlinear     = 'no';
sourceDiffIntNorm = ft_volumenormalise(cfg, sourceDiffInt);

cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'avg.pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolorlim    = [0.0 1.2];
cfg.funcolormap    = 'jet';
cfg.opacitylim     = [0.0 1.2]; 
cfg.opacitymap     = 'rampup';  
cfg.projmethod     = 'nearest'; 
cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
cfg.surfdownsample = 10;  % downsample to speed up processing
ft_sourceplot(cfg, sourceDiffIntNorm);
view ([90 0])             % rotate the object in the view

%% Analyze

    [results_delta_pre, source_delta_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 1:3);




filepath = fullfile(outpath, 'results_delta_post.mat');
if ~exist(filepath, 'file')
    [results_delta_post, source_conn_delta_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 1:3);
    save(filepath, 'results_delta_post', 'source_conn_delta_post');
    disp('Completed delta post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_theta_pre.mat');
if ~exist(filepath, 'file')
    [results_theta_pre, source_conn_theta_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 4:7);
    save(filepath, 'results_theta_pre', 'source_conn_theta_pre');
    disp('Completed theta pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_theta_post.mat');
if ~exist(filepath, 'file')
    [results_theta_post, source_conn_theta_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 4:7);
    save(filepath, 'results_theta_post', 'source_conn_theta_post');
    disp('Completed theta post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_alpha_pre.mat');
if ~exist(filepath, 'file')
    [results_alpha_pre, source_conn_alpha_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 8:15);
    save(filepath, 'results_alpha_pre', 'source_conn_alpha_pre');
    disp('Completed alpha pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_alpha_post.mat');
if ~exist(filepath, 'file')
    [results_alpha_post, source_conn_alpha_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 8:15);
    save(filepath, 'results_alpha_post', 'source_conn_alpha_post');
    disp('Completed alpha post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_gamma_pre.mat');
if ~exist(filepath, 'file')
    [results_gamma_pre, source_conn_gamma_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 25:50);
    save(filepath, 'results_gamma_pre', 'source_conn_gamma_pre');
    disp('Completed gamma pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_gamma_post.mat');
if ~exist(filepath, 'file')
    [results_gamma_post, source_conn_gamma_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 25:50);
    save(filepath, 'results_gamma_post', 'source_conn_gamma_post');
    disp('Completed gamma post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_highgamma_pre.mat');
if ~exist(filepath, 'file')
    [results_highgamma_pre, source_conn_highgamma_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 70:100);
    save(filepath, 'results_highgamma_pre', 'source_conn_highgamma_pre');
    disp('Completed highgamma pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_highgamma_post.mat');
if ~exist(filepath, 'file')
    [results_highgamma_post, source_conn_highgamma_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 70:100);
    save(filepath, 'results_highgamma_post', 'source_conn_highgamma_post');
    disp('Completed highgamma post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end







%% Stefan:
%function [source_int, source_proj] = visual_source_analysis(vol, sourcemodel, mri, data, foi, latency, tapsmofrq, trials)
 
cfg = [];
cfg.trials = trials;
selected = ft_selectdata(cfg, data);

source_proj = conn_source_reconstruction(vol, sourcemodel, selected, foi, latency, tapsmofrq);

cfg              = [];
cfg.parameter    = {'pow', 'nai'};
cfg.interpmethod = 'linear';
source_int  = ft_sourceinterpolate(cfg, source_proj, mri);


  
%% prepare leadfield

 vol = ft_convert_units(vol, 'cm');
    
    %% compute the source model
    cfg = [];
    cfg.grid.resolution = 1;
    cfg.grid.tight  = 'yes';
    cfg.inwardshift = 0;
    cfg.headmodel   = vol;
    cfg.normalize = 'yes';
    sourcemodel2 = ft_prepare_sourcemodel(cfg);











    filename = [outpath 'sourcemodel.mat'];
    if ~exist(filename, 'file')
        sourcemodel = conn_prepare_sourcemodel(segmentedmri, vol, 1, 0, 1);
        save(filename, 'sourcemodel');
    else
        load(filename);
    end

%%







end