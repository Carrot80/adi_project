function [] = adi_prepare_hdm(mniPath, path2data, outPath, subj)


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


if ~exist([mniPath 'mri_realigned.mat'], 'file')
    cfg = [];
    cfg.method = 'interactive';
    cfg.coordsys = 'ctf';
    mri = ft_read_mri([mniPath 'T1warped.mat']);
    mri_realigned = ft_volumerealign(cfg, mri);
    save ([mniPath 'mri_realigned'], 'mri_realigned');
else 
    load ([mniPath 'mri_realigned.mat'], 'mri_realigned')
end
% Segmentation works properly when the voxels of the anatomical images are homogenous (i.e. the size of the voxel is the same into each direction). 
cfg = [];
mri_reslice = ft_volumereslice(cfg, mri_realigned);

%% segment
% if ~exist([mniPath 'segmentedmri_EEG'], 'file')
    cfg           = [];
    cfg.output = {'brain','skull','scalp'} ;% für EEG relevant
%     cfg.output = {'brain','skull','scalp'};% für EE
    cfg.spmversion = 'spm12';
    segmentedmri  = ft_volumesegment(cfg, mri_reslice);
    segmentedmri.anatomy = mri_reslice.anatomy;
    save ([mniPath 'segmentedmri_EEG'], 'segmentedmri');
% else
%     load ([mniPath 'segmentedmri_EEG'], 'segmentedmri');
% end

%% Mesh: surfaces are created at the borders of the different tissue-types
if ~exist([mniPath 'bnd_EEG'], 'file')
    cfg = [];
    cfg.tissue = {'csf'};
    cfg.tissue = {'brain','skull','scalp'};
    cfg.numvertices = [3000 2000 1000];
    bnd = ft_prepare_mesh(cfg, segmentedmri);
    save ([outPath 'bnd_EEG'], 'bnd');
else
    load ([outPath 'bnd_EEG'], 'bnd');
end
disp(bnd(1))

%% Create a volume conduction model using 'openmeeg', or 'bemcp'.
if exist ([outPath 'vol.mat'], 'file')
    load([outPath 'vol.mat'], 'vol');
else
    cfg        = [];
    cfg.method = 'bemcp'; % You can also specify 'dipoli', 'bemcp', or another method.
    vol        = ft_prepare_headmodel(cfg, bnd);
    save ([outPath 'vol_EEG'], 'vol');
    disp(vol)
end

% figure;
% ft_plot_mesh(vol.bnd(3),'facecolor','none'); %scalp
% figure;
% ft_plot_mesh(vol.bnd(2),'facecolor','none'); %skull
% figure;
% ft_plot_mesh(vol.bnd(1),'facecolor','none'); %brain

figure
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05); % brain
hold on;
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4); % skull
hold on;
ft_plot_mesh(vol.bnd(3),'edgecolor','none','facecolor',[0.5 0.5 0.5], 'facealpha',0.4); % scalp


ftpath   = '\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\';
load(fullfile([ftpath 'template/sourcemodel/standard_sourcemodel3d10mm']));
template_grid = sourcemodel;
clear sourcemodel;


% source model
fileName = [outPath 'sourcemodel.mat'];
if ~exist(fileName, 'file')
    sourcemodel = conn_prepare_sourcemodel(segmentedmri, vol, 1, 0, 1);
    save(fileName, 'sourcemodel');
% else
%     load(filename);
% end


end


function sourcemodel = conn_prepare_sourcemodel(mri, vol, useTemplate, useMask, grayMatterOnly, templateResolution)

if nargin < 3
    useTemplate = 0;
end

if nargin < 4
    useMask = 0;
end

if nargin < 5
    grayMatterOnly = 0;
end

if nargin < 6
    templateResolution = 10;
end

if useTemplate
    %% Version using standard template model
    ftpath = fileparts(which('ft_defaults'));
    load(fullfile(ftpath, ['template/sourcemodel/standard_sourcemodel3d' num2str(templateResolution) 'mm']));
    template_grid = sourcemodel;
    clear sourcemodel;

    template_grid = ft_convert_units(template_grid, 'cm');
    
    if useMask
        atlas = ft_read_atlas(fullfile(ftpath, 'template\atlas\aal\ROI_MNI_V4.nii'));
        atlas = ft_convert_units(atlas,'cm');% assure that atlas and template_grid are expressed in the %same units
    
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
    end
    
    cfg                = [];
    cfg.grid.warpmni   = 'yes';
    cfg.grid.template  = template_grid;
    cfg.grid.nonlinear = 'yes'; % use non-linear normalization
    cfg.mri            = mri;
    sourcemodel        = ft_prepare_sourcemodel(cfg);
    sourcemodel.templateResolution = templateResolution;
else
    vol = ft_convert_units(vol, 'cm');
    
    %% compute the source model
    cfg = [];
    cfg.grid.resolution = 1;
    cfg.grid.tight  = 'yes';
    cfg.inwardshift = 0;
    cfg.headmodel   = vol;
    sourcemodel = ft_prepare_sourcemodel(cfg);
end

% Remove all source locations outside of gray matter if specified
if grayMatterOnly
    disp('Constraining to segmented gray matter...');
    mri = ft_convert_units(mri, 'cm');
    pos = ft_warp_apply(inv(mri.transform), sourcemodel.pos);
    pos = round(pos);
    isgray = false(1, size(pos, 2));
    for n=1:size(pos,1)
        if pos(n, 1) <=0 || pos(n, 2)<= 0 || pos(n, 3) <= 0 ...
                || pos(n, 1) > size(mri.gray, 1) || pos(n, 2) > size(mri.gray, 2) || pos(n, 3) > size(mri.gray, 3)
            isgray(n) = false;
        else
            isgray(n) = mri.gray(pos(n,1), pos(n,2), pos(n, 3)); 
        end
    end
    sourcemodel.inside = isgray;
    fprintf('%d locations inside gray matter.\n', sum(sourcemodel.inside));
end
