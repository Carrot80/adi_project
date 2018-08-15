
function adi_source_SVM (path2vol, path2data, mriPath, outPath_extdisc, freqbandname, condition)

% calculations spatial filter per Run and Subject, 
% time =  latency(1,:);

%% run 1:
if ~exist([outPath_extdisc '\run1\spatialfilter_loose_orientation_singletrials_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath,  num2str(1), freqbandname, condition, outPath_extdisc)    
end

%% run2:

if ~exist([outPath_extdisc '\run2\spatialfilter_loose_orientation_singletrials_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(2), freqbandname, condition, outPath_extdisc)
end

%% run3:

if ~exist([outPath_extdisc '\run3\spatialfilter_loose_orientation_singletrials_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(3), freqbandname, condition, outPath_extdisc)
end

end

function [virtsens] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, run, freq, condition, outPath_extdisc)
    
% load data:
index = zeros(3,1);
condition = cell(3,1);
if exist([path2data, 'like' '500_' run '.mat'], 'file')
%    condition{1} = 'data_like_bpfreq';
   index(1) = 1;
   data_like = load ([path2data, 'like' '500_' run '.mat'], 'cleanMEG_interp');
   [data_like_bpfreq] = adi_bpfilter(data_like.cleanMEG_interp, freq);
   data.like = data_like_bpfreq;
   clear data_like data_like_bpfreq
else
%     condition{1} = [];
end
if exist([path2data, 'dislike' '500_' run '.mat'], 'file')
    index(2) = 1;
%     condition{2} = 'data_dislike_bpfreq';
    data_dislike = load ([path2data, 'dislike' '500_' run '.mat'], 'cleanMEG_interp');
    [data_dislike_bpfreq] = adi_bpfilter(data_dislike.cleanMEG_interp, freq);
    data.dislike = data_dislike_bpfreq;
    clear data_dislike clear data_dislike_bpfreq
else
%      condition{2} = [];
end
if exist([path2data, 'dontcare' '500_' run '.mat'], 'file')
    index(3) = 1;
%     condition{3} = 'data_dontcare_bpfreq';
    data_dontcare = load ([path2data, 'dontcare' '500_' run '.mat'], 'cleanMEG_interp');
    [data_dontcare_bpfreq] = adi_bpfilter(data_dontcare.cleanMEG_interp, freq);
    data.dontcare = data_dontcare_bpfreq;
    clear data_dontcare data_dontcare_bpfreq
else
%        condition{3} = [];
end
fieldsnames = fields(data);
% ind_cond = find(~cellfun(@isempty,condition));
cfg = [];
switch length( find(index))
    case 3
    data_appended = ft_appenddata(cfg, data.(fieldsnames{1}), data.(fieldsnames{2}), data.(fieldsnames{3}));
    case 2 
    data_appended = ft_appenddata(cfg, data.(fieldsnames{1}), data.(fieldsnames{2}));
end   

    %% load template vol:
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
%         cfg            = [];
%         cfg.resolution = 1;
%         cfg.dim        = [256 256 256];
%         mri_resliced = ft_volumereslice(cfg, mri_realigned);
        mri_realigned_cm = ft_convert_units(mri_realigned, 'cm');
        grad = ft_convert_units(data_appended.grad, 'cm');
        data_bpfreq.grad = grad;

        cfg                = [];
        cfg.grid.warpmni   = 'yes';
        cfg.grid.template  = template_grid;
        cfg.grid.nonlinear = 'yes';
        cfg.mri            = mri_realigned_cm; % nicht sicher
        sourcemodel        = ft_prepare_sourcemodel(cfg);
    end
    if ~exist('hdm_ind', 'var')
        hdm_ind = load ([path2vol, 'vol.mat']);
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
    
    ft_plot_sens(data_bpfreq.grad,'style','g*', 'coil', 1);% plot the sensor array
    
    view ([0 -90 0])
%     
%     outPath = ([FigOutPath 'run' run '\SVM_result\' ]);
%     if ~exist(outPath, 'file')
%         mkdir (outPath)
%     end
%     
%     savefig([FigOutPath 'run' run '\SVM_result\check_vol_sourcemodel_grad.fig']);
%     fig = ([FigOutPath 'run' run '\SVM_result\check_vol_sourcemodel_grad']);
%     print('-dpng', fig); 
%     close all


    %We first create the leadfield using ft_prepare_leadfield using the individual head model from the previous step, the sensor array and the sourcemodel.
    cfg                 = [];
    cfg.channel         = data_appended.label;% ensure that rejected sensors are not present
    cfg.grad            = data_appended.grad;
    cfg.vol             = hdm_ind;
    cfg.lcmv.reducerank = 2; % default for MEG is 2, for EEG is 3
    cfg.grid = sourcemodel;
    % cfg.grid.pos=rois; 
    [leadfield] = ft_prepare_leadfield(cfg);% =grid

    %%
    
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = [0 1]; % oder von 0 bis 1?
    cfg.vartrllength = 2;
    cfg.keeptrials = 'yes';
    try
        avg_data_appended = ft_timelockanalysis(cfg, data_appended);
    catch ME
        if 1==strcmp(ME.message,'data has variable trial lengths, you specified not to accept that')
           for j = 1:length(data_appended.time)
            [ind(j)] =  isequal(length(data_appended.time{j}), 3052);
           end
           trl_ind = find(ind == false);
           for j=1:length(trl_ind)
               data_appended.time{1, trl_ind(j)}(:,1)=[];
               data_appended.trial{1,trl_ind(j)}(:,1)=[];
           end
           avg_data_appended = ft_timelockanalysis(cfg, data_appended);    
        end
        
    end
    
    % Now we make a first call to ft_sourceanalysis in order to compute the spatial filters on the basis of the entire data and keep them in the output for a later use.
    cfg = [];
    cfg.method = 'lcmv';
    cfg.grid = leadfield;  % Stefan: cfg.grid = sourcemodel
    cfg.vol = hdm_ind; % Stefan: cfg.headmodel = vol;
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectnoise = 'yes';
    cfg.channel = data_appended.label;
    cfg.lcmv.fixedori='no';
    cfg.lcmv.lamda='5%';
    try
        source_avg = ft_sourceanalysis(cfg, avg_data_appended);
    catch

    end

    outPath_extdisc_subj = [ outPath_extdisc 'run' run filesep];
    if ~exist(outPath_extdisc_subj, 'dir')
        mkdir (outPath_extdisc_subj)
    end
    
    spatialfilter_orig = source_avg.avg.filter;
    save([outPath_extdisc_subj 'spatialfilter_loose_orientation_singletrials_' freq], 'spatialfilter_orig')
    
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
    data_bpfreq.trial{1,k} = data_bpfreq.trial{1,k}(1:248,:)
end

for k = 1:length(data_bpfreq.trial)
    data_bpfreq.grad.label(249:end) = [];
    data_bpfreq.grad.chanori(249:end, :) = [];
    data_bpfreq.grad.chanpos(249:end, :) = [];
    data_bpfreq.grad.tra(249:end, :) = [];
    data_bpfreq.label(249:end) = [];
end

    
cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);


end
     


