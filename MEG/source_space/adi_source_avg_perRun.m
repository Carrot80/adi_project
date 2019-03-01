
function adi_source_SVM (path2vol, path2data, mriPath, outPath_extdisc, freqbandname, condition)

% calculations spatial filter per Run and Subject, 
% time =  latency(1,:);

%% run 1:
if ~exist([outPath_extdisc '\run1\source_avg' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath,  num2str(1), freqbandname, condition, outPath_extdisc)    
end

%% run2:

if ~exist([outPath_extdisc '\run2\source_avg' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(2), freqbandname, condition, outPath_extdisc)
end

%% run3:

if ~exist([outPath_extdisc '\run3\source_avg' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(3), freqbandname, condition, outPath_extdisc)
end

end

function [virtsens] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, run, freq, condition, outPath_extdisc)
    
% load data:
index = zeros(3,1);
condition = cell(3,1);
if exist([path2data, 'Neu_Like' '500_' run '.mat'], 'file')
%    condition{1} = 'data_like_bpfreq';
   index(1) = 1;
   data.like = load ([path2data, 'Neu_Like' '500_' run '.mat'], 'cleanMEG_interp');
   [like_bpfreq] = adi_bpfilter(data.like.cleanMEG_interp, freq);

else
    condition{1} = [];
end
if exist([path2data, 'Neu_Dislike' '500_' run '.mat'], 'file')
    index(2) = 1;
    condition{2} = 'data_dislike_bpfreq';
    data.dislike = load ([path2data, 'Neu_Dislike' '500_' run '.mat'], 'cleanMEG_interp');
    [dislike_bpfreq] = adi_bpfilter(data.dislike.cleanMEG_interp, freq);

else
     condition{2} = [];
end

if ~exist('data', 'var')
    return
end

fieldsnames = fields(data);
% ind_cond = find(~cellfun(@isempty,condition));
cfg = [];
switch length( find(index))
    case 3
    data_appended = ft_appenddata(cfg, data.(fieldsnames{1}).cleanMEG_interp, data.(fieldsnames{2}).cleanMEG_interp, data.(fieldsnames{3}).cleanMEG_interp);
    case 2 
    data_appended = ft_appenddata(cfg, dislike_bpfreq, like_bpfreq);
end   

fsample = like_bpfreq.fsample;

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
        data_appended.grad = grad; % geändert

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
    
    ft_plot_sens(data_appended.grad,'style','g*', 'coil', 1);% plot the sensor array
    
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
    %     [data_dislike_bpfreq] = adi_bpfilter(data_dislike.cleanMEG_interp, freq);
    cfg = [];
    cfg.covariance = 'yes'; % compute the covariance for single trials, then average
    cfg.covariancewindow = [0 1]; % oder von 0 bis 1?  da ich mich für like vs dislike interessiere, evtl. Zeitbereich 0-1s, bei prä-vs.poststim-intervall-Konstrast prästim einbeziehen
    cfg.removemean = 'no'; % andernfalls zieht er von den Daten den Mittelwert ab zur Berechnung der Kovarianz
    cfg.vartrllength = 2;
    cfg.keeptrials = 'yes';
%     cfg.preproc.bpfilter = 'yes';
%     cfg.preproc.bpfreq = [1 95];
%     cfg.preproc.demean = 'no';
%     cfg.preproc.baselinewindow = [-0.5 0]; % baseline correction
%     avg_data_appended = ft_timelockanalysis(cfg, data_appended);
    avg_dislike_bpfreq = ft_timelockanalysis(cfg, dislike_bpfreq);
    avg_like_bpfreq = ft_timelockanalysis(cfg, like_bpfreq);
    
    avg_cov_dislike_cov = squeeze(mean(avg_dislike_bpfreq.cov));
    avg_cov_like_cov = squeeze(mean(avg_like_bpfreq.cov));
    
%     cfg = [];
%     cfg.vartrllength = 2;
%     cfg.keeptrials = 'yes';
%     cfg.covariance = 'no';
%     avg_data_like = ft_timelockanalysis(cfg, like_bpfreq);
%     avg_data_dislike = ft_timelockanalysis(cfg, dislike_bpfreq);
%     avg_data_like.cov = avg_data_appended.cov;
%     avg_data_dislike.cov = avg_data_appended.cov;
    

    [FourRef,Fref] = fftBasic(squeeze(mean(avg_like_bpfreq.trial)), round(fsample));
    figure;
    plot(FourRef, Fref)
    outPath_extdisc_subj = [ outPath_extdisc 'run' run filesep];
    if ~exist(outPath_extdisc_subj, 'dir')
        mkdir (outPath_extdisc_subj)
    end
    savefig([outPath_extdisc '\run' run '\freqspectrum_avg_filtered_trials.fig'])
    close
    
    cfg = [];
    cfg.latency = 'prestim';
    avg_baseline_like=ft_selectdata(cfg, avg_data_like);
    avg_baseline_dislike=ft_selectdata(cfg, avg_data_dislike);
    
    
    
    % Now we make a first call to ft_sourceanalysis in order to compute the spatial filters on the basis of the entire data and keep them in the output for a later use.
    cfg = [];
    cfg.method = 'lcmv';
%     cfg.latency = [-0.5 1];
    cfg.grid = leadfield;  % Stefan: cfg.grid = sourcemodel
    cfg.vol = hdm_ind; % Stefan: cfg.headmodel = vol;
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectnoise = 'yes';
    cfg.channel = data_appended.label;
    cfg.lcmv.fixedori='no';
    cfg.lcmv.reducerank = 'yes';
    cfg.lcmv.lambda='5%';
    cfg.latency = [0 0.8];
    [source_avg_dislike_bpfreq] = ft_sourceanalysis(cfg, avg_dislike_bpfreq );
    [source_avg_like_bpfreq] = ft_sourceanalysis(cfg, avg_like_bpfreq );
 
    source_avg.(['run' run]).like = source_avg_lik_bpfreq;
    source_avg.(['run' run]).dislike = source_avg_dislike_bpfreq;

    save([outPath_extdisc_subj 'source_avg_' freq], 'source_avg')
%     
    %% plot like
    nansum(source_avg_like_bpfreq.avg.pow)
    nansum(source_avg_dislike_bpfreq.avg.pow)
    source_like_NAI = source_avg_like_bpfreq;
    source_like_NAI.avg.pow = source_avg_like_bpfreq.avg.pow ./ source_avg_like_bpfreq.avg.noise;

    cfg = [];
%     cfg.downsample = 2;
    cfg.parameter = 'avg.pow';
    sourceNAIInt_like = ft_sourceinterpolate(cfg, source_like_NAI , mri_realigned_cm);
   
    cfg = [];
    cfg.method        = 'slice';
    cfg.funparameter  = 'avg.pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.funcolorlim   = [1.0 5];
    cfg.opacitylim    = [1.0 5];
    cfg.opacitymap    = 'rampup';  
    figure
    ft_sourceplot(cfg, sourceNAIInt_like);
    
    %% plot dislike:
    
    source_dislike_NAI = source_avg_dislike_bpfreq;
    source_dislike_NAI.avg.pow = source_avg_dislike_bpfreq.avg.pow ./ source_avg_dislike_bpfreq.avg.noise;

    cfg = [];
%     cfg.downsample = 2;
    cfg.parameter = 'avg.pow';
    sourceNAIInt_dislike = ft_sourceinterpolate(cfg, source_dislike_NAI , mri_realigned_cm);

    
    cfg = [];
    cfg.method        = 'slice';
    cfg.funparameter  = 'avg.pow';
    cfg.maskparameter = cfg.funparameter;
    cfg.funcolorlim   = [1.0 5];
    cfg.opacitylim    = [1.0 5];
    cfg.opacitymap    = 'rampup';  
    figure
    ft_sourceplot(cfg, sourceNAIInt_dislike);
    
    
    %% plot diff
    sourceDiff = source_avg_like_bpfreq;
    sourceDiff.avg.pow = (source_avg_like_bpfreq.avg.pow - source_avg_dislike_bpfreq.avg.pow) ./ source_avg_dislike_bpfreq.avg.pow;
    cfg            = [];
%     cfg.downsample = 2;
    cfg.parameter  = 'avg.pow';
    sourceDiffInt  = ft_sourceinterpolate(cfg, sourceDiff , mri_realigned_cm);
    
    cfg = [];
    cfg.method        = 'slice';
    cfg.funparameter  = 'avg.pow';
    cfg.maskparameter = cfg.funparameter;
%     cfg.funcolorlim   = [-1 2];
%     cfg.opacitylim    = [-1 2];
    cfg.opacitymap    = 'rampup';  
    ft_sourceplot(cfg, sourceDiffInt);

%     nansum(source_dislike_NAI.avg.pow)
%     nansum(source_like_NAI.avg.pow)
    
    [h,p,ci,tstat] = ttest(source_dislike_NAI.avg.pow,source_like_NAI.avg.pow,'Alpha',0.001)
    [p,h,stats] = ranksum(source_dislike_NAI.avg.pow,source_like_NAI.avg.pow)
    
      %% reduce the source reconstructed data to the dominant orientation
    cfg = [];
    cfg.projectmom = 'yes';
    cfg.keepnoisemom = 'yes';
    source_proj_like = ft_sourcedescriptives(cfg, source_avg_like); % mir ist unklar, worauf dies hinausläuft
    source_proj_dislike = ft_sourcedescriptives(cfg, source_avg_dislike); 
    
    source_proj_like.avg.nai=source_proj_like.avg.nai(source_proj_like.inside==1);
    source_proj_dislike.avg.nai=source_proj_dislike.avg.nai(source_proj_dislike.inside==1);
    
    cfg = [];
    cfg.parameter = 'avg.nai';
    cfg.operation = '((x2-x1)./x1)*100';
    source_ratio = ft_math(cfg, source_proj_like, source_proj_dislike);

    source_ratio.pos = template_grid.pos;
    templatefile = fullfile('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20181213\template\anatomy\single_subj_T1.nii');
    template_mri = ft_read_mri(templatefile);
    cfg              = [];
%     cfg.voxelcoord   = 'no';
    cfg.parameter    = 'avg.nai'; 
    cfg.interpmethod = 'linear';
    source_int  = ft_sourceinterpolate(cfg, source_proj_like, mri_realigned_cm);
    source_int.coordsys = 'mni';

    cfg              = [];
%     cfg.voxelcoord   = 'no';
    cfg.parameter    = 'nai'; 
    cfg.interpmethod = 'linear';
    source_int  = ft_sourceinterpolate(cfg, source_ratio, mri_realigned);
    source_int.coordsys = 'mni';

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
%    
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
    case 'bp1_45Hz'
        bpfreq = [1 45];
    case'1_5_45Hz'
        bpfreq = [1.5 45];
    case 'bp2_45Hz'
        bpfreq = [2 45];
    case 'bp3_45Hz'
        bpfreq = [3 45];
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
    case 'bp10-45Hz'
        bpfreq = [10 45];
end

cfg = [];
% cfg.channel = 'meg';
% cfg.covariance = 'yes'; % compute the covariance for single trials, then average
% cfg.preproc.bpfilter = 'yes';
% cfg.preproc.bpfreq = [5 75];
% cfg.preproc.demean = 'yes';
% cfg.preproc.baselinewindow = [-inf 0]; % baseline correction
cfg.keeptrials = 'no';
cfg.vartrllength = 2;
cfg.trials        = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') %|| 1 == strcmp(bpname, 'bp1-45Hz')
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
  
cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);
data_bpfreq_res_sel.trialinfo=data_bpfreq.trialinfo; 

fn_data_bpfreq_res_sel = fieldnames(data_bpfreq_res_sel);
fn_filename = fieldnames(filename);

diff_fieldnames = setdiff(fn_filename, fn_data_bpfreq_res_sel);

for k=1:length(diff_fieldnames)
    data_bpfreq_res_sel.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
end

end
     


