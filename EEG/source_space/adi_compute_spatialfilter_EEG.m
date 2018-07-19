
function adi_spatialfilter_EEG (path2vol, path2data, mriPath, SubjOutPath, outPath_extdisc, freqbandname, like, dislike)

% calculations of spatial filter

%% run 1:
% if ~exist([SubjOutPath '\run1\SVM_result\virtSens\virtsens_like_vs_dislike_1_' freqbandname '.mat'], 'file')
  if ~exist([outPath_extdisc 'sourcespace\run1\spatialfilter_loose_orientation_singletrials_like_' freqbandname '.mat'], 'file')
      
    [virtsens_like, virtsens_ns_like] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, num2str(1), freqbandname, 'like', outPath_extdisc)                                                                   
    [virtsens_dislike, virtsens_ns_dislike] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, num2str(1), freqbandname, 'dislike', outPath_extdisc)

%     [Condition1vs2] = adi_crossvalidation (virtsens_like, virtsens_dislike, freqbandname, 'virtsens_like_vs_dislike', latency, num2str(1), SubjOutPath)
%     adi_figureTPRcrossval_SVM (Condition1vs2, time, like, dislike, 'virtsens_like_vs_dislike', SubjOutPath, freqbandname, num2str(1))

%     [Condition1vs2] = adi_crossvalidation (virtsens_ns_like, virtsens_ns_dislike, freqbandname, 'virtsens_ns_like_vs_dislike', latency, num2str(1), SubjOutPath)
%     adi_figureTPRcrossval_SVM (Condition1vs2, time, like, dislike, 'virtsens_ns_like_vs_dislike', SubjOutPath, freqbandname, num2str(1))

    clear virtsens_like virtsens_ns_like virtsens_dislike virtsens_ns_dislike

end

%% run2:

% if ~exist([SubjOutPath '\run2\SVM_result\virtSens\virtsens_like_vs_dislike_2_' freqbandname '.mat'], 'file')
if ~exist([outPath_extdisc 'sourcespace\run2\spatialfilter_loose_orientation_singletrials_like_' freqbandname '.mat'], 'file')
    
    [virtsens_like, virtsens_ns_like] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, num2str(2), freqbandname, 'like', outPath_extdisc)
    [virtsens_dislike, virtsens_ns_dislike] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, num2str(2), freqbandname, 'dislike', outPath_extdisc)

%     [Condition1vs2] = adi_crossvalidation (virtsens_like, virtsens_dislike, freqbandname, 'virtsens_like_vs_dislike', latency, num2str(2), SubjOutPath)
%     adi_figureTPRcrossval_SVM (Condition1vs2, time, like, dislike, 'virtsens_like_vs_dislike', SubjOutPath, freqbandname, num2str(2))
% 
%     [Condition1vs2] = adi_crossvalidation (virtsens_ns_like, virtsens_ns_dislike, freqbandname, 'virtsens_ns_like_vs_dislike', latency, num2str(2), SubjOutPath)
%     adi_figureTPRcrossval_SVM (Condition1vs2, time, like, dislike, 'virtsens_ns_like_vs_dislike', SubjOutPath, freqbandname, num2str(2))

    clear virtsens_like virtsens_ns_like virtsens_dislike virtsens_ns_dislike
end

%% run3:

% if ~exist([SubjOutPath '\run3\SVM_result\virtSens\virtsens_like_vs_dislike_3_' freqbandname '.mat'], 'file')
  if ~exist([outPath_extdisc 'sourcespace\run3\spatialfilter_loose_orientation_singletrials_like_' freqbandname '.mat'], 'file')
    [virtsens_like, virtsens_ns_like] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, num2str(3), freqbandname, 'like', outPath_extdisc)
    [virtsens_dislike, virtsens_ns_dislike] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, num2str(3), freqbandname, 'dislike', outPath_extdisc)

%     [Condition1vs2] = adi_crossvalidation (virtsens_like, virtsens_dislike, freqbandname, 'virtsens_like_vs_dislike', latency, num2str(3), SubjOutPath)
%     adi_figureTPRcrossval_SVM (Condition1vs2, time, like, dislike, 'virtsens_ns_like_vs_dislike', SubjOutPath, freqbandname, num2str(3))
% 
%     [Condition1vs2] = adi_crossvalidation (virtsens_ns_like, virtsens_ns_dislike, freqbandname, 'virtsens_ns_like_vs_dislike', latency, num2str(3), SubjOutPath)
%     adi_figureTPRcrossval_SVM (Condition1vs2, time, like, dislike, 'virtsens_ns_like_vs_dislike', SubjOutPath, freqbandname, num2str(3))

    clear virtsens_like virtsens_ns_like virtsens_dislike virtsens_ns_dislike

end

end

function [virtsens, virtsens_ns] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, SubjOutPath, run, freq, condition, outPath_extdisc)
    
% load data:
    if exist([path2data, condition '500_' run '.mat'], 'file')
        load ([path2data, condition '500_' run '.mat'])
    else
        return
    end
    
    [data_bpfreq] = adi_bpfilter_EEG_statistics(cleanEEG_interp, freq);
    
    cfg = [];
    cfg.latency = [-0.5 1];
    data_bpfreq_red = ft_selectdata(cfg, data_bpfreq)
    clear data_bpfreq
    if isfield (data_bpfreq_red, 'grad')
        data_bpfreq_red = rmfield(data_bpfreq_red, 'grad');
    end
    
    % load template vol:
    if ~exist('template_grid', 'var')
        load('\\nas-fa0efsusr1\herfurkn1\My Documents\MATLAB\Bibliotheken\fieldtrip-20180607\template\headmodel\standard_bem');
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
        mri_resliced = ft_volumereslice(cfg, mri_realigned);
        mri_resliced = ft_convert_units(mri_resliced, 'cm');
        elec = ft_convert_units(data_bpfreq_red.elec, 'cm');
        data_bpfreq_red.elec = elec;

        cfg                = [];
        cfg.grid.warpmni   = 'yes';
        cfg.grid.template  = template_grid;
        cfg.grid.nonlinear = 'yes';
        cfg.mri            = mri_resliced; % nicht sicher
        sourcemodel        = ft_prepare_sourcemodel(cfg);
    end
    if ~exist('hdm_ind', 'var')
        hdm_ind = load ([path2vol, 'vol_EEG.mat']);
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
    
    ft_plot_sens(data_bpfreq_red.elec,'style','g*', 'coil', 1);% plot the sensor array
    
    view ([0 -90 0])
    
    outPath = ([SubjOutPath 'run' run '\SVM_result\' ]);
    if ~exist(outPath, 'file')
        mkdir (outPath)
    end
    
    savefig([SubjOutPath 'run' run '\SVM_result\check_vol_sourcemodel_elec.fig']);
    fig = ([SubjOutPath 'run' run '\SVM_result\check_vol_sourcemodel_elec']);
    print('-dpng', fig); 
    close all


    %We first create the leadfield using ft_prepare_leadfield using the individual head model from the previous step, the sensor array and the sourcemodel.
    cfg                 = [];
    cfg.channel         = data_bpfreq_red.label;% ensure that rejected sensors are not present
    cfg.elec            = data_bpfreq_red.elec;
    cfg.vol             = hdm_ind;
    cfg.lcmv.reducerank = 3; % default for MEG is 2, for EEG is 3
    cfg.grid = sourcemodel;
    % cfg.grid.pos=rois; 
    [grid] = ft_prepare_leadfield(cfg);

    %%
% 
%     cfg = [];
%     cfg.toilim = prestim;
%     datapre = ft_redefinetrial(cfg, data_bpfreq_red);
%     cfg.toilim = poststim;
%     datapost = ft_redefinetrial(cfg, data_bpfreq_red);
    
    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = [-.5 1]; % oder von 0 bis 1?
    cfg.keeptrials = 'yes';
    avg_data_bpfreq_red = ft_timelockanalysis(cfg, data_bpfreq_red);
    
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
    cfg.channel = data_bpfreq_red.label;
    cfg.lcmv.fixedori='yes';
    cfg.lcmv.lamda='5%';
    source_avg = ft_sourceanalysis(cfg, avg_data_bpfreq_red);
    
%     cfg = [];
%     cfg.method = 'lcmv';
%     cfg.grid = grid;
%     cfg.grid.filter = source_avg.avg.filter;
%     cfg.headmodel = hdm_ind;
%     cfg.lcmv.projectnoise = 'yes';
%     source_pre = ft_sourceanalysis(cfg, avgpre);
%     source_post = ft_sourceanalysis(cfg, avgpst);
    
    outPath_extdisc_subj = [ outPath_extdisc 'sourcespace\run' run filesep];
    if ~exist('outPath_extdisc_subj', 'dir')
        mkdir (outPath_extdisc_subj)
    end
    
%     outPath_extdisc_subj_time = [ outPath_extdisc 'sourcespace\_run' run filesep prestim '_' poststim ];
%     if ~exist('outPath_extdisc_subj_time', 'dir')
%         mkdir (outPath_extdisc_subj_time)
%     end
%     save([outPath_extdisc_subj_time 'source_pre' ], 'source_pre')
%     save([outPath_extdisc_subj_time 'source_post' ], 'source_post')
    
    spatialfilter_orig = source_avg.avg.filter;
    save([outPath_extdisc_subj 'spatialfilter_fixed_orientation_singletrials_' condition '_' freq], 'spatialfilter_orig')
    
    clear spatialfilter_orig
    
    %% unconstrained:
    cfg = [];
    cfg.method = 'lcmv';
    cfg.grid = grid;  % Stefan: cfg.grid = sourcemodel
    cfg.vol = hdm_ind; % Stefan: cfg.headmodel = vol;
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectnoise = 'yes';
    cfg.channel = data_bpfreq_red.label;
    cfg.lcmv.fixedori='no';
    cfg.lcmv.lamda='5%';
    source_avg_unconstrained = ft_sourceanalysis(cfg, avg_data_bpfreq_red);
    
    spatialfilter_orig = source_avg_unconstrained.avg.filter;
    save([outPath_extdisc_subj 'spatialfilter_loose_orientation_singletrials_' condition '_' freq], 'spatialfilter_orig')
    clear spatialfilter_orig
    
    % On the basis of the computed filters, kept in the output, it is now possible to multiply them with the data. This operation will yield time series commonly known as virtual sensors.
    
    spatialfilter = cat(1,source_avg.avg.filter{:});
%     spatialfilter = spatialfilter(:, 1:248);
%     
%     for k = 1:length(data_bpfreq_red.trial)
%         data_bpfreq_red.trial{1,k} = data_bpfreq_red.trial{1,k}(1:248,:)
%     end
    
    virtsens = [];
    for i = 1:length(data_bpfreq_red.trial)
        virtsens.trial{i} = spatialfilter*data_bpfreq_red.trial{i};
    end
    virtsens.time=data_bpfreq_red.time;
    virtsens.fsample=data_bpfreq_red.fsample;

    for i=1:length(virtsens.trial{1}(:,1))
        virtsens.label{i}=num2str(i);
    end

    % siehe Skript von Yuval: => all will have a similar noise level
    ns = mean(abs(spatialfilter),2);

    for k = 1:length(virtsens.trial)
        virtsens_ns.trial{k} = virtsens.trial{1,k}./repmat(ns,1,size(virtsens.trial{1,k},2)); % => all will have a similar noise level
    end

    virtsens_ns.time = virtsens.time;
    virtsens_ns.fsample = virtsens.fsample;
    
    for i = 1:length(virtsens_ns.trial{1}(:,1))
        virtsens_ns.label{i} = num2str(i);
    end
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



function [Condition1vs2] = adi_crossvalidation (condition1run, condition2run, freqbandname, NameCond1vs2, latency, run, SubjOutPath)

    cfg             = [];
    cfg.parameter   = 'trial';
    cfg.keeptrials  = 'yes'; % classifiers operate on individual trials
    cfg.channel     = 'EEG';
    cfg.vartrllength = 2;
    tCondition1     = ft_timelockanalysis(cfg, condition1run);    
    tCondition2     = ft_timelockanalysis(cfg, condition2run); 

    cfg         = [];
    cfg.method  = 'crossvalidate'; % layout braucht es nicht, verschiedene layouts führen zum gleichen ergebnis
    cfg.channel = 'EEG';
    cfg.statistic = {'accuracy', 'binomial', 'contingency'};
    cfg.design  = [ones(size(condition1run.trial,2),1); 2*ones(size(condition2run.trial,2),1)]';
    cfg.resample = 'true';
    %     cfg.mva = {dml.standardizer dml.enet 'family', 'binomial', 'alpha', 0.3};

    Condition1vs2 = [];
    Condition1vs2.Accuracy=[];
    Condition1vs2.Binominal=[];
    Condition1vs2.Latency=latency;
    Condition1vs2.stats='5f-crossvalidation';

    for i = 1:length(latency)  
        cfg.latency = [latency(1,i) latency(2,i)]; 
        stat = ft_timelockstatistics(cfg, tCondition1, tCondition2);
        Condition1vs2.Accuracy(1,i) = stat.statistic.accuracy;
        Condition1vs2.Binominal(1,i) = stat.statistic.binomial;
        Condition1vs2.Contingency{1,i} = stat.statistic.contingency;
    end

    Condition1vs2.latency = latency;
    Condition1vs2.design = cfg.design;

    outPath = ([SubjOutPath 'run' run '\SVM_result\virtSens\' ]);
    if ~exist(outPath, 'dir')
        mkdir (outPath)
    end
    save ([outPath filesep NameCond1vs2 '_' run '_' freqbandname '.mat'], 'Condition1vs2'); 


end


function [] = adi_figureTPRcrossval_SVM (Condition1vs2, time, cond1, cond2, NameCond1vs2, outPath, freqbandname, run)
 
figure
plot(time, Condition1vs2.Accuracy);
hold on

indSig = []; 
indSig(1:length(Condition1vs2.Binominal)) = NaN;
sig = find(Condition1vs2.Binominal <= 0.1);
indSig(sig) = Condition1vs2.Binominal(sig);
plot(time, indSig,'r+'); 

title([NameCond1vs2, '_', freqbandname]);
xlabel('time');
ylabel('accuracy/p-value'); 
ylim([0 1]);

savefig([outPath, '\run' run '\SVM_result\virtSens\' NameCond1vs2, run, '_', freqbandname,'.fig'])
fig = ([outPath, '\run' run '\SVM_result\virtSens\' NameCond1vs2, run, '_', freqbandname]);
print('-dpng', fig); 

total_cond1 = Condition1vs2.Contingency{1,1}(1,1) + Condition1vs2.Contingency{1,1}(1,2);
total_cond2 = Condition1vs2.Contingency{1,1}(2,1) + Condition1vs2.Contingency{1,1}(2,2);

for i = 1:length(Condition1vs2.Contingency)
    tpr_cond1(i)= Condition1vs2.Contingency{1,i}(1,1)/total_cond1; % true positive rate
    tpr_cond2(i)= Condition1vs2.Contingency{1,i}(2,2)/total_cond2; % true positive rate
end

figure
plot(time, Condition1vs2.Accuracy, 'r');
hold on
plot(time, indSig,'r+');  
hold on
plot(time, tpr_cond1,'b'); 
hold on
plot(time, tpr_cond2,'k'); 
leg_cond1 = (['TPR_', cond1]);
leg_cond2 = (['TPR_', cond2]);
legend ({'total TPR', 'p-value of total TPR', leg_cond1, leg_cond2});
title([NameCond1vs2, '_', freqbandname]);
xlabel('time');
ylabel('accuracy/p-value'); 
ylim ([0 1]);
savefig([outPath '\run' run '\SVM_result\virtSens\' NameCond1vs2 '_' run '_' freqbandname '_TPR.fig']);
fig = ([outPath '\run' run '\SVM_result\virtSens\' NameCond1vs2 '_' run '_' freqbandname '_TPR']);
print('-dpng', fig); 
close all
    
end


function [data_bpfreq] = adi_bpfilter_EEG_statistics(data, bpname)


switch bpname
    case 'bp1-45Hz'
        data_bpfreq = data;
        return
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
end

cfg=[];
cfg.trials        = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') %|| 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

[data_bpfreq] = ft_preprocessing(cfg, data);

end
