
function adi_source_SVM (path2vol, path2data, mriPath, outPath_extdisc, freqbandname, condition)

% calculations spatial filter per Run and Subject, 
% time =  latency(1,:);

%% run 1:
if ~exist([outPath_extdisc '\run1\spatialfilter_fixed_orientation_singletrials_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath,  num2str(1), freqbandname, condition, outPath_extdisc)    
end

%% run2:

if ~exist([outPath_extdisc '\run2\spatialfilter_fixed_orientation_singletrials_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(2), freqbandname, condition, outPath_extdisc)
end

%% run3:

if ~exist([outPath_extdisc '\run3\spatialfilter_fixed_orientation_singletrials_' freqbandname '.mat'], 'file')
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
   load ([path2data, 'Neu_Like' '500_' run '.mat'], 'cleanMEG_interp');
   if strcmp(freq, 'bp1_95Hz') 
       data.like = cleanMEG_interp;
   else
       [data_like_bpfreq] = adi_bpfilter(cleanMEG_interp, freq);
        data.like = data_like_bpfreq;
   end
   clear cleanMEG_interp 
else
    condition{1} = [];
end
if exist([path2data, 'Neu_Dislike' '500_' run '.mat'], 'file')
    index(2) = 1;
    condition{2} = 'data_dislike_bpfreq';
    load ([path2data, 'Neu_Dislike' '500_' run '.mat'], 'cleanMEG_interp');
    if strcmp(freq, 'bp1_95Hz') 
        data.dislike = cleanMEG_interp;
    else
    [data_dislike_bpfreq] = adi_bpfilter(cleanMEG_interp, freq);
         data.dislike = data_dislike_bpfreq;
    end
    clearvars cleanMEG_interp 
else
     condition{2} = [];
end
% if exist([path2data, 'dontcare' '500_' run '.mat'], 'file')
%     index(3) = 1;
% %     condition{3} = 'data_dontcare_bpfreq';
%     data_dontcare = load ([path2data, 'dontcare' '500_' run '.mat'], 'cleanMEG_interp');
%     [data_dontcare_bpfreq] = adi_bpfilter(data_dontcare.cleanMEG_interp, freq);
%     data.dontcare = data_dontcare_bpfreq;
%     clear data_dontcare data_dontcare_bpfreq
% else
%        condition{3} = [];
% end

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
    data_appended = ft_appenddata(cfg, data.(fieldsnames{1}), data.(fieldsnames{2}));
end   

fsample=data.(fieldsnames{1}).fsample;

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
        data_appended.grad = grad; % ge�ndert

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
  
    cfg = [];
    cfg.covariance = 'yes'; % compute the covariance for single trials, then average
    cfg.covariancewindow = [0 1]; % oder von 0 bis 1?  da ich mich f�r like vs dislike interessiere, evtl. Zeitbereich 0-1s, bei pr�-vs.poststim-intervall-Konstrast pr�stim einbeziehen
    cfg.removemean = 'no'; % andernfalls zieht er von den Daten den Mittelwert ab zur Berechnung der Kovarianz
    cfg.vartrllength = 2;
    cfg.keeptrials = 'yes';
%     cfg.preproc.bpfilter = 'yes';
%     cfg.preproc.bpfreq = [1 95];
%     cfg.preproc.demean = 'no';
%     cfg.preproc.baselinewindow = [-0.5 0]; % baseline correction
    cfg.vartrllength = 2;

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
    cfg = [];
    avg = ft_timelockanalysis(cfg, data_appended);
    [FourRef,Fref] = fftBasic(avg.avg, round(fsample));
    figure;
    plot(FourRef, Fref)
    PathFig = [outPath_extdisc '\run' run '\'];
    if ~exist(PathFig, 'dir')
        mkdir(PathFig)
    end
    savefig([outPath_extdisc '\run' run '\freqspectrum_avg_filtered_trials.fig'])
    close
    
    % Now we make a first call to ft_sourceanalysis in order to compute the spatial filters on the basis of the entire data and keep them in the output for a later use.
    cfg = [];
    cfg.method = 'lcmv';
%     cfg.latency = [-0.5 1];
    cfg.grid = leadfield;  % Stefan: cfg.grid = sourcemodel
    cfg.vol = hdm_ind; % Stefan: cfg.headmodel = vol;
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectnoise = 'yes';
    cfg.channel = data_appended.label;
    cfg.keeptrials = 'yes';
    cfg.lcmv.fixedori='yes';
    cfg.lcmv.lambda='5%';
    try
        source_avg = ft_sourceanalysis(cfg, avg_data_appended);
    catch

    end

    outPath_extdisc_subj = [ outPath_extdisc 'run' run filesep];
    if ~exist(outPath_extdisc_subj, 'dir')
        mkdir (outPath_extdisc_subj)
    end
    
    
    spatialfilter_orig = source_avg.avg.filter;
%     save([outPath_extdisc_subj 'spatialfilter_loose_orientation_singletrials_' freq], 'spatialfilter_orig')
    save([outPath_extdisc_subj 'source_avg_appended_conditions_' freq], 'source_avg')
    
    
end


function [data_bpfreq_res_sel] = adi_bpfilter(filename, bpname)


switch bpname
    case 'bp1_95Hz'
       data_bpfreq_res_sel =  filename;
        return
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
cfg.keeptrials = 'yes';
cfg.vartrllength = 2;

cfg.trials  = 'all'; 
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





end
     


