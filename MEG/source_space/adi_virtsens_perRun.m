
function adi_source_SVM (path2vol, path2data, mriPath, outPath_extdisc, freqbandname, condition)

% calculations spatial filter per Run and Subject, 
% time =  latency(1,:);

%% run 1:
if ~exist([outPath_extdisc 'virtsens\run1\virtsens_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath,  num2str(1), freqbandname, condition, outPath_extdisc)    
end

%% run2:

if ~exist([outPath_extdisc 'virtsens\run2\virtsens_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(2), freqbandname, condition, outPath_extdisc)
end

%% run3:

if ~exist([outPath_extdisc 'virtsens\run3\virtsens_' freqbandname '.mat'], 'file')
    adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, num2str(3), freqbandname, condition, outPath_extdisc)
end

end

function [virtsens] = adi_source_reconstruction_virt_sens(path2vol, path2data, mriPath, run, freq, condition, outPath_extdisc)
    
% load data:
file_name = ['Neu_Like' '500_' run];
dir_run = dir([path2data '*.mat']);
for i=1:length(dir_run)
    if 1==contains(dir_run(i).name, file_name)
        load ([path2data dir_run(i).name], 'cleanMEG_interp'); 
        [data_like_bpfreq] = adi_bpfilter(cleanMEG_interp, freq);
        clear cleanMEG_interp 
    end

end

file_name = ['Neu_Dislike' '500_' run];
for i=1:length(dir_run)
    if 1==contains(dir_run(i).name, file_name)
        load ([path2data dir_run(i).name], 'cleanMEG_interp'); 
        [data_dislike_bpfreq] = adi_bpfilter(cleanMEG_interp, freq);
        clear cleanMEG_interp 
    end

end

if ~exist('data_dislike_bpfreq', 'var')
    return
end

data_all_conditions.trial = [data_like_bpfreq.trial  data_dislike_bpfreq.trial];
data_all_conditions.time = [data_like_bpfreq.time  data_dislike_bpfreq.time];
data_all_conditions.response = [data_like_bpfreq.trialinfo.response data_dislike_bpfreq.trialinfo.response];
data_all_conditions.response_label = [data_like_bpfreq.trialinfo.response_label data_dislike_bpfreq.trialinfo.response_label];
data_all_conditions.triggerlabel = [data_like_bpfreq.trialinfo.triggerlabel data_dislike_bpfreq.trialinfo.triggerlabel];
data_all_conditions.balldesign = [data_like_bpfreq.trialinfo.balldesign data_dislike_bpfreq.trialinfo.balldesign];
data_all_conditions.fsample = data_like_bpfreq.fsample;
data_all_conditions.grad = data_like_bpfreq.grad;
data_all_conditions.label = data_like_bpfreq.label;
  
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
        grad = ft_convert_units(data_all_conditions.grad, 'cm');
        data_all_conditions.grad = grad; % geändert

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
    
    ft_plot_sens(data_all_conditions.grad,'style','g*', 'coil', 1);% plot the sensor array
    
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
    cfg.channel         = data_all_conditions.label;% ensure that rejected sensors are not present
    cfg.grad            = data_all_conditions.grad;
    cfg.vol             = hdm_ind;
    cfg.lcmv.reducerank = 2; % default for MEG is 2, for EEG is 3
    cfg.grid = sourcemodel;
    % cfg.grid.pos=rois; 
    [leadfield] = ft_prepare_leadfield(cfg);% =grid

    %%
  
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
    cfg.vartrllength = 2;

    try
        avg_data_all_conditions = ft_timelockanalysis(cfg, data_all_conditions);
    catch ME
        if 1==strcmp(ME.message,'data has variable trial lengths, you specified not to accept that')
           for j = 1:length(data_all_conditions.time)
            [ind(j)] =  isequal(length(data_all_conditions.time{j}), 3052);
           end
           trl_ind = find(ind == false);
           for j=1:length(trl_ind)
               data_all_conditions.time{1, trl_ind(j)}(:,1)=[];
               data_all_conditions.trial{1,trl_ind(j)}(:,1)=[];
           end
           avg_data_all_conditions = ft_timelockanalysis(cfg, data_all_conditions);     
        end
    end
    cfg = [];
    avg = ft_timelockanalysis(cfg, data_all_conditions);
    [FourRef,Fref] = fftBasic(avg.avg, round(data_all_conditions.fsample));
    figure;
    plot(FourRef, Fref)
    PathFig = [outPath_extdisc 'source_avg\run' run '\'];
    if ~exist(PathFig, 'dir')
        mkdir(PathFig)
    end
    savefig([PathFig '\freqspectrum_avg_filtered_trials.fig'])
    close
    
    % Now we make a first call to ft_sourceanalysis in order to compute the spatial filters on the basis of the entire data and keep them in the output for a later use.
    cfg = [];
    cfg.method = 'lcmv';
%     cfg.latency = [-0.5 1];
    cfg.grid = leadfield;  % Stefan: cfg.grid = sourcemodel
    cfg.vol = hdm_ind; % Stefan: cfg.headmodel = vol;
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.projectnoise = 'yes';
    cfg.channel = data_all_conditions.label;
    cfg.keeptrials = 'yes';
    cfg.lcmv.fixedori='yes';
    cfg.lcmv.lambda='5%';
    try
        source_avg = ft_sourceanalysis(cfg, avg_data_all_conditions);
    catch

    end

    outPath_extdisc_subj = [ outPath_extdisc filesep 'source_avg' filesep 'run' run filesep];
   
    save([outPath_extdisc_subj 'source_avg_appended_conditions_' freq], 'source_avg')
    spatialfilter = cat(1,source_avg.avg.filter{:});   
    
        % laut fieldtrip discussion list zuerst weights mit sensor data multiplizieren: https://mailman.science.ru.nl/pipermail/fieldtrip/2017-July/011661.html    
    virtsens = [];
    for k = 1:length(data_all_conditions.trial)
        virtsens.trial{k} = spatialfilter*data_all_conditions.trial{k};
    end
    
    label = [1:size(virtsens.trial{1,1},1)]';
    for k = 1:length(label)
        virtsens.label{k,1} = num2str(label(k));
    end

    
    fn_data_all_conditions = fieldnames(data_all_conditions);
    fn_virtsens = fieldnames(virtsens);

    diff_fieldnames = setdiff(fn_data_all_conditions, fieldnames(virtsens));

    for k=1:length(diff_fieldnames)
        virtsens.(diff_fieldnames{k}) = data_all_conditions.(diff_fieldnames{k});
    end

       
    [pxx,f] = pwelch(virtsens.trial{1,1}',[],[],[],data_all_conditions.fsample); 
    figure
    semilogy(f, pxx)
    Pathvirtsens =  [ outPath_extdisc filesep 'virtsens' filesep 'run' run filesep ];
    if ~exist(Pathvirtsens, 'dir')
        mkdir(Pathvirtsens)
    end
    
    savefig([Pathvirtsens 'freqspectrum_virtsens_' freq '.fig'])
    close
    
    save([Pathvirtsens 'virtsens_' freq '.mat'], 'virtsens')
     % sanity check:
    
    cfg = [];
    cfg.vartrllength = 2;
    avg = ft_timelockanalysis(cfg, virtsens);
    figure;
    plot(virtsens.time{1,1}, avg.avg)
    sanity_path = [Pathvirtsens '\sanity_check\'];
    if ~exist(sanity_path, 'dir')
        mkdir (sanity_path)
    end
    savefig([sanity_path 'virtsens_all_conditions_run' num2str(run) '_' freq '.fig'])
    close
 
%% noise normalization nach Gespräch mit Stefan:
%     abs_timevector =  abs(data_all_conditions.time{1,1});
%     value_closest_zero = min(abs_timevector);
%     [~, ind_closest_zero] = find(abs_timevector == value_closest_zero);
%     
%     virtsens_ns = virtsens;
%     virtsens_ns = rmfield(virtsens_ns, 'trial');
%     for k = 1:length(virtsens.trial)
%         for p = 1:size(virtsens.trial{k},1)
%             virtsens_ns.trial{k}(p,:) = virtsens.trial{k}(p,:)./mean(virtsens.trial{k}(p,1:ind_closest_zero));
%         end
%     end
%     clearvars virtsens
    
%      % sanity check nr. 2:
%     cfg = [];
%     cfg.vartrllength = 2;
%     avg_ns = ft_timelockanalysis(cfg, virtsens_ns);
%     figure;
%     plot(virtsens_ns.time{1,1}, avg_ns.avg)
%     savefig([sanity_path 'virtsens_ns_all_conditions_run' num2str(run) '_' freq '.fig'])
%     close
%     
%     
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


for k = 1:length(data_bpfreq_res_sel.trial)
    data_bpfreq_res_sel.grad.label(249:end) = [];
    data_bpfreq_res_sel.grad.chanori(249:end, :) = [];
    data_bpfreq_res_sel.grad.chanpos(249:end, :) = [];
    data_bpfreq_res_sel.grad.tra(249:end, :) = [];
    data_bpfreq_res_sel.label(249:end) = [];
end





end
     


