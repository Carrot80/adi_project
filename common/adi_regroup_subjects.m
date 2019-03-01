
function  [vs_allRuns] = adi_appendvirtsensMEG(vs_allRuns, path2sensordata, path2spatialfilter, dir_out, subject, pattern, trigger, freq, run)

if ~exist ([path2spatialfilter 'run' num2str(run)], 'dir')
    return
end

list_run = dir(fullfile([path2sensordata '*' run '*.mat']));

for k=1:length(list_run)
    data(k) = load ([path2sensordata list_run(k).name]);
end

switch length(data)
    case 3
        data_all_conditions.trial = [data(1).cleanMEG_interp.trial  data(2).cleanMEG_interp.trial data(3).cleanMEG_interp.trial];
        data_all_conditions.time = [data(1).cleanMEG_interp.time  data(2).cleanMEG_interp.time data(3).cleanMEG_interp.time];
        data_all_conditions.response = [data(1).cleanMEG_interp.trialinfo.response data(2).cleanMEG_interp.trialinfo.response data(3).cleanMEG_interp.trialinfo.response];
        data_all_conditions.response_label = [data(1).cleanMEG_interp.trialinfo.response_label data(2).cleanMEG_interp.trialinfo.response_label data(3).cleanMEG_interp.trialinfo.response_label];
        data_all_conditions.triggerlabel = [data(1).cleanMEG_interp.trialinfo.triggerlabel data(2).cleanMEG_interp.trialinfo.triggerlabel data(3).cleanMEG_interp.trialinfo.triggerlabel];
        data_all_conditions.balldesign = [data(1).cleanMEG_interp.trialinfo.balldesign data(2).cleanMEG_interp.trialinfo.balldesign data(3).cleanMEG_interp.trialinfo.balldesign];
    case 2
        data_all_conditions.trial = [data(1).cleanMEG_interp.trial  data(2).cleanMEG_interp.trial];
        data_all_conditions.time = [data(1).cleanMEG_interp.time  data(2).cleanMEG_interp.time];
        data_all_conditions.response = [data(1).cleanMEG_interp.trialinfo.response data(2).cleanMEG_interp.trialinfo.response];
        data_all_conditions.response_label = [data(1).cleanMEG_interp.trialinfo.response_label data(2).cleanMEG_interp.trialinfo.response_label];
        data_all_conditions.triggerlabel = [data(1).cleanMEG_interp.trialinfo.triggerlabel data(2).cleanMEG_interp.trialinfo.triggerlabel];
        data_all_conditions.balldesign = [data(1).cleanMEG_interp.trialinfo.balldesign data(2).cleanMEG_interp.trialinfo.balldesign];
    case 1
        data_all_conditions.trial = data(1).cleanMEG_interp.trial;
        data_all_conditions.time = data(1).cleanMEG_interp.time;
        data_all_conditions.response = [data(1).cleanMEG_interp.trialinfo.response];
        data_all_conditions.response_label = [data(1).cleanMEG_interp.trialinfo.response_label];
        data_all_conditions.triggerlabel = [data(1).cleanMEG_interp.trialinfo.triggerlabel];
        data_all_conditions.balldesign = [data(1).cleanMEG_interp.trialinfo.balldesign];
end
    data_all_conditions.label = data(1).cleanMEG_interp.label;
    data_all_conditions.fsample = data(1).cleanMEG_interp.fsample;
    data_all_conditions.dimord = data(1).cleanMEG_interp.dimord;
    
    clearvars data
    
    
    [vs_allRuns] = adi_appendvirtsensMEG(vs_allRuns, data_all_conditions, path2sensordata, path2spatialfilter, dir_out, subject, pattern, trigger, freq, run);

    
    
    
end

function [vs_allRuns] = adi_appendvirtsensMEG(vs_allRuns, data_all_conditions, path2sensordata, path2spatialfilter, dir_out, subject, pattern, trigger, freq, run)


% if ~exist([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'file')

    % spatial filter wird geladen
    % data_all_runs werden aufgeteilt in Muster volley, space und
    % soccerball(=football)    

%     [data_all_conditions] = regroup_trials(data_all_conditions, pattern, trigger);
%     [data_all_conditions_bpfreq] = adi_bpfilter(data_all_conditions, freq, dir_out, run);
    data_all_conditions_bpfreq = data_all_conditions;
    clearvars data_all_conditions     
    
    load ([path2spatialfilter 'run' num2str(run) '\source_avg_appended_conditions_' freq '.mat'], 'source_avg');
    spatialfilter = cat(1,source_avg.avg.filter{:});   

    
%%
    % laut fieldtrip discussion list zuerst weights mit sensor data multiplizieren: https://mailman.science.ru.nl/pipermail/fieldtrip/2017-July/011661.html    
    virtsens = [];
    for k = 1:length(data_all_conditions_bpfreq.trial)
        virtsens.trial{k} = spatialfilter*data_all_conditions_bpfreq.trial{k};
    end

    for p = 1:length(virtsens.trial)
        euclid_norm = zeros(length(virtsens.trial{1,p})/3, size(virtsens.trial{1,p}(1,:),2));
        n = 1;
        for k = 1:3:length(virtsens.trial{1,p})
            euclid_norm(n,:) = sqrt(virtsens.trial{1,p}(k,:).^2+virtsens.trial{1,p}(k+1,:).^2+virtsens.trial{1,p}(k+2,:).^2);
            n = n+1;
        end  
        vs_euclid_norm.trial{1,p} = euclid_norm;
        clear euclid_norm
    end
    virtsens.trial = vs_euclid_norm.trial;
    clearvars vs_euclid_norm

    virtsens.label = cell(length(virtsens.trial{1,1}(:,1)), 1);
    for n = 1:length(virtsens.trial{1,1}(:,1))
        virtsens.label{n,1}=num2str(n);
    end
  
    
    
    
    %%
    fn_data_all_conditions_bpfreq = fields(data_all_conditions_bpfreq);
    fn_virtsens = fields(virtsens);
    diff_fieldnames = setdiff(fn_data_all_conditions_bpfreq, fn_virtsens);
    
    for k=1:length(diff_fieldnames)
        virtsens.(diff_fieldnames{k}) = data_all_conditions_bpfreq.(diff_fieldnames{k}); 
    end
  
    % sanity check:
    
    cfg = [];
    cfg.vartrllength = 2;
    avg = ft_timelockanalysis(cfg, virtsens);
    figure;
    plot(virtsens.time{1,1}, avg.avg)
    sanity_path = [dir_out '\sanity_check\'];
    if ~exist(sanity_path, 'dir')
        mkdir (sanity_path)
    end
    savefig([sanity_path 'virtsens_all_conditions_run' num2str(run) '_' freq '.fig'])
    close
 
%% noise normalization nach Gespräch mit Stefan:
    virtsens_ns = virtsens;
    virtsens_ns = rmfield(virtsens_ns, 'trial');
    for k = 1:length(virtsens.trial)
        for p = 1:size(virtsens.trial{k},1)
            virtsens_ns.trial{k}(p,:) = virtsens.trial{k}(p,:)/mean(virtsens.trial{k}(p,1:129));
        end
    end
    clearvars virtsens
    
     % sanity check nr. 2:
    cfg = [];
    cfg.vartrllength = 2;
    avg_ns = ft_timelockanalysis(cfg, virtsens_ns);
    figure;
    plot(virtsens_ns.time{1,1}, avg_ns.avg)
    savefig([sanity_path 'virtsens_ns_all_conditions_run' num2str(run) '_' freq '.fig'])
    close
    
    vs_allRuns.(['run' num2str(run)]) = virtsens_ns;
    clearvars virtsens_ns
    
    fsample = vs_allRuns.(['run' num2str(run)]).fsample;
    [FourRef,Fref]=fftBasic(avg_ns.avg, round(fsample));
    figure;
    plot(FourRef, Fref)
    savefig([dir_out 'avg_virtsens_freq_spectrum_run' run '.fig'])
    close
    

      %% svd function aus connectivity tutorial:
   clearvars source_avg 
      
     virtsens = [];
    for k = 1:length(data_all_conditions_bpfreq.trial)
        virtsens.trial{k} = spatialfilter*data_all_conditions_bpfreq.trial{k};
    end  
    
    virtualchanneldata = [];
    virtualchanneldata.label = {'cortex'};
    virtualchanneldata.time = data_all_conditions_bpfreq.time;
    
for k=1:length(virtsens.trial)    
    n=1;
    for p=1:3:length(virtsens.trial{k})
        vs_timeseries = cat(2, virtsens.trial{k}(p:p+2,:));
        [u, s, v] = svd(vs_timeseries, 'econ'); 
        timeseriesmaxproj = u(:,1)' * vs_timeseries;
        virtualchanneldata.trial{1,k}(n,:) = timeseriesmaxproj;
        clear timeseriesmaxproj u s v vs_timeseries
        n=n+1;
    end
end    

  % sanity check:
    
    figure;
    plot(virtualchanneldata.time{1,1}, abs(virtualchanneldata.trial{1,1}))

    fsample = vs_allRuns.(['run' num2str(run)]).fsample;
    [FourRef,Fref]=fftBasic(virtualchanneldata.trial{1,1}, round(fsample));
    figure;
    plot(FourRef, Fref)
    
    
    virtualchanneldata_noisecorr = virtualchanneldata;
    virtualchanneldata_noisecorr = rmfield(virtualchanneldata_noisecorr, 'trial');
    for k = 1:length(virtualchanneldata.trial)
        for p = 1:size(virtualchanneldata.trial{k},1)
            virtualchanneldata_noisecorr.trial{k}(p,:) = virtualchanneldata.trial{k}(p,:)/mean(virtualchanneldata.trial{k}(p,1:129));
        end
    end
    
    % sanity check no 2:
    
    figure;
    plot(virtualchanneldata_noisecorr.time{1,1}, abs(virtualchanneldata_noisecorr.trial{1,1}))
    

    

end


function [data_bpfreq_sel_res] = adi_bpfilter(filename, bpname, dir_out, run)

for k = 1:length(filename.trial)
    filename.trial{k}(249:end,:)=[];
end  

% for k = 1:length(filename.trial)
% %     data_bpfreq.grad.label(249:end) = [];
% %     data_bpfreq.grad.chanori(249:end, :) = [];
% %     data_bpfreq.grad.chanpos(249:end, :) = [];
% %     data_bpfreq.grad.tra(249:end, :) = [];
%     filename.label(249:end) = [];
% end

% for m = 1:length(filename.trial)
%     if 3053 == length(filename.time{1,m}) 
%         filename.trial{1,m}(:,1) = [];
%         filename.time{1,m}(:,1) = [];
%     end
% end

switch bpname
    case 'bp1_45Hz'
        bpfreq = [1 45];
    case '1_5_45Hz'
        bpfreq = [1.5 45];
    case 'bp2_45Hz'
        bpfreq = [2 45];
    case 'delta'
        bpfreq = [1 4];
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


% trl_nan = zeros(length(filename.trial));
% for m = 1:length(filename.trial)
%  if isnan(filename.trial{m}(1,1))
%     trl_nan(m) = m;
%     [ind] = sum(isnan(filename.trial{m}(1,:)));
%     filename.trial{m}(:,1:ind) = [];
%     filename.time{m}(1:ind) = [];
%  end
% end

cfg=[];
cfg.trials        = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta')% || 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

[data_bpfreq] = ft_preprocessing(cfg, filename);

cfg = [];
avg = ft_timelockanalysis(cfg, data_bpfreq);
figure;
plot(avg.time, avg.avg)
path_avg = [dir_out '\sanity_check\avg\'];
if ~exist ( path_avg, 'dir')
    mkdir (path_avg)
end
savefig([path_avg 'avg_' bpname '_run' run '.fig'])
 

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_sel = ft_selectdata(cfg, data_bpfreq);

den Bereich nochmal anschauen:
sRate = 1017.25;
[FourRef,Fref]=fftBasic(data_bpfreq_sel.trial{1,1},round(sRate));
figure
plot(FourRef, Fref)


cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'yes';
[data_bpfreq_sel_res] = ft_resampledata(cfg, data_bpfreq_sel);

fn_data_bpfreq_sel_res = fieldnames(data_bpfreq_sel_res);
fn_filename = fieldnames(filename);

diff_fieldnames = setdiff(fn_filename, fn_data_bpfreq_sel_res);

for k=1:length(diff_fieldnames)
    data_bpfreq_sel_res.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
end

       
end
     
