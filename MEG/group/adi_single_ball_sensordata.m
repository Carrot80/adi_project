
function  [sensordata_ball_freq_allsubj] = append_sensordata(sensordata_ball_freq_allsubj, path2sensor_data, subject, pattern, trigger, freq, balldesign, delete_runs, subjnum)

dir_runs = dir([path2sensor_data '*mat']);
if isempty(delete_runs)
    for k = 1:length(dir_runs)
        sensor_data(k)= load ([path2sensor_data  dir_runs(k).name], 'cleanMEG_interp');
        for p = 1:length(sensor_data(k).cleanMEG_interp.trial)
            sensor_data(k).cleanMEG_interp.run{p} = dir_runs(k).name(end-4);
            sensor_data(k).cleanMEG_interp.subject{p} = subject;
        end
    end
% for i = 1:length(delete_runs)
%     if 1==strcmp(delete_runs{i,1}, subject)
%         for k = 1:size(dir_runs,1)
%             temp=contains(dir_runs(k).name, ['500_' delete_runs{i,2}]);
%             ind_(k)=temp;
%             clearvars temp
%         end
%             ind =find(ind_);
%             dir_runs(ind)=[];
%             clearvars ind ind_
%     end
% end
else

    z=1;
    for k = 1:length(dir_runs)
        if 0 == delete_runs.(subject).(['run' dir_runs(k).name(end-4)]).(balldesign)
            sensor_data(k)= load ([path2sensor_data  dir_runs(k).name], 'cleanMEG_interp');

            for p = 1:length(sensor_data(k).cleanMEG_interp.trial)
                sensor_data(k).cleanMEG_interp.run{p} = dir_runs(k).name(end-4);
                sensor_data(k).cleanMEG_interp.subject{p} = subject;
            end
        else
            files_2_delete(z) = k;
            z = z+1;
        end
    end

    if ~exist('sensor_data', 'var')
        return
    end

    if exist('files_2_delete', 'var')
        while files_2_delete(end) > length(sensor_data)
              files_2_delete(end) =[];
        end
        sensor_data(files_2_delete) = [];
    end

end

sensordata_allruns = [];
sensordata_allruns.trial = sensor_data(1).cleanMEG_interp.trial;
sensordata_allruns.time = sensor_data(1).cleanMEG_interp.time;
sensordata_allruns.response_label = sensor_data(1).cleanMEG_interp.trialinfo.response_label;
sensordata_allruns.balldesign_short = sensor_data(1).cleanMEG_interp.trialinfo.balldesign_short;
sensordata_allruns.run = sensor_data(1).cleanMEG_interp.run;
sensordata_allruns.subject = sensor_data(1).cleanMEG_interp.subject;

for k=2:length(sensor_data)
    sensordata_allruns.trial = cat(2,sensordata_allruns.trial, sensor_data(k).cleanMEG_interp.trial);
    sensordata_allruns.time = cat(2,sensordata_allruns.time, sensor_data(k).cleanMEG_interp.time);
    sensordata_allruns.response_label = cat(2,sensordata_allruns.response_label, sensor_data(k).cleanMEG_interp.trialinfo.response_label);
    sensordata_allruns.balldesign_short = cat(2,sensordata_allruns.balldesign_short, sensor_data(k).cleanMEG_interp.trialinfo.balldesign_short); 
    sensordata_allruns.run = cat(2, sensordata_allruns.run, sensor_data(k).cleanMEG_interp.run);
    sensordata_allruns.subject = cat(2, sensordata_allruns.subject, sensor_data(k).cleanMEG_interp.subject);
end

sensordata_allruns.label = sensor_data(1).cleanMEG_interp.label;
sensordata_allruns.fsample = sensor_data(1).cleanMEG_interp.fsample;
sensordata_allruns.grad = sensor_data(1).cleanMEG_interp.grad;

% clearvars sensor_data

for k = 1:length(sensordata_allruns.balldesign_short)
    balldesign_trials(k) = sensordata_allruns.balldesign_short{k};
end

ind_balldesign = find(strcmp(balldesign_trials, balldesign));
if ~isempty(ind_balldesign)
    sensordata_allruns.trial =  sensordata_allruns.trial(ind_balldesign);
    sensordata_allruns.time =  sensordata_allruns.time(ind_balldesign);
    sensordata_allruns.response_label =  sensordata_allruns.response_label(ind_balldesign);
    sensordata_allruns.balldesign_short =  sensordata_allruns.balldesign_short(ind_balldesign);
    sensordata_allruns.run =  sensordata_allruns.run(ind_balldesign);
    sensordata_allruns.subject =  sensordata_allruns.subject(ind_balldesign);
end

% unique(balldesign_trials)

%% delete balldesigns which were not clearly rated:

if ~isempty(delete_runs)

    delete_trials = zeros(1, length(sensordata_allruns.trial));
    for p=1:length(sensordata_allruns.trial)
    %     balldesign_trial=sensordata_allruns.balldesign_short{1,p}{1,1};
        % check if trial should be deleted:
        if 1==delete_runs.(subject).(['run' sensordata_allruns.run{p}]).(sensordata_allruns.balldesign_short{p}{1,1})
           delete_trials(p) = 1;
        end
    end

    trials2keep=find(~delete_trials);
    sensordata_allruns.trial =  sensordata_allruns.trial(trials2keep);
    sensordata_allruns.time =  sensordata_allruns.time(trials2keep);
    sensordata_allruns.response_label =  sensordata_allruns.response_label(trials2keep);
    sensordata_allruns.balldesign_short =  sensordata_allruns.balldesign_short(trials2keep);
    sensordata_allruns.run =  sensordata_allruns.run(trials2keep);
    sensordata_allruns.subject =  sensordata_allruns.subject(trials2keep);
    
    
    % for k=1:length(sensordata_allruns.trial)
    %    temp(k)=strcmp(sensordata_allruns.balldesign_short{k}, balldesign);
    % end
    % ind=find(temp);
    
end

if ~isempty(ind_balldesign)
    [sensordata_allruns_ball_bpfreq] = adi_bpfilter(sensordata_allruns, freq);
else
    clearvars ind temp sensordata_allruns
    return
end

for kk = 1:length(sensordata_allruns_ball_bpfreq.trial)
    sensordata_allruns_ball_bpfreq.trl_no(kk) = kk; 
end

sensordata_ball_freq_allsubj(subjnum)=sensordata_allruns_ball_bpfreq;
clearvars ind temp sensordata_allruns_ball_bpfreq sensordata_allruns
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

try
    [data_bpfreq] = ft_preprocessing(cfg, filename); 
    [warnMsg, warnID] = lastwarn;
    if ~isempty(warnMsg)
       warnMsg
    end
catch
    
    warnMsg
end

cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);

for k = 1:length(data_bpfreq_res_sel.trial)
    data_bpfreq_res_sel.grad.label(249:end) = [];
    data_bpfreq_res_sel.grad.chanori(249:end, :) = [];
    data_bpfreq_res_sel.grad.chanpos(249:end, :) = [];
    data_bpfreq_res_sel.grad.tra(249:end, :) = [];
    data_bpfreq_res_sel.label(249:end) = [];
end

fn_data_bpfreq_sel_res = fieldnames(data_bpfreq_res_sel);
fn_filename = fieldnames(filename);

diff_fieldnames = setdiff(fn_filename, fn_data_bpfreq_sel_res);

for k=1:length(diff_fieldnames)
    data_bpfreq_res_sel.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
end
fn_filename{end+1}='cfg';
data_bpfreq_res_sel = orderfields(data_bpfreq_res_sel, fn_filename);
clearvars filename data_bpfreq data_bpfreq_res 

end
     


