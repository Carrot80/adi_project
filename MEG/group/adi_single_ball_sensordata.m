


function  [sensordata_ball_freq_allsubj] = append_sensordata(sensordata_ball_freq_allsubj, path2sensor_data, subject, pattern, trigger, freq, balldesign, delete_runs, subjnum)

dir_runs = dir(path2sensor_data);
dir_runs(1:2)=[];

for i = 1:length(delete_runs)
    if 1==strcmp(delete_runs{i,1}, subject)
        for k = 1:size(dir_runs,1)
            temp=contains(dir_runs(k).name, ['500_' delete_runs{i,2}]);
            ind_(k)=temp;
            clearvars temp
        end
            ind =find(ind_);
            dir_runs(ind)=[];
            clearvars ind ind_
    end
end

for k = 1:length(dir_runs)
    sensor_data(k)= load ([path2sensor_data  dir_runs(k).name], 'cleanMEG_interp');
end

sensordata_allruns = [];
sensordata_allruns.trial = sensor_data(1).cleanMEG_interp.trial;
sensordata_allruns.time = sensor_data(1).cleanMEG_interp.time;
sensordata_allruns.response_label = sensor_data(1).cleanMEG_interp.trialinfo.response_label;
sensordata_allruns.balldesign_short = sensor_data(1).cleanMEG_interp.trialinfo.balldesign_short;

for k=2:length(sensor_data)
    sensordata_allruns.trial = cat(2,sensordata_allruns.trial, sensor_data(k).cleanMEG_interp.trial);
    sensordata_allruns.time = cat(2,sensordata_allruns.time, sensor_data(k).cleanMEG_interp.time);
    sensordata_allruns.response_label = cat(2,sensordata_allruns.response_label, sensor_data(k).cleanMEG_interp.trialinfo.response_label);
    sensordata_allruns.balldesign_short = cat(2,sensordata_allruns.balldesign_short, sensor_data(k).cleanMEG_interp.trialinfo.balldesign_short); 
end

sensordata_allruns.label = sensor_data(1).cleanMEG_interp.label;
sensordata_allruns.fsample = sensor_data(1).cleanMEG_interp.fsample;
sensordata_allruns.grad = sensor_data(1).cleanMEG_interp.grad;
 
clearvars sensor_data
    
for k=1:length(sensordata_allruns.trial)
   temp(k)=strcmp(sensordata_allruns.balldesign_short{k}, balldesign);
end

ind=find(temp);

sensordata_allruns.trial =  sensordata_allruns.trial(ind);
sensordata_allruns.time =  sensordata_allruns.time(ind);
sensordata_allruns.response_label =  sensordata_allruns.response_label(ind);
sensordata_allruns.balldesign_short =  sensordata_allruns.balldesign_short(ind);

if ~isempty(sensordata_allruns.trial)
    [sensordata_allruns_ball_bpfreq] = adi_bpfilter(sensordata_allruns, freq);
else
    clearvars ind temp sensordata_allruns
    return
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
catch
    
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
     


