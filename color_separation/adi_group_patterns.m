
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
        try
            data_all_conditions.triggerlabel = [data(1).cleanMEG_interp.trialinfo.sampleinfo(:,5)' data(2).cleanMEG_interp.trialinfo.sampleinfo(:,5)' data(3).cleanMEG_interp.trialinfo.sampleinfo(:,5)'];
        catch
            data_all_conditions.triggerlabel = [data(1).cleanMEG_interp.trialinfo.triggerlabel data(2).cleanMEG_interp.trialinfo.triggerlabel data(3).cleanMEG_interp.trialinfo.triggerlabel];
        end
    case 2
        data_all_conditions.trial = [data(1).cleanMEG_interp.trial  data(2).cleanMEG_interp.trial];
        data_all_conditions.time = [data(1).cleanMEG_interp.time  data(2).cleanMEG_interp.time];
        try
            data_all_conditions.triggerlabel = [data(1).cleanMEG_interp.trialinfo.sampleinfo(:,5)' data(2).cleanMEG_interp.trialinfo.sampleinfo(:,5)'];
        catch
            data_all_conditions.triggerlabel = [data(1).cleanMEG_interp.trialinfo.triggerlabel data(2).cleanMEG_interp.trialinfo.triggerlabel];
        end
    case 1
        data_all_conditions.trial = data(1).cleanMEG_interp.trial;
        data_all_conditions.time = data(1).cleanMEG_interp.time;
        data_all_conditions.triggerlabel = data(1).cleanMEG_interp.trialinfo.sampleinfo(:,5)';
end
    data_all_conditions.label = data(1).cleanMEG_interp.label;
    data_all_conditions.fsample = data(1).cleanMEG_interp.fsample;
    data_all_conditions.dimord = data(1).cleanMEG_interp.dimord;
    
    [vs_allRuns] = adi_appendvirtsensMEG(vs_allRuns, data_all_conditions, path2sensordata, path2spatialfilter, dir_out, subject, pattern, trigger, freq, run);

end

function [vs_allRuns] = adi_appendvirtsensMEG(vs_allRuns, data_all_conditions, path2sensordata, path2spatialfilter, dir_out, subject, pattern, trigger, freq, run)


% if ~exist([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'file')

    % spatial filter wird geladen
    % data_all_runs werden aufgeteilt in Muster volley, space und
    % soccerball(=football)    

    [data_all_conditions] = regroup_trials(data_all_conditions, pattern, trigger);
    [data_all_conditions_bpfreq] = adi_bpfilter(data_all_conditions, freq);
    clearvars data_all_conditions     
    
    load ([path2spatialfilter 'run' num2str(run) '\source_avg_appended_conditions_' freq '.mat'], 'source_avg');
    spatialfilter = cat(1,source_avg.avg.filter{:});   


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
    
    fn_data_all_conditions_bpfreq = fields(data_all_conditions_bpfreq);
    fn_virtsens = fields(virtsens);
    diff_fieldnames = setdiff(fn_data_all_conditions_bpfreq, fn_virtsens);
    
    for k=1:length(diff_fieldnames)
        virtsens.(diff_fieldnames{k}) = data_all_conditions_bpfreq.(diff_fieldnames{k}); 
    end
  
    % sanity check:
    cfg = [];
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
    
    cfg = [];
    avg_ns = ft_timelockanalysis(cfg, virtsens_ns);
    figure;
    plot(virtsens.time{1,1}, avg_ns.avg)
    
    % sanity check nr. 2:
    cfg = [];
    avg_ns = ft_timelockanalysis(cfg, virtsens_ns);
    figure;
    plot(virtsens.time{1,1}, avg_ns.avg)
    savefig([sanity_path 'virtsens_ns_all_conditions_run' num2str(run) '_' freq '.fig'])
    close
    
    vs_allRuns.(['run' num2str(run)]) = virtsens_ns;
   

end

function [data_all_conditions] = regroup_trials(data_all_conditions, pattern, trigger) 

 % Volleyballdesign:
    for k = 1:size(pattern.volley.labels,1)
        pattern.volley.trigger(k) = trigger.triggerchannel(find(strcmp(trigger.labels_short,pattern.volley.labels(k))));% attern.(conditions{j}).labels; %(volleyball) = j=1;
        volley.trialnumbers{k,:} = find(data_all_conditions.triggerlabel == pattern.volley.trigger(k));
        pattern.volley.balldesign{k} = trigger.labels(find(strcmp(trigger.labels_short,pattern.volley.labels(k))));
        volley.balldesign{k,:} = find(data_all_conditions.triggerlabel == pattern.volley.trigger(k)); % hier weitermachen
    end
    
    volley.trialnumbers = sort(cat(2,volley.trialnumbers{:}));
    volley.trial = data_all_conditions.trial(volley.trialnumbers);
    volley.triggers = data_all_conditions.triggerlabel(volley.trialnumbers);
    volley.time = data_all_conditions.time(volley.trialnumbers);
    volley.ballpattern = cell(1, length(volley.trial));
    volley.ballpattern(:) = {'Volley'};
    volley.balldesign = cell(1, length(volley.trial));
    for k = 1:length(volley.triggers)
        volley.balldesign(k) = trigger.labels(find(trigger.triggerchannel == volley.triggers(k)));
    end
    
    %Spaceballdesign:
    for k = 1:size(pattern.space.labels,1)
        pattern.space.trigger(k) = trigger.triggerchannel(find(strcmp(trigger.labels_short,pattern.space.labels(k))));% attern.(conditions{j}).labels; %(spaceball) = j=1;
        space.trialnumbers{k,:} = find(data_all_conditions.triggerlabel == pattern.space.trigger(k));
        pattern.space.balldesign{k} = trigger.labels(find(strcmp(trigger.labels_short,pattern.space.labels(k))));
        space.balldesign{k,:} = find(data_all_conditions.triggerlabel == pattern.space.trigger(k)); % hier weitermachen
    end
    
    space.trialnumbers = sort(cat(2,space.trialnumbers{:}));
    space.trial = data_all_conditions.trial(space.trialnumbers);
    space.triggers = data_all_conditions.triggerlabel(space.trialnumbers);
    space.time = data_all_conditions.time(space.trialnumbers);
    space.ballpattern = cell(1, length(space.trial));
    space.ballpattern(:) = {'Space'};
    space.balldesign = cell(1, length(space.trial));
    for k = 1:length(space.triggers)
        space.balldesign(k) = trigger.labels(find(trigger.triggerchannel == space.triggers(k)));
    end
    
    
    %Soccerballdesign:
    for k = 1:size(pattern.soccer.labels,1)
        pattern.soccer.trigger(k) = trigger.triggerchannel(find(strcmp(trigger.labels_short,pattern.soccer.labels(k))));% attern.(conditions{j}).labels; %(soccerball) = j=1;
        soccer.trialnumbers{k,:} = find(data_all_conditions.triggerlabel == pattern.soccer.trigger(k));
        pattern.soccer.balldesign{k} = trigger.labels(find(strcmp(trigger.labels_short,pattern.soccer.labels(k))));
        soccer.balldesign{k,:} = find(data_all_conditions.triggerlabel == pattern.soccer.trigger(k)); % hier weitermachen
    end
    
    soccer.trialnumbers = sort(cat(2,soccer.trialnumbers{:}));
    soccer.trial = data_all_conditions.trial(soccer.trialnumbers);
    soccer.triggers = data_all_conditions.triggerlabel(soccer.trialnumbers);
    soccer.time = data_all_conditions.time(soccer.trialnumbers);
    soccer.ballpattern = cell(1, length(soccer.trial));
    soccer.ballpattern(:) = {'Soccer'};
    soccer.balldesign = cell(1, length(soccer.trial));
    for k = 1:length(soccer.triggers)
        soccer.balldesign(k) = trigger.labels(find(trigger.triggerchannel == soccer.triggers(k)));
    end
    
    fields2remove = {'trial', 'time', 'triggerlabel'};
    data_all_conditions = rmfield(data_all_conditions, fields2remove);
    
    fieldnames_all = fieldnames(volley);
    
    for k=1:size(fieldnames_all,1)
        data_all_conditions.(fieldnames_all{k}) = cat(2, volley.(fieldnames_all{k}), space.(fieldnames_all{k}), soccer.(fieldnames_all{k})) ; 
    end
     
end



function [data_bpfreq_sel_res] = adi_bpfilter(filename, bpname)


for k = 1:length(filename.trial)
    filename.trial{k}(249:end,:)=[];
end  

for k = 1:length(filename.trial)
%     data_bpfreq.grad.label(249:end) = [];
%     data_bpfreq.grad.chanori(249:end, :) = [];
%     data_bpfreq.grad.chanpos(249:end, :) = [];
%     data_bpfreq.grad.tra(249:end, :) = [];
    filename.label(249:end) = [];
end

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

for m = 1:length(filename.trial)
 if isnan(filename.trial{m}(1,1))
    [ind] = sum(isnan(filename.trial{m}(1,:)));
    filename.trial{m}(:,1:ind) = [];
    filename.time{m}(1:ind) = [];
 end
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

for m = 1:length(data_bpfreq.trial)
    if 3053 == length(data_bpfreq.time{1,m}) 
        data_bpfreq.trial{1,m}(:,1) = [];
        data_bpfreq.time{1,m}(:,1) = [];
    end
end
cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_sel = ft_selectdata(cfg, data_bpfreq);

cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'no';
[data_bpfreq_sel_res] = ft_resampledata(cfg, data_bpfreq_sel);

fn_data_bpfreq_sel_res = fieldnames(data_bpfreq_sel_res);
fn_filename = fieldnames(filename);

diff_fieldnames = setdiff(fn_filename, fn_data_bpfreq_sel_res);

for k=1:length(diff_fieldnames)
    data_bpfreq_sel_res.(diff_fieldnames{k}) = filename.(diff_fieldnames{k});
end

       
end
     
