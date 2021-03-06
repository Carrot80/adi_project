
function  [virtsens_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, cfg_virtsens, condition)


 [virtsens_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, cfg_virtsens, condition);
 
 
end

function [vs_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, cfg_virtsens, condition)

% append Runs 
% Channels   --> Configuration of Brainstorm MEG Channels    
% if ~exist([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'file')
    vs_allRuns = [];
    source_avg_appended_conditions = {};

    for j = 1:size(condition,2)
        fileName = ([condition{j} '*.mat']);
        files = dir(fullfile(path2data, fileName));
        size_files = size(files);

        for i = 1:size_files(1,1)
            load ([path2data files(i).name]);

            [data_bpfreq] = adi_bpfilter(cleanMEG_interp, freqname);
         
            if ~isfield(source_avg_appended_conditions, (['run' num2str(i)]))
                source_avg_appended_conditions.(['run' num2str(i)]) = load ([outPath 'MEG\sourcespace\spatialfilter\run' num2str(i) '\source_avg_appended_conditions_' freqname '.mat']);
                source_avg_appended_conditions.(['run' num2str(i)]).spatialfilter = cat(1,source_avg_appended_conditions.(['run' num2str(i)]).source_avg.avg.filter{:});   
            end

            % laut fieldtrip discussion list zuerst weights mit sensor data multiplizieren: https://mailman.science.ru.nl/pipermail/fieldtrip/2017-July/011661.html    
            virtsens = [];
            for k = 1:length(data_bpfreq.trial)
                virtsens.trial{k} = source_avg_appended_conditions.(['run' num2str(i)]).spatialfilter*data_bpfreq.trial{k};
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
            
            virtsens.label = cell(length(virtsens.trial{1,1}(:,1)), 1);
            for n = 1:length(virtsens.trial{1,1}(:,1))
                virtsens.label{n,1}=num2str(n);
            end
            
            
            fn_data_bpfreq = fields(data_bpfreq);
            fn_virtsens = fields(virtsens);
            diff_fieldnames = setdiff(fn_data_bpfreq, fn_virtsens);

            for k=1:length(diff_fieldnames)
                virtsens.(diff_fieldnames{k}) = data_bpfreq.(diff_fieldnames{k}); 
            end
             
            % sanity check:
            cfg = [];
            avg = ft_timelockanalysis(cfg, virtsens);
            figure;
            plot(virtsens.time{1,1}, avg.avg)
            sanity_path = [outPath 'MEG\sourcespace\noROIs\sanity_check\'];
            if ~exist(sanity_path, 'dir')
                mkdir (sanity_path)
            end
            savefig([sanity_path 'virtsens_' files(i).name '_' freqname '.fig'])
            close

        %% noise normalization nach Gespr�ch mit Stefan:
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
        avg_ns = ft_timelockanalysis(cfg, virtsens_ns);
        figure;
        plot(virtsens_ns.time{1,1}, avg_ns.avg)
        savefig([sanity_path 'virtsens_ns_all_conditions_run' num2str(run) '_' freq '.fig'])
        close
  
%%
            files2append.(condition{j}).(['run' num2str(i)]) = virtsens_ns;
            clear virtsens_ns data_bpfreq cleanMEG_interp vs_euclid_norm

        end
    end 
    condition = fields(files2append);
    
    for k = 1:size(fields(files2append),1)
        switch size(fields(files2append.(condition{k})),1)
            case 3 
                vs_allRuns.(condition{k}).trial = [files2append.(condition{k}).(['run' num2str(1)]).trial files2append.(condition{k}).(['run' num2str(2)]).trial files2append.(condition{k}).(['run' num2str(3)]).trial];
                vs_allRuns.(condition{k}).time = [files2append.(condition{k}).(['run' num2str(1)]).time files2append.(condition{k}).(['run' num2str(2)]).time files2append.(condition{k}).(['run' num2str(3)]).time];
                vs_allRuns.(condition{k}).label =  files2append.(condition{k}).(['run' num2str(1)]).label;
            case 2
                vs_allRuns.(condition{k}).trial = [files2append.(condition{k}).(['run' num2str(1)]).trial files2append.(condition{k}).(['run' num2str(2)]).trial];
                vs_allRuns.(condition{k}).time = [files2append.(condition{k}).(['run' num2str(1)]).time files2append.(condition{k}).(['run' num2str(2)]).time];
                vs_allRuns.(condition{k}).label =  files2append.(condition{k}).(['run' num2str(1)]).label;
        end
    end
    
    % sanity check Nr. 3:
   
    cfg = [];
    for m = 1:size(fields(vs_allRuns),1)
        avg = ft_timelockanalysis(cfg, vs_allRuns.(condition{m}));
        figure;
        plot(avg.time, mean(avg.avg))
        savefig([sanity_path 'mean_virtsens_allRuns_' (condition{m}) '_' freqname '.fig'])
        close
    end
    
    pathAppended = [outPath '\MEG\sourcespace\noROIs\runs_appended\virtsens\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end
   
    save ([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');

% else
%     vs_allRuns = [];
    %load ([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'vs_allRuns');
% end

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


for k = 1:length(data_bpfreq.trial)
    data_bpfreq.trial{1,k} = data_bpfreq.trial{1,k}(1:248,:);
end    

for k = 1:length(data_bpfreq.trial)
    data_bpfreq.grad.label(249:end) = [];
    data_bpfreq.grad.chanori(249:end, :) = [];
    data_bpfreq.grad.chanpos(249:end, :) = [];
    data_bpfreq.grad.tra(249:end, :) = [];
    data_bpfreq.label(249:end) = [];
end

for m = 1:length(data_bpfreq.trial)
    if 3053 == length(data_bpfreq.time{1,m}) 
        data_bpfreq.trial{1,m}(:,1) = [];
        data_bpfreq.time{1,m}(:,1) = [];
    end
end
cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);
if isfield(filename, 'trialinfo')
    data_bpfreq_res_sel.trialinfo = filename.trialinfo;
end
       
end
     
