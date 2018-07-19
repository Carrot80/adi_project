
function  [virtsens_dislike_allRuns, virtsens_like_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, cfg_virtsens)


 [virtsens_like_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, 'like', cfg_virtsens)
 [virtsens_dislike_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, 'dislike', cfg_virtsens)
 
end

function [vs_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, condition, cfg_virtsens)

% append Runs 
% Channels   --> Configuration of Brainstorm MEG Channels    
vs_allRuns = [];

if ~exist([outPath '\MEG\sourcespace\runs_appended\virtsens\' cfg_virtsens '_' condition '_allRuns_', freqname, '.mat'], 'file')
    fileName = ([condition '*.mat']);
    files = dir(fullfile(path2data, fileName));
    size_files = size(files);
    
    for i = 1:(size_files(1,1))
        load ([path2data files(i).name]);
        
        [data_bpfreq] = adi_bpfilter(cleanMEG_interp, freqname);
        
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

        switch i
            case 1
                load ([outPath 'MEG\sourcespace\run1\spatialfilter_loose_orientation_singletrials_' condition '_' freqname '.mat']);
            case 2 
                load ([outPath 'MEG\sourcespace\run2\spatialfilter_loose_orientation_singletrials_' condition '_' freqname '.mat']);
            case 3
                load ([outPath 'MEG\sourcespace\run3\spatialfilter_loose_orientation_singletrials_' condition '_' freqname '.mat']);
        end
        spatialfilter = cat(1,spatialfilter_orig{:});        
        virtsens = [];
        for k = 1:length(data_bpfreq.trial)
            virtsens.trial{k} = spatialfilter*data_bpfreq.trial{k};
        end
        virtsens.time = data_bpfreq.time;
        virtsens.fsample = data_bpfreq.fsample;
        
         % siehe Skript von Yuval: => all will have a similar noise level
        
        if 1 == strcmp (cfg_virtsens, 'virtsens_ns')
            ns = mean(abs(spatialfilter),2);
            
            for k = 1:length(virtsens.trial)
                virtsens_ns.trial{k} = virtsens.trial{1,k}./repmat(ns,1,size(virtsens.trial{1,k},2)); % => all will have a similar noise level
            end

            virtsens_ns.time = virtsens.time;
            virtsens_ns.fsample = virtsens.fsample;

            for k = 1:length(virtsens_ns.trial{1}(:,1))
                virtsens_ns.label{k} = num2str(k);
            end
            virtsens = virtsens_ns;
            clear virtsens_ns
        end
        
        % timesamples reduzieren, um Dateigr��e ebenfalls zu reduzieren:
        cfg = [];
        cfg.latency = [-0.5 1]; 
        virtsens = ft_selectdata(cfg, virtsens);

        files2append(:,i) = virtsens;
        clear spatialfilter_orig spatialfilter virtsens data_bpfreq cleanMEG_interp

    end
   
    switch size_files(1)
        case 3
            vs_allRuns.trial = [files2append(1).trial files2append(2).trial files2append(3).trial];
            vs_allRuns.time = [files2append(1).time files2append(2).time files2append(3).time];
        case 2
            vs_allRuns.trial = [files2append(1).trial files2append(2).trial];
            vs_allRuns.time = [files2append(1).time files2append(2).time];
    end
    
    pathAppended = [outPath '\MEG\sourcespace\runs_appended\virtsens\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end
    save ([outPath 'MEG\sourcespace\runs_appended\virtsens\' cfg_virtsens '_' condition '_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');
end


end


function [data_bpfreq] = adi_bpfilter(filename, bpname)


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

       
end
     