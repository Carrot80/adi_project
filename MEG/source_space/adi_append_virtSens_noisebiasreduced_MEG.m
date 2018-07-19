
function  [vs_ns_dislike_allRuns, vs_ns_like_allRuns] = adi_append_virtSensMEG_noisebiasreduced(inPath, outPath, freqname)
% wird nicht mehr verwendet, kann gelöscht werden
 [vs_ns_like_allRuns] = append_virtSensMEG_noisebiasreduced(inPath, outPath, freqname, 'like')
 [vs_ns_dislike_allRuns] = append_virtSensMEG_noisebiasreduced(inPath, outPath, freqname, 'dislike')

end


function [vs_ns_allRuns] = append_virtSensMEG_noisebiasreduced(inPath, outPath, freqname, condition)
% append Runs 
% Channels   --> Configuration of Brainstorm MEG Channels    
    
vs_ns_allRuns = [];

if ~exist([outPath '\runs_appended\virtsens\vs_ns_' condition '_allRuns_', freqname, '.mat'], 'file');
    fileName = ([condition '*' freqname '.mat']);
    files = dir(fullfile(inPath, fileName));
    size_files = size(files);

    for i = 1:(size_files(1,1))
        load ([inPath files(i).name]);
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
                load ([outPath 'sourcespace\run1\spatialfilter_fixed_orientation_singletrials_' condition '_' freqname '.mat']);
            case 2 
                load ([outPath 'sourcespace\run2\spatialfilter_fixed_orientation_singletrials_' condition '_' freqname '.mat']);
            case 3
                load ([outPath 'sourcespace\run3\spatialfilter_fixed_orientation_singletrials_' condition '_' freqname '.mat']);
        end
        spatialfilter = cat(1, spatialfilter_orig{:});        
        virtsens = [];
        for k = 1:length(data_bpfreq.trial)
            virtsens.trial{k} = spatialfilter*data_bpfreq.trial{k};
        end;
        virtsens.time = data_bpfreq.time;
        virtsens.fsample = data_bpfreq.fsample;

        % siehe Skript von Yuval: => all will have a similar noise level
        ns = mean(abs(spatialfilter),2);

        for k = 1:length(virtsens.trial)
            virtsens_ns.trial{k} = virtsens.trial{1,k}./repmat(ns,1,size(virtsens.trial{1,k},2)); % => all will have a similar noise level
        end

        virtsens_ns.time = virtsens.time;
        virtsens_ns.fsample = virtsens.fsample;

        for k = 1:length(virtsens_ns.trial{1}(:,1))
            virtsens_ns.label{k} = num2str(k);
        end
               
        % timesamples reduzieren, um dateigröße ebenfalls zu reduzieren:
        cfg = [];
        cfg.latency = [-0.5 1]; 
        virtsens = ft_selectdata(cfg, virtsens);

        files2append(:,i) = virtsens;
        clear data_bpfreq spatialfilter spatialfilter_orig
    end


    switch size_files(1)
        case 3
            vs_ns_allRuns.trial = [files2append(1).trial files2append(2).trial files2append(3).trial];
            vs_ns_allRuns.time = [files2append(1).time files2append(2).time files2append(3).time];
        case 2
            vs_ns_allRuns.trial = [files2append(1).trial files2append(2).trial];
            vs_ns_allRuns.time = [files2append(1).time files2append(2).time];
    end
    
    pathAppended = [outPath '\runs_appended\virtsens\'];
    if ~exist (pathAppended, 'dir')
        mkdir (pathAppended)
    end
    save ([outPath '\runs_appended\virtsens\vs_ns_' condition '_allRuns_', freqname, '.mat'], 'vs_ns_allRuns');
    clear files2append virtsens
end

end