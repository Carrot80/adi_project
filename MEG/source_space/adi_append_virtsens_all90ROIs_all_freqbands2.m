
 function []= main(path2data, outPath, freqname, condition)

 [session] = mk_SVM_struct([], path2data, outPath_extdisc, freqname, condition);

 
%   save ([outPath 'MEG\sourcespace\runs_appended\virtsens\' cfg_virtsens '_' condition '_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');
 end
 
 
 function [session] = mk_SVM_struct(session, path2data, outPath, freqname, condition)
 
 fileName = [outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\virtsens_all_conditions_allRuns_', freqname, '.mat'];
 if ~exist(fileName, 'file')
    [vs_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, condition);
 else 
    load (fileName)
 end

 
num_conditions = size(fields(vs_allRuns),1);
fieldnames =   fields(vs_allRuns);       

if 3 == num_conditions
    data = cat(2, vs_allRuns.(condition{1}).trial, vs_allRuns.(condition{2}).trial, vs_allRuns.(condition{3}).trial);
    session.labels = cat(2, ones(1,length(vs_allRuns.(condition{1}).trial)), 2*ones(1, length(vs_allRuns.(condition{2}).trial)), 3*ones(1, length(vs_allRuns.(condition{3}).trial)));
elseif 2 == num_conditions
    data = cat(2, vs_allRuns.(fieldnames{1}).trial, vs_allRuns.(fieldnames{2}).trial);
    switch (fieldnames{1})
        case 'like'
            ind_field_1 = 1;
        case 'dislike' 
            ind_field_1 = 2;
        case 'dontcare'
            ind_field_1 = 3;
    end
    switch (fieldnames{2})
        case 'like'
            ind_field_2 = 1;
        case 'dislike' 
            ind_field_2 = 2;
        case 'dontcare'
            ind_field_2 = 3;
    end       
    session.labels = cat(2, ind_field_1*ones(1, length(vs_allRuns.(fieldnames{1}).trial)), ind_field_2*ones(1, length(vs_allRuns.(fieldnames{2}).trial)));
end
for i = 1:length(data)
    temp = data{1,i};
    session.data(i,:,:) = temp;
end
session.time = vs_allRuns.(fieldnames{1}).time;

end
  


function [vs_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, condition, atlas)

% append Runs 
% Channels   --> Configuration of Brainstorm MEG Channels    
if ~exist([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\virtsens_all_conditions_allRuns_', freqname, '.mat'], 'file')
    vs_allRuns = [];
    spatialfilter = {};

    for j = 1:size(condition,2)
        fileName = ([condition{j} '*.mat']);
        files = dir(fullfile(path2data, fileName));
        size_files = size(files);

        for i = 1:size_files(1,1)
            load ([path2data files(i).name]);

            [data_bpfreq] = adi_bpfilter(cleanMEG_interp, freqname);
         
            if ~isfield(spatialfilter, (['run' num2str(i)]))
                spatialfilter.(['run' num2str(i)]) = load ([outPath 'MEG\sourcespace\spatialfilter\run' num2str(i) '\spatialfilter_loose_orientation_singletrials_' freqname '.mat']);
                spatialfilter.(['run' num2str(i)]) = cat(1,spatialfilter.(['run' num2str(i)]).spatialfilter_orig{:});   
            end

                
            virtsens = [];
            for k = 1:length(data_bpfreq.trial)
                virtsens.trial{k} = spatialfilter.(['run' num2str(i)])*data_bpfreq.trial{k};
            end
        
            virtsens.time = data_bpfreq.time;
            virtsens.fsample = data_bpfreq.fsample;
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
            
            % sanity check:
            cfg=[];
            avg = ft_timelockanalysis(cfg, virtsens);
            figure;
            plot(virtsens.time{1,1}, mean(avg.avg))
            sanity_path = [outPath 'MEG\sourcespace\all_voxels_without_cerebellum\sanity_check\'];
            if ~exist(sanity_path, 'dir')
                mkdir (sanity_path)
            end
            savefig([sanity_path 'mean_virtsens_' files(i).name '_' freqname '.fig'])
            close
            
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

            files2append.(condition{j}).(['run' num2str(i)]) = virtsens;
            clear virtsens data_bpfreq cleanMEG_interp vs_euclid_norm

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
    
    % sanity check Nr. 2:
   
    cfg = [];
    for m = 1:size(fields(vs_allRuns),1)
        avg = ft_timelockanalysis(cfg, vs_allRuns.(condition{m}));
        figure;
        plot(avg.time, mean(avg.avg))
        savefig([sanity_path 'mean_virtsens_allRuns_' (condition{m}) '_' freqname '.fig'])
        close
    end
    
    pathAppended = [outPath '\MEG\sourcespace\all_voxels_without_cerebellum\runs_appended\virtsens\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end
   
    save ([outPath 'MEG\sourcespace\all_voxels_without_cerebellum\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');

else
    vs_allRuns = [];
    %load ([outPath 'MEG\sourcespace\all_voxels_without_cerebellum\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'vs_allRuns');
end


% session.data = session.data(:, 92:end,:);
% session.atlas = atlas;
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
if isfield(filename, 'trialinfo')
    data_bpfreq.trialinfo = filename.trialinfo;
end

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

       
end
     
