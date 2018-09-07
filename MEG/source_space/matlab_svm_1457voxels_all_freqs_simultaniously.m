
 function []= main(path2data, outPath_extdisc, freqname, condition)
 sessions = [];
 [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, condition);
 [sessions] = concatenate_sessions(sessions, freqname);
 adi_matlab_svm (sessions, outPath_extdisc)
 
%   save ([outPath 'MEG\sourcespace\runs_appended\virtsens\' cfg_virtsens '_' condition '_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');
 end
 
 
 function [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, condition)
 
 
 for k = 1:length(freqname)
     fileName = [outPath_extdisc 'MEG\sourcespace\noROIs\runs_appended\virtsens\virtsens_all_conditions_allRuns_', freqname{k}, '.mat'];
     if ~exist(fileName, 'file')
        [vs_allRuns] = adi_appendvirtsensMEG(path2data, outPath_extdisc, freqname{k}, condition);
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
    if strcmp((freqname{k}), 'bp1-45Hz')
        sessions.bp1_45Hz = session;
    else
        sessions.(freqname{k}) = session;
    end
    clear session vs_allRuns
 end
 

 
end
  


     
 %%
 function [] = adi_matlab_svm (sessions, outPath_extdisc)
 


um_conditions = length(unique(session.labels));
data = sessions.data(:,:,150);
data = data *10^11;
data(:,451)=labels';
isfloat(data)
data_table = array2table(data);
data_table(:,91:382)=[];
labels = sessions.labels;
table_labels = array2table(labels')
data_table(:, 451)= table_labels;
switch num_conditions
    case 3
        % like vs dislike:
        ind_like = find(session.labels == 1);
        ind_dislike = find(session.labels == 2);
        ind_dontcare = find(session.labels == 3);
        X = session.data([ind_like ind_dislike],:,:);
        y = session.labels([ind_like ind_dislike]);
        rand_num = randperm(size(X,1));
        X_train = X(rand_num(1:floor((0.8*size(X,1)))),:,:);
        y_train = y(rand_num(1:floor((0.8*size(X,1)))));

        X_test = X(rand_num(floor((0.8*size(X,1)))+1:end),:,:);
        y_test = y(rand_num(floor((0.8*size(X,1)))+1:end));
        c = cvpartition(y_train, 'k', 5);
    case 2
end

for i = 1:size(session.data,3)
 X_train_sample = X_train(:,:,i);
 X_test_sample = X_test(:,:,i);

opts = statset('display', 'iter');
fun = @(train_data, train_labels, test_data,  test_labels)...
sum(predict(fitcsvm(train_data, train_labels, 'KernelFunction', 'rbf'), test_data) ~= test_labels)
[fs, history] = sequentialfs(fun, X_train_sample, y_train', 'cv', c, 'options', opts, 'nfeatures', 5);

%% Best hyperparameters => which features are useful for training and the classification 
X_train_w_best_features = X_train(:,fs);   
Md1 = fitcsvm(X_train_w_best_features, y_train, 'KernelFunction', 'rbf', 'OptimizeHyperparameters', 'all');
 

Md1 = fitcsvm(X_train_w_best_features, y_train, 'KernelFunction', 'rbf', 'OptimizeHyperparameters', 'auto', ...
    'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', 'ShowPlots', true));
    %                         cv = cvpartition(200,'KFold',10);
%                         opts = struct('Optimizer','bayesopt','ShowPlots',true,'CVPartition',cv,...
%                         'AcquisitionFunctionName','expected-improvement-plus');
%                         model_matlab_svmtoolbox = fitcsvm(data_train, y_train,'KernelFunction','rbf',...
%                         'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts)
%                         model_matlab_svmtoolbox = fitcsvm(data_train, y_train, 'CacheSize', 'maximal', 'CrossVal', 'on', 'KernelFunction', 'linear'); % Matlab machine learning toolbox: Default: 'linear' for two-class
% %                       learning and 'gaussian' (or 'rbf') for one-class learning
%% test

X_test_w_best_feature = X_test(:, fs);
accuracy = sum(predict(Md1, X_test_w_best_feature) == y_test)/length(y_test)*100;

%%
figure;

hgscatter = gscatter(X_train_w_best_features(:,1), X_train_w_best_features(:,2), y_train);
hold on;
j_sv = plot(Md1.SupportVectors(:,1), Md1.SupportVectors(:,2), 'ko', 'markersize', 8) ;
 
 
 end
 end
%% 

function [sessions_all_freqs] = concatenate_sessions(sessions, freqname)
if strcmp((freqname{1}), 'bp1-45Hz')
    sessions_all_freqs.data = sessions.bp1_45Hz.data;
    sessions_all_freqs.labels =  sessions.bp1_45Hz.labels;
    sessions_all_freqs.time =  sessions.bp1_45Hz.time{1,1};
    
else
    sessions_all_freqs.data = sessions.(freqname{1}).data;
    sessions_all_freqs.labels = sessions.(freqname{1}).labels;
    sessions_all_freqs.labels =  sessions.bp1_45Hz.time{1,1};
end 

for i = 1:length(freqname)
    for k = 1:size(sessions.(freqname{2}).data,2)
        freqs{k,i} = (freqname{i});
    end
end
freqs = reshape(freqs, (size(freqs,1)*size(freqs,2)), 1);
for i = 1:length(freqname)
    for k = 1:size(sessions.(freqname{2}).data,2)
        voxels{k,i} = k;
    end
end
voxels =  reshape(voxels, (size(freqs,1)*size(freqs,2)), 1); %repmat(1:size(sessions.(freqname{2}).data,2),1 ,5);
sessions_all_freqs.freq_voxels = [voxels freqs];

for k = 2:length(freqname)
    sessions_all_freqs.data = cat(2, sessions_all_freqs.data, sessions.(freqname{k}).data);
end

% save('sessions_all_freqs', 'sessions', '-v7.3')

end

 %%
 
 function [sigma] = kh_sigma_memory_friendly(sigma_conditions, sigma_, n_conditions, n_time, n_sensors, Xpseudo_train)
 
time_sample_memory_friendly = 5;

poolobj = parpool('local',4); 
 parfor c = 1:n_conditions
  % compute sigma for each time point, then average across time
  sigma_temp = nan(7, n_sensors, n_sensors);
    p = 1;
    i = 1;
    a = 1;
    for k = time_sample_memory_friendly : time_sample_memory_friendly : n_time
        tmp_ = nan(time_sample_memory_friendly, n_sensors, n_sensors);
        z = 1;
        for t = p:k
            tmp_(z, :, :) = cov1para(Xpseudo_train(sigma_conditions == c, :, t));
            z = z+1;
        end
        temp_mean_1(i, :, :) = mean(tmp_, 1);
        if i == 11
            sigma_temp(a,:,:) = mean(temp_mean_1, 1); 
            clear temp_mean_1
            temp_mean_1 = nan(11, n_sensors, n_sensors);
            a = a+1;
            i = 1;
        else
            i = i+1;
        end
        clear tmp_

        p = k + 1; 
    end

     sigma_(c, :, :) = mean(sigma_temp, 1);
     clear sigma_temp
 end
 
delete(gcp('nocreate'))
sigma = squeeze(mean(sigma_, 1));
          
 end
 
function [vs_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, cfg_virtsens, condition)

% append Runs 
% Channels   --> Configuration of Brainstorm MEG Channels    
if ~exist([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'file')
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
            sanity_path = [outPath 'MEG\sourcespace\noROIs\sanity_check\'];
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
    
    pathAppended = [outPath '\MEG\sourcespace\noROIs\runs_appended\virtsens\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end
   
    save ([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');

else
    vs_allRuns = [];
    %load ([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'vs_allRuns');
end

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
     
