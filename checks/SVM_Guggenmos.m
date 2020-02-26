
 function [virtsens_allRuns] = main(path2data, outPath, outPath_extdisc, freqname, like, dislike, dontcare, cfg_virtsens)

 sessions(1) = struct('data', [], 'labels', []);
 sessions(2) = struct('data', [], 'labels', []);
 sessions(3) = struct('data', [], 'labels', []);
 
 [sessions] = mk_SVM_struct(sessions, path2data, outPath, outPath_extdisc, freqname, like, dislike, dontcare, num2str(1), cfg_virtsens)
 [sessions] = mk_SVM_struct(sessions, path2data, outPath, outPath_extdisc, freqname, like, dislike, dontcare, num2str(2), cfg_virtsens)
 [sessions] = mk_SVM_struct(sessions, path2data, outPath, outPath_extdisc, freqname, like, dislike, dontcare, num2str(3), cfg_virtsens)
 
 SVM_Guggenmos_tutorial(sessions)
 
%   save ([outPath 'MEG\sourcespace\runs_appended\virtsens\' cfg_virtsens '_' condition '_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');
 end
 
 
 function [sessions] = mk_SVM_struct(sessions, path2data, outPath, outPath_extdisc, freqname, like, dislike, dontcare, run, cfg_virtsens)
  
  fileName = (['*' run '*.mat']);
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

        switch run
            case '1'
                load ([outPath_extdisc 'MEG\sourcespace\run1\spatialfilter_loose_orientation_singletrials_' files(i).name(1:end-9) '_' freqname '.mat']);
            case '2' 
                load ([outPath_extdisc 'MEG\sourcespace\run2\spatialfilter_loose_orientation_singletrials_' files(i).name(1:end-9) '_' freqname '.mat']);
            case '3'
                load ([outPath_extdisc 'MEG\sourcespace\run3\spatialfilter_loose_orientation_singletrials_' files(i).name(1:end-9) '_' freqname '.mat']);
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
        
        % timesamples reduzieren, um Dateigröße ebenfalls zu reduzieren:
        cfg = [];
        cfg.latency = [-0.5 1]; 
        virtsens = ft_selectdata(cfg, virtsens);
        size_sessions_data = size(sessions(str2double(run)).data, 1);
        trl_num = length(virtsens.trial);
        start = 1+size_sessions_data;
        end_ = start + trl_num -1;
        for n=start:end_
            sessions(str2double(run)).data(n,:,:) = virtsens.trial{1,n-start+1};
        end
        
        switch files(i).name(1:end-9)
            case 'dislike'
                condition = 2;
            case 'like'
                condition = 1;
            case 'dontcare'
                condition = 3;
        end
        
        sessions(str2double(run)).labels(start:end_) = condition.* ones(1, length(virtsens.trial));
        clear spatialfilter_orig spatialfilter virtsens data_bpfreq cleanMEG_interp
  
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
     
 %%
function [] = SVM_Guggenmos_tutorial(sessions)
% We set a seed, in order to make analyses reproducible:
rng(10);

% Now we set some parameters. Only the number of permutations and 
% the number of pseudo-trials are free parameters. The number of conditions, 
% sensors, time points and sessions are derived from the data (i.e., from the sessions variable above).

% Parameters
n_perm = 20;  % number of permutations
n_pseudo = 5;  % number of pseudo-trials
n_conditions = length(unique(sessions(1).labels));
n_sensors = size(sessions(1).data, 2);
n_time = size(sessions(1).data, 3);
n_sessions = length(sessions);

% The analytic logic is contained in a nested for loop, with loops for the number of sessions, 
% number of permutations, number of timepoints, number of conditions, and number of conditions again. 
% Overall, the logic contains 4 crucial steps:

% 1. Compute pseudo-trials for the training and test data
% 2. Whiten the training data (here using the Epoch method, which is recommended in our manuscript)
% 3. Fit the classifier to the training data
% 4. Compute classification accuracy on test data
% 

% pre-load mechanism, for convenience
% preload_result = true; % for recomputing the decoding analyses, set to false
% if preload_result
%     load(fullfile(pwd, 'result_decoding.mat'))
% else
    % We define three classifiers that will be compared: Support Vector Machine, Weighted Robust Distance,
    % Gaussian Naive Bayes
    clfs = {'svm', 'gnb', 'weird'};
    for c = 1:length(clfs)
        result.(clfs{c}) = nan(n_sessions, n_perm, n_conditions, n_conditions, n_time);
    end
    for s = 1:n_sessions

        fprintf('sessions %g / %g\n', s, n_sessions)

        X = sessions(s).data;
        y = sessions(s).labels;
        
        conditions = unique(y);
        n_trials = histc(y, conditions);

        for f = 1:n_perm
            fprintf('\tPermutation %g / %g\n', f, n_perm)
            
            % precompute permutations
            ind_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
            ind_pseudo_test = nan(n_conditions, n_conditions, 2);
            labels_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
            labels_pseudo_test = nan(n_conditions, n_conditions, 2);
            for c1 = 1:n_conditions
                range_c1 = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
                for c2 = 1:n_conditions
                    range_c2 = (c2-1)*(n_pseudo-1)+1:c2*(n_pseudo-1);
                    ind_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = [range_c1 range_c2];
                    ind_pseudo_test(c1, c2, :) = [c1 c2];
                    labels_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = ...
                        [conditions(c1)*ones(1,n_pseudo-1) conditions(c2)*ones(1,n_pseudo-1)];
                    labels_pseudo_test(c1, c2, :) = conditions([c1 c2]);
                end
            end              
            train_indices = cell(1, n_conditions*(n_pseudo-1));
            test_indices = cell(1, n_conditions);
            for c1 = 1:n_conditions  % separate permutation for each condition
                prm_ = randperm(n_trials(c1));                
                prm = cell(1, n_pseudo);
                splitsize = n_trials(c1) / n_pseudo;
                for i = 1:n_pseudo
                    idxs = floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1;
                    prm{i} = prm_(idxs + 1);
                end                                
                ind = cellfun(@(x)x+sum(n_trials(1:c1-1)), prm, 'UniformOutput', 0);
                xrange = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
                for i = 1:length(xrange)
                    train_indices{xrange(i)} = ind{i}; % train und test trial indices für einzelne Bedingung und run
                end
                test_indices{c1} = ind{end};
            end                                

            % 1. Compute pseudo-trials for training and test
            Xpseudo_train = nan(length(train_indices), n_sensors, n_time);
            Xpseudo_test = nan(length(test_indices), n_sensors, n_time);
            for i = 1:length(train_indices)
                Xpseudo_train(i, :, :) = mean(X(train_indices{i}, :, :), 1);
            end
            for i = 1:length(test_indices)
                Xpseudo_test(i, :, :) = mean(X(test_indices{i}, :, :), 1);
            end


            % 2. Whitening using the Epoch method
            sigma_conditions = reshape(squeeze(labels_pseudo_train(1,:,n_pseudo:end))',1,[]);
            sigma_ = nan(n_conditions, n_sensors, n_sensors);
            for c = 1:n_conditions
                % compute sigma for each time point, then average across time
                tmp_ = nan(n_time, n_sensors, n_sensors);
                for t = 1:n_time
                    tmp_(t, :, :) = cov1para(Xpseudo_train(sigma_conditions==c, :, t));
                end
                sigma_(c, :, :) = mean(tmp_, 1);
            end
            sigma = squeeze(mean(sigma_, 1));  % average across conditions
            sigma_inv = sigma^-0.5;
            for t = 1:n_time
                Xpseudo_train(:, :, t) = squeeze(Xpseudo_train(:, :, t)) * sigma_inv;
                Xpseudo_test(:, :, t) = squeeze(Xpseudo_test(:, :, t)) * sigma_inv;
            end

            for t = 1:n_time
                for c1 = 1:n_conditions-1
                    for c2 = c1+1:n_conditions
                        % 3. Fit the classifier using training data
                        data_train = Xpseudo_train(ind_pseudo_train(c1, c2, :), :, t);
                        y_train = squeeze(labels_pseudo_train(c1, c2, :));
                        model_svm = svmtrain(y_train, data_train, '-c 1 -q 0 -t 0');
                        model_weird = weirdtrain(y_train, data_train);
                        model_gnb = gnbtrain(y_train, data_train);

                        % 4. Compute and store classification accuracies
                        data_test = Xpseudo_test(ind_pseudo_test(c1, c2, :), :, t);
                        y_train = squeeze(labels_pseudo_test(c1, c2, :));
                        result.svm(s, f, c1, c2, t) = ...
                            mean(svmpredict(y_train,data_test,model_svm,'-q 0 -t 0')==y_train)-0.5;
                        result.weird(s, f, c1, c2, t) = ...
                            mean(weirdpredict(y_train,data_test,model_weird)==y_train)-0.5;
                        result.gnb(s, f, c1, c2, t) = ...
                            mean(gnbpredict(y_train,data_test,model_gnb)==y_train)-0.5;
                    end
                end
            end
        end
    end
    % average across permutations
    for c = 1:length(clfs)
        result_.(clfs{c}) = nan(n_sessions, n_perm, n_conditions, n_conditions, n_time);
    end 
    result_.svm = squeeze(nanmean(result.svm, 2));
    result_.gnb = squeeze(nanmean(result.gnb, 2));
    result_.weird = squeeze(nanmean(result.weird, 2));
    result = result_;
    save(fullfile(pwd, 'result_decoding.mat'), 'result')
% end

% Now we plot the average classification accuracy time course by collapsing across sessions and conditions:

hold on
plot(-100:10:1001, 100*squeeze(nanmean(nanmean(nanmean(result.svm, 1), 2), 3))+50, 'linewidth', 2)
plot(-100:10:1001, 100*squeeze(nanmean(nanmean(nanmean(result.weird, 1), 2), 3))+50, 'linewidth', 2)
plot(-100:10:1001, 100*squeeze(nanmean(nanmean(nanmean(result.gnb, 1), 2), 3))+50, 'linewidth', 2)
plot([-100 1000], [50 50], 'k-')
xlim([-100 1000])
xlabel('Time [ms]')
ylabel('Classification accuracy')
legend('SVM', 'WeiRD', 'GNB')

% Already for one participant and a reduced data set (10 insteada of 92 conditions), these results look 
% like canonical decoding time courses. Note that the period -100ms to 0ms is the baseline phase and stimulus onset is at 0ms. In this example, Support Vector Machine and WeiRD outperform Gaussian Naive Bayes.

 end
 %%
 

