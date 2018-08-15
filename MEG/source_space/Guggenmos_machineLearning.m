
 function []= main(path2data, outPath_extdisc, freqname, like, dislike, dontcare, meanmax)

 [sessions] = mk_SVM_struct([], path2data, outPath_extdisc, freqname, like, dislike, dontcare, num2str(1), meanmax)
 [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, like, dislike, dontcare, num2str(2), meanmax)
 [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, like, dislike, dontcare, num2str(3), meanmax)
 [sessions] = kh_concatenate_runs(sessions)
 SVM_Guggenmos_tutorial(sessions, outPath_extdisc)
 
%   save ([outPath 'MEG\sourcespace\runs_appended\virtsens\' cfg_virtsens '_' condition '_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');
 end
 
 
 function [sessions] = mk_SVM_struct(sessions, path2data, outPath_extdisc, freqname, like, dislike, dontcare, run, meanmax)
  
  fileName = (['vs_' meanmax '*' freqname '.mat']);
  files = dir(fullfile([path2data 'run' run filesep fileName]));
  size_files = size(files);
    
  for i = 1:(size_files(1,1))
      virtsens = load ([path2data 'run' run filesep files(i).name]);
        
        if 1 == contains(files(i).name, '_dislike_')
            condition = 2;
            fn = [meanmax '_all_rois_dislike' ];
            for k = 1:size(virtsens.(fn).trial,2)
                for p = 1:size(virtsens.(fn).trial(:,k),1)
                    data(str2double(run)).dislike(k,p,:) = virtsens.(fn).trial{p,k};
                end 
            end
            labels(str2double(run)).dislike = 2*ones(1,size(virtsens.(fn).trial,2));
        elseif 1 == contains(files(i).name, '_like_')
            condition = 1;
            fn = [meanmax '_all_rois_like' ];
            for k = 1:size(virtsens.(fn).trial,2)
                for p = 1:size(virtsens.(fn).trial(:,k),1)
                    data(str2double(run)).like(k,p,:) = virtsens.(fn).trial{p,k};
                end 
            end
            labels(str2double(run)).like = ones(1,size(virtsens.(fn).trial,2));
             
        elseif 1 == contains(files(i).name, '_dontcare_')
            condition = 3;
             fn = [meanmax '_all_rois_dontcare' ];
             for k = 1:size(virtsens.(fn).trial,2)
                 for p = 1:size(virtsens.(fn).trial(:,k),1)
                    data(str2double(run)).dontcare(k,p,:) = virtsens.(fn).trial{p,k};
                 end 
             end
             labels(str2double(run)).dontcare = 3*ones(1,size(virtsens.(fn).trial,2));
        end 
          
%         size_sessions_data = size(sessions(str2double(run)).data, 1);
%         trl_num = size(virtsens.(fn).trial,2);
%         start = 1+size_sessions_data;
%         end_ = start + trl_num -1;
%         for n=start:end_
%             sessions(str2double(run)).data(n,:,:) = virtsens.(fn).trial{1,n-start+1};
%         end

  end
  
num_conditions = size(fields(data),1);
fieldnames =   fields(data);  
switch num_conditions
    case 2
        sessions(str2double(run)).data = cat(1, data.(fieldnames{1}), data.(fieldnames{2}));
        sessions(str2double(run)).labels = [labels.(fieldnames{1}) labels.(fieldnames{2})];
    case 3
        sessions(str2double(run)).data = cat(1, data.(fieldnames{1}), data.(fieldnames{2}), data.(fieldnames{3}));
        sessions(str2double(run)).labels = cat(2, labels.(fieldnames{1}), labels.(fieldnames{2}), labels.(fieldnames{3}));
end


%   sessions(str2double(run)).labels(start:end_) = condition.* ones(1, length(virtsens.trial));
  
   
end

function  [sessions] = kh_concatenate_runs(sessions)

sessions2.data = cat(1, sessions(1).data, sessions(2).data, sessions(2).data);
sessions2.labels = cat(2, sessions(1).labels, sessions(2).labels, sessions(2).labels);

end



     
 %%
function [] = SVM_Guggenmos_tutorial(sessions, outPath)
% We set a seed, in order to make analyses reproducible:
rng(10);

% Now we set some parameters. Only the number of permutations and 
% the number of pseudo-trials are free parameters. The number of conditions, 
% sensors, time points and sessions are derived from the data (i.e., from the sessions variable above).

% Parameters
n_perm = 20;  % number of permutations
n_pseudo = 5;  % number of pseudo-trials % evtl. erhöhen auf 10?
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
                    train_indices{xrange(i)} = ind{i}; % train indices: 80% aller trials beider Bedigungen pro run: die ersten 4 trials sind aus Bedingung A, die letzten 4 aus Bedingung B(jeweils 80%)
                end
                test_indices{c1} = ind{end}; % 20 % der trials beider Bedingungen gemischt
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
                        model_svm = svmtrain(y_train, data_train, '-c 1 -t 0'); % libSVM: -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"; -q : quiet mode (no outputs)\; "-t kernel_type : set type of kernel function (default 2)\n"   
                        model_svm_mltoolbox = fitcsvm(data_train, y_train, 'CacheSize', 'maximal', 'CrossVal', 'on', 'KernelFunction', 'linear'); % Matlab machine learning toolbox: Default: 'linear' for two-class
%                        learning and 'gaussian' (or 'rbf') for one-class learning
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
figure; axis tight
hold on
plot(max_all_rois_dislike.time, 100*squeeze(nanmean(nanmean(nanmean(result.svm, 1), 2), 3))+50, 'linewidth', 1) % nanmean = mean ignoring nans
plot(max_all_rois_dislike.time, 100*squeeze(nanmean(nanmean(nanmean(result.weird, 1), 2), 3))+50, 'linewidth', 1)
plot(max_all_rois_dislike.time, 100*squeeze(nanmean(nanmean(nanmean(result.gnb, 1), 2), 3))+50, 'linewidth', 1)
plot([-0.5 1], [50 50], 'k-')
xlim([-0.5 1])
xlabel('Time [s]')
ylabel('Classification accuracy')
legend('SVM', 'WeiRD', 'GNB')

% Already for one participant and a reduced data set (10 insteada of 92 conditions), these results look 
% like canonical decoding time courses. Note that the period -100ms to 0ms is the baseline phase and stimulus onset is at 0ms. In this example, Support Vector Machine and WeiRD outperform Gaussian Naive Bayes.

 end
 %%
 

