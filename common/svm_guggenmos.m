 %%
function [] = SVM_Guggenmos_tutorial(sessions, outPath, freq, condition, atlas)

% We set a seed, in order to make analyses reproducible:
rng('default');
rng(10);
% 

for k=1:length(sessions.label)
    switch sessions.label{k}
        case 'like' % 'Volley'
            sessions.labels(k) = 1;
        case 'dislike' % 'Space'
            sessions.labels(k) = 2;
        case 'dontcare' % 'Soccer'
            sessions.labels(k) = 3;
    end
end

% for k=1:length(sessions.label)
%     switch sessions.label{k}
%         case 'Volley'
%             sessions.labels(k) = 1;
%         case 'Space'
%             sessions.labels(k) = 2;
%         case 'Soccer'
%             sessions.labels(k) = 3;
%     end
% end
% Now we set some parameters. Only the number of permutations and 
% the number of pseudo-trials are free parameters. The number of conditions, 
% sensors, time points and sessions are derived from the data (i.e., from the sessions variable above).

% Parameters
n_perm = 10;  % number of permutations
n_pseudo = 5;  % number of pseudo-trials % evtl. erhöhen auf 10?

n_conditions = length(unique(sessions(1).response_label));
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
%             sigma_conditions = reshape(squeeze(labels_pseudo_train(1,:,n_pseudo:end))',1,[]);
%             sigma_ = nan(n_conditions, n_sensors, n_sensors);
%             for c = 1:n_conditions
%                 % compute sigma for each time point, then average across time
%                 tmp_ = nan(n_time, n_sensors, n_sensors);
%                 for t = 1:n_time
%                     tmp_(t, :, :) = cov1para(Xpseudo_train(sigma_conditions==c, :, t));
%                 end
%                 sigma_(c, :, :) = mean(tmp_, 1);
%             end
%             sigma = squeeze(mean(sigma_, 1));  % average across conditions
%             sigma_inv = sigma^-0.5;
%             for t = 1:n_time
%                 Xpseudo_train(:, :, t) = squeeze(Xpseudo_train(:, :, t)) * sigma_inv;
%                 Xpseudo_test(:, :, t) = squeeze(Xpseudo_test(:, :, t)) * sigma_inv;
%             end

            for t = 1:n_time
                for c1 = 1:n_conditions-1
                    for c2 = c1+1:n_conditions
                        % 3. Fit the classifier using training data
                        data_train = Xpseudo_train(ind_pseudo_train(c1, c2, :), :, t);
                        y_train = squeeze(labels_pseudo_train(c1, c2, :));
                        model_svm = svmtrain(y_train, data_train, '-c 1 -q 0 -t 0'); % libSVM: -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"; -q : quiet mode (no outputs)\; "-t kernel_type : set type of kernel function (default 2)\n"                       
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
    result = result_; % result = 4D matrix; Beispiel: 2 *9*9*111 => 2 = Anzahl sessions, 9 = Anzahl an Bedingungen, 111 = Anzahl an Zeitpunkten
    result.time = sessions.time;
    result.volley = '1';
    result.space = '2';
    result.soccer = '3';
    
    if ~exist (outPath, 'dir')
        mkdir (outPath)
    end
    save([outPath 'result_decoding_' freq], 'result')
    
%% Now we plot the average classification accuracy time course by collapsing across sessions and conditions:
if 3 == length(conditions) 
    figure; axis tight
    hold on
    plot(sessions.time{1,1}, 100*squeeze(result.svm(1,2,:))+50, 'b-', 'linewidth', 1, 'DisplayName','Volley vs Space'); %legend('like vs dislike');
    plot(sessions.time{1,1}, 100*squeeze(result.svm(1,3,:))+50, 'r-', 'linewidth', 0.5, 'DisplayName', 'Volley vs Soccer');
    plot(sessions.time{1,1}, 100*squeeze(result.svm(2,3,:))+50, 'k-', 'linewidth',0.5, 'DisplayName', 'Space vs Soccer');
    legend('show')
    plot([-0.5 1], [50 50], 'k-')
    xlim([-0.5 1])
    xlabel('Time [s]');
    ylabel('Classification accuracy svm');
    savefig([outPath 'decoding_result_' freq '_svm.fig'])
    close
    figure; axis tight
    hold on
    plot(sessions.time{1,1}, 100*squeeze(result.weird(1,2,:))+50, 'b-', 'linewidth', 1)
    plot(sessions.time{1,1}, 100*squeeze(result.weird(1,3,:))+50, 'r-', 'linewidth', 0.5)
    plot(sessions.time{1,1}, 100*squeeze(result.weird(2,3,:))+50, 'k-', 'linewidth', 0.5)
    plot([-0.5 1], [50 50], 'k-')
    xlim([-0.5 1])
    xlabel('Time [s]')
    ylabel('Classification accuracy weird')
    legend('Volley vs Space', 'Volley vs Soccer', 'Soccer vs Soccer')
    savefig([outPath 'decoding_result_' freq '_weird.fig'])
    close
    figure; axis tight
    hold on
    plot(sessions.time{1,1}, 100*squeeze(result.gnb(1,2,:))+50, 'b-', 'linewidth', 1)
    plot(sessions.time{1,1}, 100*squeeze(result.gnb(1,3,:))+50, 'r-', 'linewidth', 0.5)
    plot(sessions.time{1,1}, 100*squeeze(result.gnb(2,3,:))+50, 'k-', 'linewidth', 0.5)
    plot([-0.5 1], [50 50], 'k-')
    xlim([-0.5 1])
    xlabel('Time [s]')
    ylabel('Classification accuracy gnb')
      legend('Volley vs Space', 'Volley vs Soccer', 'Soccer vs Soccer')
    savefig([outPath 'decoding_result_' freq '_gnb.fig'])
    close
else
    figure; axis tight
    hold on
    plot(sessions.time{1,1}, 100*squeeze(nanmean(nanmean(result.svm, 1), 2))+50, 'linewidth', 1) % nanmean = mean ignoring nans
    plot(sessions.time{1,1}, 100*squeeze(nanmean(nanmean(result.weird, 1), 2))+50, 'linewidth', 1)
    plot(sessions.time{1,1}, 100*squeeze(nanmean(nanmean(result.gnb, 1), 2))+50, 'linewidth', 1)
    plot([-0.5 1], [50 50], 'k-')
    xlim([-0.5 1])
    xlabel('Time [s]')
    ylabel('Classification accuracy')
    legend('SVM', 'WeiRD', 'GNB')
    savefig([outPath 'MEG\sourcespace\noROIs\Guggenmos_decoding_results\decoding_result_' freq 'SVM_WeiRD_BNB.fig'])
end

% Already for one participant and a reduced data set (10 insteada of 92 conditions), these results look 
% like canonical decoding time courses. Note that the period -100ms to 0ms is the baseline phase and stimulus onset is at 0ms. In this example, Support Vector Machine and WeiRD outperform Gaussian Naive Bayes.

 end
 %%
