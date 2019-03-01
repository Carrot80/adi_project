 %% ohne noise-whitenting
function [] = SVM_Guggenmos_tutorial(sessions, outPath, freq)

avg=squeeze(mean(sessions.data));
sRate=256;    
[FourRef,Fref]=fftBasic(avg,round(sRate));
figure
plot(FourRef, Fref)
if ~exist (outPath, 'dir')
    mkdir (outPath)
end
% savefig ([outPath filesep 'freqspectrum_' freq '.fig'])
close

% We set a seed, in order to make analyses reproducible:
rng('default');
% rng(10);


n_conditions = length(unique(sessions(1).labels));
if 1==isequal(n_conditions, 3) && 0==all(sessions(1).labels)
    ind=find(sessions(1).labels==0);
    sessions(1).labels(ind)=[];
    sessions(1).data(ind,:,:)=[];
end

% labels_rand = sessions.labels(randperm(length(sessions.labels)));
% sessions.labels=labels_rand;

n_conditions = length(unique(sessions(1).labels));
n_sensors = size(sessions(1).data, 2);
n_time = size(sessions(1).data, 3);
n_sessions = length(sessions);
n_perm = 10;
ncond=length(unique(sessions.labels));
labels_cond = unique(sessions.labels);
ntrials_cond = histc(sessions.labels, unique(sessions.labels));

% The analytic logic is contained in a nested for loop, with loops for the number of sessions, 
% number of permutations, number of timepoints, number of conditions, and number of conditions again. 
% Overall, the logic contains 4 crucial steps:

% 2. Whiten the training data (here using the Epoch method, which is recommended in our manuscript)
% 3. Fit the classifier to the training data
% 4. Compute classification accuracy on test data
% 

    % We define three classifiers that will be compared: Support Vector Machine, Weighted Robust Distance,
    % Gaussian Naive Bayes
    clfs = {'svm', 'gnb', 'weird'};
    for c = 1:length(clfs)
        result.(clfs{c}) = nan(n_sessions, n_perm, n_time);
    end
    
    for s = 1:n_sessions
        
        fprintf('sessions %g / %g\n', s, n_sessions)  
        
        for f = 1:n_perm
            fprintf('\tPermutation %g / %g\n', f, n_perm)
             % precompute permutations
            
             
             
            rand_cond1 = randperm(ntrials_cond(1));
            rand_cond2 = randperm(ntrials_cond(2))+ ntrials_cond(1); 

            trials_train = [rand_cond1(1:round(ntrials_cond(1)*0.8)) rand_cond2(1:round(ntrials_cond(2)*0.8))];
            trials_test = [rand_cond1(round(ntrials_cond(1)*0.8)+1:ntrials_cond(1)) rand_cond2(round(ntrials_cond(2)*0.8)+1:ntrials_cond(2))];

            data_train = sessions.data(trials_train, :,:);
            data_test = sessions.data(trials_test, :,:);

            labels_train = sessions.labels(trials_train);
            labels_test = sessions.labels(trials_test);
                     
            for t = 1:n_time
                        % 3. Fit the classifier using training data
                        model_svm = svmtrain(labels_train', data_train(:, :, t), '-c 1 -q 0 -t 0'); % libSVM: -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"; -q : quiet mode (no outputs)\; "-t kernel_type : set type of kernel function (default 2)\n"                       
                        model_weird = weirdtrain(labels_train', data_train(:, :, t));
                        model_gnb = gnbtrain(labels_train', data_train(:, :, t));

                        % 4. Compute and store classification accuracies
                        result.svm(s, f, t) = ...
                            mean(svmpredict(labels_test',data_test(:, :, t), model_svm,'-q 0 -t 0')==labels_test');
                        result.weird(s, f, t) = ...
                            mean(weirdpredict(labels_test',data_test(:, :, t), model_weird)==labels_test');
                        result.gnb(s, f, t) = ...
                            mean(gnbpredict(labels_test',data_test(:, :, t), model_gnb)==labels_test');
                    clearvars model_svm model_weird model_gnb
            end
            clearvars data_train data_test rand_cond1 rand_cond2
        end
    end
    % average across permutations
    for c = 1:length(clfs)
        result_.(clfs{c}) = nan(n_sessions, n_perm, n_time);
    end 
    result_.svm = squeeze(nanmean(result.svm, 2));
    result_.gnb = squeeze(nanmean(result.gnb, 2));
    result_.weird = squeeze(nanmean(result.weird, 2));
%     result_ = result_; % result = 4D matrix; Beispiel: 2 *9*9*111 => 2 = Anzahl sessions, 9 = Anzahl an Bedingungen, 111 = Anzahl an Zeitpunkten
    result_.time = sessions.time;
    result_.like = '1';
    result_.dislike = '2';
    save([outPath 'result_decoding_' freq], 'result_')
    
%% Now we plot the average classification accuracy time course by collapsing across sessions and conditions:

figure; axis tight
hold on
plot(result_.time{1,1}, 100*result_.svm, 'linewidth', 1)
plot(result_.time{1,1}, 100*result_.weird, 'linewidth', 1)
plot(result_.time{1,1}, 100*result_.gnb, 'linewidth', 1)
plot([-0.5 1], [50 50], 'k-')
xlim([-0.5 1])
ylim([0 100])
xlabel('Time [s]')
ylabel('Classification accuracy')
title('svm')
legend('SVM', 'WeiRD', 'GNB')
% legend('SVM', 'WeiRD')
savefig([outPath 'decoding_result_' freq 'SVM_WeiRD_BNB.fig'])


% Already for one participant and a reduced data set (10 insteada of 92 conditions), these results look 
% like canonical decoding time courses. Note that the period -100ms to 0ms is the baseline phase and stimulus onset is at 0ms. In this example, Support Vector Machine and WeiRD outperform Gaussian Naive Bayes.

 end
 %%
