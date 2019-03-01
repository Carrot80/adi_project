 %%
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
n_perm = 5;
ncond=length(unique(sessions.labels));
labels_cond = unique(sessions.labels);
ntrials_cond = histc(sessions.labels, unique(sessions.labels));
smaller_cond = find(ntrials_cond==min(ntrials_cond));
ntrials_diff = ntrials_cond(1)- ntrials_cond(2);
% n_rep = floor(ntrials_diff/5);

% ntrl = 1:5:ntrials_diff;
% n=1;
% for i=ntrl
%    data_eq_cond(n,:,:) = mean(sessions.data(ntrl:ntrl+5,:,:));
%    n=n+1; 
%             
% end            
% 
% sessions.data_eq_cond = cat(1,data_eq_cond, sessions.data(5*size(data_eq_cond,1)+1:size(sessions.data,1),:,:));
% sessions.labels_=sessions.labels(5*size(data_eq_cond,1)+1-size(data_eq_cond,1):length(sessions.labels));
% sessions.labels = sessions.labels_;


    % We define three classifiers that will be compared: Support Vector Machine, Weighted Robust Distance,
    % Gaussian Naive Bayes
     
    clfs = {'svm', 'gnb', 'weird'};
    for c = 1:length(clfs)
        result.(clfs{c}) = nan(n_sessions, n_perm, n_time);
    end
    
  for s = 1:n_sessions
        
        fprintf('sessions %g / %g\n', s, n_sessions)  
        
%         for f = 1:n_perm
%             fprintf('\tPermutation %g / %g\n', f, n_perm)
             % precompute permutations
%             n1(pp).test = [];
%             n1.test = cell(n_perm,1);  
%             rand_cond1 = randperm(ntrials_cond(1));
            rand_cond2 = randperm(ntrials_cond(2))+ ntrials_cond(1); 
            
            n_test = round(ntrials_cond(smaller_cond)./n_perm);
            n_train = ntrials_cond(smaller_cond) - n_test;
            
            for i=1:ncond
                
                cond_perm(1).(['cond' num2str(i)]).trls_excl_test = 1:ntrials_cond(i);
                cond_perm(1).(['cond' num2str(i)]).test = randperm(ntrials_cond(i),n_test);
               
                for op = 1:length(cond_perm(1).(['cond' num2str(i)]).test)
                    ind_ = find(1:ntrials_cond(i) == cond_perm(1).(['cond' num2str(i)]).test(op)); 
                    ind(op)=ind_;
                    clear ind_
                end
                cond_perm(1).(['cond' num2str(i)]).trls_excl_test(ind) = [];
                clear ind
                rand_ = randperm(length(cond_perm(1).(['cond' num2str(i)]).trls_excl_test), n_train);
                cond_perm(1).(['cond' num2str(i)]).train = cond_perm(1).(['cond' num2str(i)]).trls_excl_test(rand_);

                 for f=2:n_perm 
                     cond_perm(f).(['cond' num2str(i)]).trls_excl_test = cond_perm(f-1).(['cond' num2str(i)]).trls_excl_test;   
                     rand_t = randperm(length(cond_perm(f).(['cond' num2str(i)]).trls_excl_test), n_test);
                     cond_perm(f).(['cond' num2str(i)]).test = cond_perm(f).(['cond' num2str(i)]).trls_excl_test(rand_t);
                     clear rand_t

                     for op = 1:length(cond_perm(f).(['cond' num2str(i)]).test)
                        ind_ = find(cond_perm(f).(['cond' num2str(i)]).trls_excl_test == cond_perm(f).(['cond' num2str(i)]).test(op)); 
                        ind(op)=ind_;
                        clear ind_
                     end
                     cond_perm(f).(['cond' num2str(i)]).trls_excl_test(ind) = [];
                     clear ind
                     %train_indices:
                     cond_perm(f).(['cond' num2str(i)]).trls_excl_test2 = 1:ntrials_cond(i);
                     for op = 1:length(cond_perm(f).(['cond' num2str(i)]).test)
                        ind_ = find(cond_perm(f).(['cond' num2str(i)]).trls_excl_test2 == cond_perm(f).(['cond' num2str(i)]).test(op)); 
                        ind(op)=ind_;
                        clear ind_
                     end
                     cond_perm(f).(['cond' num2str(i)]).trls_excl_test2(ind)=[];
                     rand_ = randperm(length(cond_perm(f).(['cond' num2str(i)]).trls_excl_test2), n_train);
                     cond_perm(f).(['cond' num2str(i)]).train = cond_perm(f).(['cond' num2str(i)]).trls_excl_test2(rand_);
                     clear rand_
                 end
            end
            
        for f=1:n_perm 
            data_train = [sessions.data(cond_perm(f).(['cond' num2str(1)]).train, :,:); sessions.data(cond_perm(f).(['cond' num2str(2)]).train+ntrials_cond(1),:,:)]; 
            data_test = [sessions.data(cond_perm(f).(['cond' num2str(1)]).test, :,:); sessions.data(cond_perm(f).(['cond' num2str(2)]).test+ntrials_cond(1),:,:)]; 

            labels_train = [sessions.labels(cond_perm(f).(['cond' num2str(1)]).train) sessions.labels(cond_perm(f).(['cond' num2str(2)]).train+ntrials_cond(1))];
            labels_test = [sessions.labels(cond_perm(f).(['cond' num2str(1)]).test) sessions.labels(cond_perm(f).(['cond' num2str(2)]).test+ntrials_cond(1))];
                     

            %% 2. Whitening using the Epoch method
            sigma_conditions = labels_train;
            sigma_ = nan(n_conditions, n_sensors, n_sensors);
            for c = 1:n_conditions
                % compute sigma for each time point, then average across time
                tmp_ = nan(n_time, n_sensors, n_sensors);
                for t = 1:n_time
                    tmp_(t, :, :) = cov1para(data_train(sigma_conditions==c, :, t));
                end
                sigma_(c, :, :) = mean(tmp_, 1);
            end
            sigma = squeeze(mean(sigma_, 1));  % average across conditions
            sigma_inv = sigma^-0.5;
            clearvars sigma
            for t = 1:n_time
                data_train_sigma(:, :, t) = squeeze(data_train(:, :, t)) * sigma_inv;
                data_test_sigma(:, :, t) = squeeze(data_test(:, :, t)) * sigma_inv;
            end 
            clearvars sigma_inv data_train data_test
            
            for t = 1:n_time
                        % 3. Fit the classifier using training data
                        model_svm = svmtrain(labels_train', data_train_sigma(:, :, t), '-c 1 -q 0 -t 0'); % libSVM: -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"; -q : quiet mode (no outputs)\; "-t kernel_type : set type of kernel function (default 2)\n"                       
                        model_weird = weirdtrain(labels_train', data_train_sigma(:, :, t));
                        model_gnb = gnbtrain(labels_train', data_train_sigma(:, :, t));

                        % 4. Compute and store classification accuracies
                        result.svm(s, f, t) = ...
                            mean(svmpredict(labels_test',data_test_sigma(:, :, t), model_svm,'-q 0 -t 0')==labels_test');
                        result.weird(s, f, t) = ...
                            mean(weirdpredict(labels_test',data_test_sigma(:, :, t), model_weird)==labels_test');
                        result.gnb(s, f, t) = ...
                            mean(gnbpredict(labels_test',data_test_sigma(:, :, t), model_gnb)==labels_test');
                    clearvars model_svm model_weird model_gnb
            end
            clearvars data_train_sigma data_test_sigma data_train data_test 
        end
  end
  
    % average across permutations
    for c = 1:length(clfs)
        result_.(clfs{c}) = nan(n_sessions, n_perm, n_time);
    end 
    result_.svm = squeeze(nanmean(result.svm, 2));
    result_.gnb = squeeze(nanmean(result.gnb, 2));
    result_.weird = squeeze(nanmean(result.weird, 2));
    result_final = result_; % result = 4D matrix; Beispiel: 2 *9*9*111 => 2 = Anzahl sessions, 9 = Anzahl an Bedingungen, 111 = Anzahl an Zeitpunkten
    result_final.time = sessions.time;
    result_final.like = '1';
    result_final.dislike = '2';
    save([outPath 'result_decoding_' freq], 'result_final')
    
%% Now we plot the average classification accuracy time course by collapsing across sessions and conditions:

figure; axis tight
hold on
plot(result_final.time{1,1}, 100*result_.svm, 'linewidth', 1)
plot(result_final.time{1,1}, 100*result_.weird, 'linewidth', 1)
plot(result_final.time{1,1}, 100*result_.gnb, 'linewidth', 1)
plot([-0.5 1], [50 50], 'k-')
xlim([-0.5 1])
xlabel('Time [s]')
ylabel('Classification accuracy')
title('svm')
legend('SVM', 'WeiRD', 'GNB')
savefig([outPath 'decoding_result_' freq 'SVM_WeiRD_BNB.fig'])


% Already for one participant and a reduced data set (10 insteada of 92 conditions), these results look 
% like canonical decoding time courses. Note that the period -100ms to 0ms is the baseline phase and stimulus onset is at 0ms. In this example, Support Vector Machine and WeiRD outperform Gaussian Naive Bayes.




 %%
end