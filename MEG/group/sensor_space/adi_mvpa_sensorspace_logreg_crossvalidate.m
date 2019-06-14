function stats(like, dislike, cv_like, cv_dislike, time)

if isempty(time)
    fprintf('no time input, computing classification across time...')

end
%%% Train and test a classifier "by hand", i.e. without the
%%% crossvalidation and classification across time functions provided by
%%% MVPA-Light


%% Let's have a look at the data first: Calculate and plot ERP for attended and unattended deviants

% ERP for each condition
cfg = [];
cfg.keeptrials = 'yes';
avg_like_training = ft_timelockanalysis(cfg, like);
avg_dislike_training = ft_timelockanalysis(cfg, dislike);
avg_like_test = ft_timelockanalysis(cfg, cv_like);
avg_dislike_test = ft_timelockanalysis(cfg, cv_dislike);


% Plot ERP: 
close
figure
h1 = plot(avg_like_training.time, squeeze(mean(avg_like_training.trial)), 'r'); hold on
h2 = plot(avg_dislike_training.time, squeeze(mean(avg_dislike_training.trial)), 'b');
grid on
xlabel('Time [s]'),ylabel('MEG amplitude')
title('ERP')
legend([h1(1),h2(1)],{'like_training', 'dislike_training'})


% Plot testdata:

figure
h1 = plot(avg_like_test.time, squeeze(mean(avg_like_test.trial)), 'r'); hold on
h2 = plot(avg_dislike_test.time, squeeze(mean(avg_dislike_test.trial)), 'b');
grid on
xlabel('Time [s]'),ylabel('MEG amplitude')
title('ERP')
legend([h1(1),h2(1)],{'like_test', 'dislike_test'})

%% Train and test classifier

Xtrain = cat(1, avg_like_training.trial, avg_dislike_training.trial);
clabel_train = [ones(1, size(like.trial,2)) 2*ones(1, size(dislike.trial,2))];

Xtest = cat(1, avg_like_test.trial, avg_dislike_test.trial);
clabel_test = [ones(1, size(cv_like.trial,2)) 2*ones(1, size(cv_dislike.trial,2))];



% Get default hyperparameters for the LDA classifier
param = mv_get_classifier_param('logreg');
param.reg = 'l2';
param.lambda = [0.0001 0.0006 0.0036 0.026 0.13 0.1 0.2 0.3 0.4 0.5 0.6 0.77 1 4.6416 27.82559 166.81 1000];
param.plot = 0;
param.normalise = 'none'; % ==> macht er an der Stelle nicht 

% Train an logreg classifier
for tt = 1:size(Xtrain,3)
    
    cf(tt) = train_logreg(param, Xtrain(:,:,tt), clabel_train);

    % Test classifier on the same data: the function gives the predicted
    % labels (predlabel), the decision values (dval) which represent the
    % distance to the hyperplane, and the class probability for the sample
    % belonging to class 1 (prob)
    [predlabel, dval, prob] = test_logreg(cf(tt), Xtest(:,:,tt));
    % Calculate AUC
    auc(tt) = mv_calculate_performance('auc', 'dval', dval, clabel_test);
    acc(tt) = mv_calculate_performance('acc', 'dval', dval, clabel_test);
end

figure
plot(like.time{1}, acc)
figure
plot(like.time{1}, auc)


% To calculate classification accuracy, compare the predicted labels to
% the true labels and take the mean
fprintf('Classification accuracy: %2.2f\n', mean(predlabel==clabel))



% Look at the distribution of the decision values. dvals should be positive
% for clabel 1 (attended deviant) and negative for clabel 2 (unattended
% deviant). dval = 0 is the decision boundary
figure(2)
clf
boxplot(dval, clabel)
hold on
plot(xlim, [0 0],'k--')
ylabel('Decision values')
xlabel('Class')
title('Distribution of decision values')

% Look at the distribution of the probabilities. prob should be higher
% for clabel 1 (attended deviant) than clabel 2 (unattended
% deviant)
figure(3)
clf
boxplot(prob, clabel)
hold on
plot(xlim, [0.5 0.5],'k--')
ylabel('Probability for class 1')
xlabel('Class')
title('Distribution of class probabilities')

end