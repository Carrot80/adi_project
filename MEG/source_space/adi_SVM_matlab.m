



function SVM_matlab(path2data, outPath_extdisc, freqname, condition)

[session] = mk_SVM_struct([], path2data, outPath_extdisc, freqname, condition);
num_conditions = length(unique(session.labels));
% ionosphere = array2table(X);
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

%%



%%
load fisheriris
species_num = grp2idx(species);
X = randn(100,10);
X(:, [1, 3, 5, 7]) = meas(1:100,:); 
y = species_num(1:100,:);

% 80:20
rand_num = randperm(100);
X_train = X(rand_num(1:80),:);
y_train = y(rand_num(1:80));

X_test = X(rand_num(81:end),:);
y_test = y(rand_num(81:end));

%% CV partition

c = cvpartition(y_train, 'k', 5);

%% feature selection
opts = statset('display', 'iter');
fun = @(train_data, train_labels, test_data,  test_labels)...
sum(predict(fitcsvm(train_data, train_labels, 'KernelFunction', 'rbf'), test_data) ~= test_labels)
[fs, history] = sequentialfs(fun, X_train, y_train, 'cv', c, 'options', opts, 'nfeatures', 2);
X_train_w_best_features = X_train(:,fs);   

Md1 = fitcsvm(X_train_w_best_features, y_train, 'KernelFunction', 'rbf', 'OptimizeHyperparameters', 'auto', ...
    'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName', 'expected-improvement-plus', 'ShowPlots', true));


%% test





model_matlab_svmtoolbox = fitcsvm(data_train, y_train)
SVMModel = fitcsvm(X,y);

%                         cv = cvpartition(200,'KFold',10);
                        opts = struct('Optimizer','bayesopt','ShowPlots',true,'CVPartition',cv,...
                        'AcquisitionFunctionName','expected-improvement-plus');
                        model_matlab_svmtoolbox = fitcsvm(data_train, y_train,'KernelFunction','rbf',...
                        'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts)
                        model_matlab_svmtoolbox = fitcsvm(data_train, y_train, 'CacheSize', 'maximal', 'CrossVal', 'on', 'KernelFunction', 'linear'); % Matlab machine learning toolbox: Default: 'linear' for two-class
%                       learning and 'gaussian' (or 'rbf') for one-class learning
                           result.svm_matlab_svmtoolbox(s, f, c1, c2, t) = ...
                            mean(svmpredict(y_train,data_test,model_matlab_svmtoolbox)==y_train)-0.5;   


end


 function [session] = mk_SVM_struct(session, path2data, outPath_extdisc, freqname, condition)
 
 fileName = [outPath_extdisc 'MEG\sourcespace\noROIs\runs_appended\virtsens\virtsens_all_conditions_allRuns_', freqname, '.mat'];
 if ~exist(fileName, 'file')
    [vs_allRuns] = adi_appendvirtsensMEG(path2data, outPath_extdisc, freqname, condition);
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
    session.labels = cat(1, ind_field_1*ones(1, length(vs_allRuns.(fieldnames{1}).trial)), ind_field_2*ones(1, length(vs_allRuns.(fieldnames{2}).trial)));
end
for i = 1:length(data)
    temp = data{1,i};
    session.data(i,:,:) = temp;
end
session.time = vs_allRuns.(fieldnames{1}).time;

end
  
