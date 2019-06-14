function [] = adi_glm (data_table)

%% simple logistic regression:

y = data_table.response;
X = table2array(data_table(:,3:end-1));

ind_like = find(y==1);
ind_dislike = find(y==2);
y(ind_dislike) = 0; % Werte auf 0 und 1 setzen
ind_dislike = find(y==0);

figure
plot(X(ind_like,:), 'kx', 'MarkerSize', 5 )
hold on
plot(X(ind_dislike,:), 'ko', 'MarkerSize', 5, 'Color', 'r' )
legend('like', 'dislike')
title('like_vs_dislike')


% compute cost (least square error) and gradient:
[m, n] = size(X);
initialTheta = zeros((n+1),1);
X = [ ones(m,1), X];
[J, grad]=computeCost(initialTheta, X, y); 
% theta = fminunc(@(t)computeCost(t, X, y), initialTheta, options);
% run the function optimization algorithm:
lambda = 0.0;  % lambda prevents overfitting, change  (zwischen 0 und 1?)
options = optimset('GradObj', 'on', 'MaxIter', 8000); %  set the GradObj option to on, which tells fminunc that our function returns both the cost and the gradient.

% Run fminunc to obtain the optimal theta, % This function will return theta and the cost
% the learning rate alpha does not need to be set, this does fminunc automatically 
% [optTheta, cost, exitFlag] = fminunc(@(t)(costFunction(t, X, y)), initialTheta, options);
[optTheta, cost, exitFlag] = fminunc(@(t)(costFunctionReg(t, X, y, lambda)), initialTheta, options);

% help fminunc % fminunc means function minimization unconstrained
% check the accuracy of predictions:
predictions = adi_predict(optTheta,X);
accuracy = mean(double(predictions == y) * 100);

[sort_theta, ind_sort] = sort(abs(optTheta), 'descend')
X(ind_sort)
data_table.Properties.VariableNames(2:end-1)
cellstr(data_table.Properties.VariableNames(2:end-1))


for k=1:length(data_table.Properties.VariableNames)
    features(k)=string(data_table.Properties.VariableNames(k));
end
sort_features = features(ind_sort);

for k = 1:length(data_table.response)

    switch data_table.response(k) 
        case 1

            y_cat(k)='like';

        case 2
            y_cat(k)='dislike';
    end
end



[B_cat,dev_cat,stats_cat] = mnrfit(X,y_cat) % y darf nicht aus 0 und 1 bestehen, sondern aus 1 und 2;


%% sortiere einzelne Komponenten:

Table_sorted = sortrows(data_table(2:end-1),column)

log_regression.all_comp.feature_names = data_table.Properties.VariableNames(2:end-1);
log_regression.all_comp.p = stats.p;
log_regression.all_comp.B = B;

log_regression.comp1.feature_names = data_table.Properties.VariableNames([2 3:3:end-1]);
log_regression.comp1.p = stats.p([1 2:3:end]);
log_regression.comp1.B = B([1 2:3:end]);

log_regression.comp2.feature_names = data_table.Properties.VariableNames([2 4:3:end-1]);
log_regression.comp2.p = stats.p([1 3:3:end]);
log_regression.comp2.B = B([1 3:3:end]);

log_regression.comp3.feature_names = data_table.Properties.VariableNames([2 5:3:end-1]);
log_regression.comp3.p = stats.p([1 4:3:end]);
log_regression.comp3.B = B([1 4:3:end]);

[log_regression.comp1.sort_p, ind_sort_p] = sort(log_regression.comp1.p, 'ascend');
log_regression.comp1.sort_features = transpose(log_regression.comp1.feature_names(ind_sort_p));
log_regression.comp1.sort_B = log_regression.comp1.B(ind_sort_p);

[log_regression.comp2.sort_p, ind_sort_p] = sort(log_regression.comp2.p, 'ascend');
log_regression.comp2.sort_features = transpose(log_regression.comp2.feature_names(ind_sort_p));
log_regression.comp2.sort_B = log_regression.comp2.B(ind_sort_p);

[log_regression.comp3.sort_p, ind_sort_p] = sort(log_regression.comp3.p, 'ascend');
log_regression.comp3.sort_features = transpose(log_regression.comp3.feature_names(ind_sort_p));
log_regression.comp3.sort_B = log_regression.comp3.B(ind_sort_p);


%% log. Regression getrennt für einzelne Komponenten:

[B_comp1,dev_comp1,stats_comp1] = mnrfit(X(:,1:3:end),y_cat);
[B_comp2,dev_comp2,stats_comp2] = mnrfit(X(:,2:3:end),y_cat);
[B_comp3,dev_comp3,stats_comp3] = mnrfit(X(:,3:3:end),y_cat);

[comp1.sort_p ind_sort_p] = sort(stats_comp1.p(:,3), 'ascend');
comp1.features_sorted = transpose(log_regression.comp1.feature_names(ind_sort_p))
comp1.B_sorted = B_comp1(ind_sort_p);

[comp2.sort_p ind_sort_p] = sort(stats_comp2.p(:,3), 'ascend');
comp2.features_sorted = transpose(log_regression.comp2.feature_names(ind_sort_p));
comp2.B_sorted = B_comp2(ind_sort_p);

[comp3.sort_p ind_sort_p] = sort(stats_comp3.p(:,3), 'ascend');
comp3.features_sorted = transpose(log_regression.comp3.feature_names(ind_sort_p));
comp3.B_sorted = B_comp3(ind_sort_p,3)';


%%




plotDecisionBoundary2(optTheta (1:3), X(:,2:3), y) % funktioniert nicht
%% non-linear logistic regression mixed model using SimBiology:

[beta,PSI,stats,B] = nlmefit(X,y,group,V,fun,beta0)


%%
% fixed effects: balldesign(?), ROI
% randomn effect: Subject, time, 

% hwo to deal with unbalanced data sets:
unstacked_panel = unstack(panel, 'log_GDP', 'State', 'Grouping_variable', 'YR')
unstacked_panel = sortrows(unstacked, 'YR', 'ascend');

% start with a simple model and improve it based on its results:


%1. fitlm
lm = fitlm(data_table, 'response ~ 1 + Temporal_Mid_L_comp1') % geht nicht,da response dichotom ist
% 2. linear mixed effect model:
lme = fitlme('log_GDP ~ 1 + Ycentered + (1+YRcentered | State)') % ==> (1+YRcentered | State): means we specified random effects here: we introduced a randomn effect for the intercept (1) as well as the slope (YRcentered) with state as grouping variable
lme = fitlme(data_table, 'response ~ 1 + Temporal_Mid_L_comp1 + (1+Temporal_Mid_L_comp1 | subject)') 

%2. fitglm
mdl =  fitglm(X,y,'Distribution','binomial') % nicht signifikant
mdl =  fitglm(data_table, 'Distribution','binomial') % => nicht signifikant


glme = fitglme(data_table,'response ~ 1 + (RANDOM_1 | Grp_1) + ... + (RANDOM_R | Grp_R)',Name,Value)
 % => mixed model with repeated measurement

% y ~ FIXED + (RANDOM_1 | Grp_1) + ... + (RANDOM_R | Grp_R)'


%%
compare(lm, lme)
[~,~, RE] = randomEffects(lme)
[~,~,stats]=covarianceParameters(lme)
disp(stats{1})
fitlmematrix


glme = fitglme(tbl,'response ~ 1 + (RANDOM_1 | Grp_1) + ... + (RANDOM_R | Grp_R)',Name,Value) % => mixed model with repeated measurement
glme = fitglme(tbl,'response ~ 1 + (RANDOM_1 | Grp_1) + ... + (RANDOM_R | Grp_R)',Name,Value) % => mixed model with repeated measurement


glme = fitglme(mfr,'defects ~ 1 + newprocess + time_dev + temp_dev + supplier + (1|factory)', 'Distribution','binomial','Link','log','FitMethod','Laplace', ...
    'DummyVarCoding','effects');

%%


options = optimset('GradObj', 'on', 'MaxIter', 100); % 100 iterations
initialTheta= zeros(2,1);
[optTheta, functionVal, exitFlag] = fminunc(@costFunction, initialTheta, options)
help fminunc % fminunc means function minimization unconstrained



















end



function [] = adi_glm (data_table)

% Create New Figure
figure; hold on;

% ====================== YOUR CODE HERE ======================
% Instructions: Plot the positive and negative examples on a
%               2D plot, using the option 'k+' for the positive
%               examples and 'ko' for the negative examples.
%


% Find Indices of Positive and Negative Examples
pos = find(y==1); neg = find(y == 0);
% Plot Examples
figure;
plot(X(pos, 1), X(pos, 2), 'k+','LineWidth', 2, ...
'MarkerSize', 7);
hold on
plot(X(neg, 1), X(neg, 2), 'ko', 'MarkerFaceColor', 'y', ...
'MarkerSize', 7);
legend('admitted', 'not admitted')
xlabel('exam1 score');
ylabel('exam2 score');

end


function [jVal, gradient] = costFunction(theta)


jVal = (theta(1)-5)^2 + (theta(2)-5)^2; % log einfügen

gradient = zeros(2,1);
gradient(1) = 2*(theta(1)-5);
gradient(2) = 2*(theta(2)-5); % ab theta 2 lambda einfügen
% gradient(n).. 

end


