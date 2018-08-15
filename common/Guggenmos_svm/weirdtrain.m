% function model = weirdtrain(y, X, exponential_scaling)
%
% This functions computes the 'model' of the WeiRD algorithm on a 
% training dataset X.
%
% INPUT:
%   y: label vector (n_samples x 1)
%   X: training data (n_samples x n_features)
%   exponential_scaling (optional): if true, ttest2-based feature 
%       importances are scaled exponentially (default: false)
%
% OUTPUT:
%   the returned 'model' of the weird algorithm is a struct with the 
%   following fields:
%         model.x1: centroid of class 1
%         model.x2: centroid of class 2
%         model.feature_importances_: corresponds to ttest2-values for each
%           feature, exponentially scaled if exponential_scaling==true
%         model.classes: unique entries of y (aka classes)
%
% WeiRD stands for "Weighted Robust Distance" and is a fast and simple
% classification algorithm that assigns class labels based on the distance
% to class prototypes. The distance is the euclidian distance between a
% current sample and a prototype in a space, in which each feature
% dimension is scaled by the two-sample t-value of the respective feature 
% in the training data. Class prototypes correspond to the averages of
% each feature in the training data. The current implementation works for
% two-class problems only.
%
% by Matthias Guggenmos 19/09/16
% _________________________________________________________________________
% M. Guggenmos, K. Schmack and P. Sterzer, "WeiRD - a fast and performant 
% multivoxel pattern classifier," 2016 International Workshop on Pattern 
% Recognition in Neuroimaging (PRNI), Trento, Italy, 2016, pp. 1-4.
% doi: 10.1109/PRNI.2016.7552349

function model = weirdtrain(y, X, exponential_scaling)

    if nargin == 2 || isempty(exponential_scaling)
        exponential_scaling = false;
    end

    model = struct();
    model.classes = unique(y);
    x1 = X(y==model.classes(1), :);
    x2 = X(y==model.classes(2), :);
    n1 = size(x1, 1);
    n2 = size(x2, 1);
    model.x1 = sum(x1, 1) ./ n1;
    model.x2 = sum(x2, 1) ./ n2;
    
    % two-sample t-test with unequal sample size
    var_x1 = sum(abs(bsxfun(@minus, x1, model.x1)).^2, 1) ./ (n1 - 1);
    var_x2 = sum(abs(bsxfun(@minus, x2, model.x2)).^2, 1) ./ (n1 - 1);
    gsd = (((n1 - 1) * var_x1 + (n2 - 1) * var_x2) / (n1 + n2 - 2)).^0.5;
    t = (model.x1 - model.x2) ./ (gsd * (1 / n1 + 1 / n2).^0.5);
    model.feature_importances_ = abs(t);
    
    if exponential_scaling
        model.feature_importances_ = exp(model.feature_importances_);
    end


end