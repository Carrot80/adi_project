% function predictions = weirdpredict(y, X, model, verbose)
%
% This functions computes predictions for a dataset X given a model 
% computed by weirdtrain.
%
% INPUT:
%   y: label vector (n_samples x 1)
%   X: testing data (n_samples x n_features)
%   model: struct as returned by weirdtrain
%   verbose: if true, prints accuracy (default: true)
%
% OUTPUT:
%   predictions: predicted label for each sample in X
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
% by Matthias Guggenmos 06/07/16
% _________________________________________________________________________
% Matthias Guggenmos, Katharina Schmack and Philipp Sterzer, "WeiRD - a fast
% and performant multivariate pattern classifier", 2016, International Workshop
% on Pattern Recognition in Neuroimaging (PRNI), Trento, Italy, 2016, pp. 1-4.
% doi: 10.1109/PRNI.2016.7552349


function predictions = weirdpredict(y, X, model, robust, verbose)

    if nargin <= 3 || isempty(robust)
        robust = false;
    end

    if nargin <= 4 || isempty(verbose)
        verbose = false;
    end

    % compute votes of each feature based on distance to centroid
    votes = abs(X - repmat(model.x1, size(X, 1), 1)) - ...
            abs(X - repmat(model.x2, size(X, 1), 1));

    % compute decision values as a weighted sum of votes and feature
    % importances
    fi_matrix = repmat(model.feature_importances_, size(votes, 1), 1);
    dec = dot(votes, fi_matrix, 2) / size(votes, 2);
    
    % compute predictions based on the sign of the decision value
    predictions = arrayfun(@(x)model.classes(x+1), dec > 0);
    
    if verbose
        fprintf('Accuracy: %.1f%%\n', 100*mean(predictions == y))
    end

end