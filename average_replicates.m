function [X_avg, F_unique] = average_replicates(X, F)
% Averages rows in X that correspond to identical rows in F (replicates)
% 
% Inputs:
%   X - (n_samples x n_features) matrix (e.g., spectra)
%   F - (n_samples x n_factors) matrix (e.g., DoE factors and interactions)
%
% Outputs:
%   X_avg - (n_unique x n_features) matrix with averaged replicates
%   F_unique - (n_unique x n_factors) matrix of unique factor combinations

    % Convert F to string rows for comparison
    [F_unique, ~, idx] = unique(F, 'rows', 'stable');
    n_unique = size(F_unique, 1);
    X_avg = nan(n_unique, size(X,2));

    for i = 1:n_unique
        X_avg(i,:) = mean(X(idx == i, :), 1, 'omitnan');
    end
end