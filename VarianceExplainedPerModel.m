function [VarExpM, models] = VarianceExplainedPerModel(X,scl,mult)
% VarianceExplainedPerModel:
% Computes the variance explained for all possible combinations of factors 
% (excluding the multivariate mode if present) using GEMANOVA models.
%
% INPUTS:
%   X     - N-dimensional response hypercube (multi-way array)
%   scl   - cell array with axis levels for each mode
%   mult  - Is it a multivariate signal? (1 = Yes, 0 = No)
%
% OUTPUTS:
%   VarExpM   - vector with variance explained for each tested model
%   models    - cell array containing the factor index sets for each model

    % Total sum of squares
    SStot = sum(X(:).^2);

    % Generate all combinations of factors (excluding multivariate mode)
    models = [];
    for i=1:length(scl)-mult
        M = nchoosek(1:length(scl)-mult,i);
        models = cat(1,models,mat2cell(M, ones(size(M,1),1), size(M,2)));
    end

    % Loop over all candidate models
    h = waitbar(0, 'Trying models...');
    m = length(models);
    for studied_model=1:m
        F=1; % Number of components
        Fix(:,1) = ones(length(scl) ,1); % factors to consider (None)
        if mult == 1
            Fix(end,1) = 0;
        end
        Fix(models{studied_model}) = 0;
        cross=0;
        show=0;

        % Fit GEMANOVA for this combination of factors
        model_scaled = gemanova(X,F,Fix,scl,cross,show);

        % Compute explained variance for this model
        VarExpM(studied_model) = (SStot - sum(model_scaled.res(:).^2))/SStot;

        waitbar(studied_model/m, h);
    end
    close(h)
end