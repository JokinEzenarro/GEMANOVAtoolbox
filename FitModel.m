function [model, loadings, VarExp] = FitModel(X, scl, IncludedFactors, plot)
% FitModel:
% Fits a 1-component GEMANOVA model and returns the model, factor loadings,
% and total explained variance. Optionally produces basic diagnostic plots.
%
% INPUTS:
%   X               - N-dimensional response hypercube (e.g., N1 x N2 x ... x Nk)
%   scl             - cell array with axis levels for each mode (length = ndims(X))
%   IncludedFactors - vector of factor indices to be estimated (set free in Fix)
%   plot            - 1 to generate plots; 0 to skip
%
% OUTPUTS:
%   model    - GEMANOVA model struct returned by gemanova(...)
%   loadings - cell array of factor effects/loadings (model.effects)
%   VarExp   - total explained variance (SS_model / SS_total)

    % Define number of components (currently fixed at 1)
    F = 1;  
    
    % Create Fix vector: by default, all factors fixed (=1).
    % Set the selected IncludedFactors to 0 so they are estimated.
    Fix(:,1) = ones(length(scl),1); 
    Fix(IncludedFactors,1) = 0;
    
    cross = 0;
    show = 0;
    
    % Fit GEMANOVA model
    model = gemanova(X, F, Fix, scl, cross, show);
    loadings = model.effects;
    
    % Compute explained variance
    SStot   = sum(X(:).^2);                     % total sum of squares
    SSmodel = SStot - sum(model.res(:).^2);     % model sum of squares
    VarExp  = SSmodel / SStot;                  % explained variance ratio
    
    % Plot results (optional)
    if plot == 1
        figure
        
        % Plot factor loadings for each factor except the multivariate one
        for i = 1:length(scl)-multivariate   % NOTE: multivariate undefined here
            nexttile
            plot(scl{i}, loadings{i})
            xlabel(sprintf('Factor %d scale', i))
            ylabel('Loading')
            title(sprintf('Loadings for Factor %d', i))
        end
        
        % Plot loadings for multivariate factor (if present)
        if multivariate
            nexttile
            plot(loadings{end})
            title('Multivariate loadings')
        end
        
        % Reconstruct X from loadings (approximation)
        s = size(model.res);
        Xrec = 1;
        for i = 1:size(loadings,2)
            Xrec = tensorprod(Xrec, loadings{i}); % requires tensorprod
        end
        Xrec = squeeze(squeeze(Xrec));
        
        % Scatterplot of measured vs reconstructed data (univariate only)
        if multivariate == 0
            nexttile
            plot(X, reshape(Xrec,[prod(s(1:end)) 1]), 'o')
            xlabel('Observed')
            ylabel('Reconstructed')
            title('Observed vs. Reconstructed')
        end
    end
end
