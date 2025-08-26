function [model,loadings,VarExp] = FitModel(X,scl,IncludedFactors)
%% Fit GEMANOVA model
% X = response hypercube
% f = number of factors

    F=1; % Number of components
    
    Fix(:,1) = ones(length(scl),1); % factors to consider (None)
    Fix(IncludedFactors,1) = 0;
    
    cross=0;
    show=0;
    
    model = gemanova(X,F,Fix,scl,cross,show);
    loadings = model.effects;
    
    SStot = sum(X(:).^2);
    
    SSmodel = SStot - sum(model.res(:).^2);
    VarExp = SSmodel/SStot;
