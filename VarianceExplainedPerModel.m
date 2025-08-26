function [VarExpM, models, BestModel] = VarianceExplainedPerModel(X,scl,mult)
%% Variance explained by the factor
    SStot = sum(X(:).^2);

    models = [];
    for i=1:length(scl)-mult
        M = nchoosek(1:length(scl)-mult,i);
        models = cat(1,models,mat2cell(M, ones(size(M,1),1), size(M,2)));
    end

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

        model_scaled = gemanova(X,F,Fix,scl,cross,show);

        VarExpM(studied_model) = (SStot - sum(model_scaled.res(:).^2))/SStot;

        waitbar(studied_model/m, h);
    end
    close(h)

    [~,n] = max(VarExpM);
    BestModel = models{n};
end