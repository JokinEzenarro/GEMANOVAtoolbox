%% Input data
X = NIRbanco.data(:,1:10:end); % Response matrix
F = DoE_NIRbanco.data(:,[1 4:6]); % DoE matrix
F(:,2:4) = F(:,2:4)+1;
X = X(F==6,:);
F = F(F==6,2:end);
% X needs to be a data matrix (double) with samples in rows and variables
% in columns
% F needs to be a DoE matrix (double) with samples in rows and factors in
% columns. Levels coded from 1 to n, continuously.

%% Check for replicates and average if there are
[X, F] = average_replicates(X, F);

%% Preprocessing of signal
X=X-mean(X,1,'omitnan'); % Mean center

%% Convert matrix to hypercube
[hypercube, idx_cell, scl, multivariate] = MatrixToHypercube(X,F);

%% Check best model (Optional)
[VarExpM, models, BestModel] = VarianceExplainedPerModel(hypercube,scl,multivariate);
ModelNames = cellfun(@mat2str, models, 'UniformOutput', false);

figure
bar(VarExpM*100) 
xticks(1:numel(VarExpM))
xticklabels(ModelNames)
ylim([0 100])
ylabel('Explained variance (%)')

%% Fit GEMANOVA model
IncludedFactors = [1:length(scl)]; %Full-way model
[model,loadings,VarExp] = FitModel(hypercube,scl,IncludedFactors);
disp("Explained variance: "+round(VarExp*100,2)+"%"+newline)

%% Plot results
figure
for i=1:length(scl)-multivariate
    nexttile
    plot(scl{i},loadings{i})
end
if multivariate
    nexttile
    plot(loadings{end})
end

s = size(model.res);
Xrec=1;
for i=1:size(loadings,2)
    Xrec=tensorprod(Xrec,loadings{i});
end
Xrec=squeeze(squeeze(Xrec));
if multivariate == 0
    nexttile
    plot(X,reshape(Xrec,[prod(s(1:end)) 1]),'o')
end

%% Permutation test for full factors
p = 1000; %Number of permutations
plots = 1; %Show results
pValues = FullFactorPermutation(hypercube, scl, p, model, plots, multivariate);

%% Pair-wise permutation test for levels in factors
p = 1000; %Number of permutations
plots = 1; %Show results
for f=1:length(scl)
    out = PairwiseFactorPermutation(hypercube, scl, p, model, f, plots);
end