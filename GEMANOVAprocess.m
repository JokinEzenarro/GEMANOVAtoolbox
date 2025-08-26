%% Input data
X = data; % Response matrix
F = DoE; % DoE matrix

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

%% Create factor contribution chart
[VarExpM, models] = VarianceExplainedPerModel(hypercube,scl,multivariate);
ModelNames = cellfun(@mat2str, models, 'UniformOutput', false);

figure
bar(VarExpM*100) 
xticks(1:numel(VarExpM))
xticklabels(ModelNames)
ylim([0 100])
ylabel('Explained variance (%)')

%% Fit GEMANOVA model
plots = 1;
IncludedFactors = [1:length(scl)]; %Full-way model (change if neccesary)
[model,loadings,VarExp] = FitModel(hypercube,scl,IncludedFactors,plots);
disp("Explained variance: "+round(VarExp*100,2)+"%"+newline)

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