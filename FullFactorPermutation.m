function pValues = FullFactorPermutation(X, scl, p, realModel, plots, mult)
% FullFactorPermutation:
% For each factor (dimension except the last one), performs a permutation
% test by permuting only that factor and computing the explained variance.
%
% INPUTS:
%   X     - N-dimensional data array
%   scl   - cell array containing the levels in each factor
%   p     - Number of permutations
%   realModel - GEMANOVA model to test against the permutations
%   plots - 1 to show histograms of variance explained
%   mult  - Is it a multivariate signal? (1=Yes 0=No)
%
% OUTPUT:
%   pValues - vector of p-values corresponding to each factor

    if nargin < 5
        plots = 0;
    end

    s = size(X);              % Dimensions of the data cube
    nFactors = length(scl);   % Number of factors        

    % Parameters for the GEMANOVA function
    F = 1;
    Fix(:,1) = zeros(1,nFactors);

    % Original variance explained, for reference
    SStot = sum(X(:).^2); % Total variance of original data
    SSmodel = SStot - sum(realModel.res(:).^2);
    realVarExp = SSmodel/SStot;

    pValues = zeros(1, nFactors);          % p-values for each factor
    permutedVarExp = zeros(p, nFactors);   % Store explained var for each perm

    waitHandle = waitbar(0, 'Running factor-wise permutations...');

    for f = 1:nFactors-mult  % for each factor
        dims = setdiff(1:nFactors, f);  % other factors
        levelGrids = cell(1, numel(dims));

        % Create index grids for all combinations of other factor levels
        for d = 1:numel(dims)
            levelGrids{d} = 1:s(dims(d));
        end
        [levelGrids{:}] = ndgrid(levelGrids{:});
        combos = cell2mat(cellfun(@(x) x(:), levelGrids, 'UniformOutput', false));

        for i = 1:p
            Xp = X;  % Copy original data

            for row = 1:size(combos, 1)
                % Build fixed index for each other factor
                fixed_idx = combos(row, :);  % levels of other factors
                idx_template = repmat({':'}, 1, nFactors);  % full index

                for l = 1:numel(dims)
                    idx_template{dims(l)} = fixed_idx(l);  % fix value
                end

                % Permute current factor `f` within this fixed combination
                perm_order = randperm(s(f));
                for pos = 1:s(f)
                    idx_template{f} = perm_order(pos);  % current permuted value
                    idx_template_target = idx_template;
                    idx_template_target{f} = pos;       % write into correct location

                    % Copy the slice
                    Xp(idx_template_target{:}) = X(idx_template{:});
                end
            end

            % Evaluate permuted model
            modelP = gemanovaInit(Xp, F, Fix, scl, realModel);
            permutedVarExp(i, f) = (SStot - sum(modelP.res(:).^2)) / SStot;

            waitbar(((f-1)*p+i) / (nFactors*p), waitHandle);
        end

        % Compute p-value
        pValues(f) = (sum(permutedVarExp(:, f) >= realVarExp)+1) / p;
    end

    close(waitHandle);

    % Optional plotting
    if plots == 1
        figure;
        for f = 1:nFactors-mult
            nexttile;
            histogram(permutedVarExp(:, f)*100,25);
            % histogram(permutedVarExp(:, f)*100,[0:4:100]);
            hold on;
            % xlim([0 100])
            % ylim([0 p])
            xline(realVarExp*100, 'r', 'LineWidth', 2)
            title(['Factor ' num2str(f) ' (p=' num2str(pValues(f)) ')'])
            xlabel('Explained variance (%)')
            ylabel('Count')
        end
    end
end