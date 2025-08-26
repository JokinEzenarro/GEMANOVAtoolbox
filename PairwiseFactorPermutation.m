function out = PairwiseFactorPermutation(X, scl, p, realModel, f, plots)
% PairwiseFactorPermutation:
% For ONE chosen factor f, perform pairwise permutation tests between ALL
% its levels (pairwise). For each pair (l1,l2), we permute ONLY those two
% levels within each stratum defined by the other factors (design-preserving),
% refit GEMANOVA, and use as test statistic the explained variance.
%
% INPUTS
%   X         : N-dimensional data array (same shape used to fit realModel)
%   scl       : cell array with levels per factor (length = nFactors)
%   p         : number of permutations
%   realModel : fitted GEMANOVA model (must contain fields: effects{f}, res, F, Fix)
%   f         : index of the factor to test (1-based)
%   plots     : 1 to produce a level x level grid of histograms; 0 otherwise
%
% OUTPUT (struct 'out')
%   out.factor_index  : f
%   out.levels        : vector of level indices for factor f
%   out.pairs         : Kx2 array with tested pairs (rows are [l1 l2])
%   out.pvals         : Kx1 vector of p-values (upper-tail)
%   out.stat_obs      : Kx1 vector of observed |u_f(l1)-u_f(l2)|
%   out.stats_perm    : cell Kx1; each cell is a p-by-1 vector of permuted stats
%   out.note          : description of the test statistic and H0
%   out.table         : summary table

    if nargin < 6, plots = 0; end

    % Basic checks
    s = size(X);
    nFactors = numel(scl);

    % Levels and pairs
    levels = scl{f}(:)';              % row vector of level indices (1..Lf or labels you use)
    Lf = numel(levels);
    pairs = nchoosek(levels, 2);      % K x 2
    K = size(pairs, 1);

    % Precompute strata (combinations of other factors)
    dims_other = setdiff(1:nFactors, f);
    if ~isempty(dims_other)
        levelGrids = arrayfun(@(d) 1:s(dims_other(d)), 1:numel(dims_other), 'UniformOutput', false);
        [levelGrids{:}] = ndgrid(levelGrids{:});
        combos = cell2mat(cellfun(@(x) x(:), levelGrids, 'UniformOutput', false)); % (#combos x numel(dims_other))
    else
        combos = zeros(1,0); % no other factors
    end

    % Observed statistics from the fitted model (Explained variance)
    SStot = sum(X(:).^2); % Total variance of original data
    SSmodel = SStot - sum(realModel.res(:).^2);
    stat_obs = SSmodel/SStot;

    % Storage
    stats_perm = cell(K,1);
    pvals      = nan(K,1);

    % Global waitbar across all pairs x permutations
    totalIters = K * p;
    iterCount  = 0;
    wb = waitbar(0, sprintf('Pairwise permutations (factor %d)', f), ...
                 'Name','PairwiseFactorPermutation');


    % Main loop over pairs
    for k = 1:K
        l1 = pairs(k,1);
        l2 = pairs(k,2);
        stats = zeros(p,1);

        for it = 1:p
            Xp = X; % copy

            % Within each stratum of the OTHER factors, randomly swap l1/l2
            for row = 1:size(combos,1)
                idx_template = repmat({':'}, 1, nFactors);
                % Fix other-factor indices for this stratum
                for kk = 1:numel(dims_other)
                    idx_template{dims_other(kk)} = combos(row, kk);
                end
                % Indices for the two levels to be permuted
                idx_l1 = idx_template; idx_l1{f} = l1;
                idx_l2 = idx_template; idx_l2{f} = l2;

                % Random swap with probability 0.5 (design-preserving permutation)
                if rand < 0.5
                    tmp = Xp(idx_l1{:});
                    Xp(idx_l1{:}) = Xp(idx_l2{:});
                    Xp(idx_l2{:}) = tmp;
                end
            end

            % Refit GEMANOVA on permuted data
            F = 1; 
            Fix = zeros(nFactors,1);
            modelP = gemanovaInit(Xp, F, Fix, scl, realModel);

            % Test statistic from permuted model (Explained variance)
            stats(it) = (SStot - sum(modelP.res(:).^2)) / SStot;

            % Update waitbar
            iterCount = iterCount + 1;
            if mod(iterCount, max(1, floor(totalIters/100))) == 0
                waitbar(iterCount/totalIters, wb, ...
                    sprintf('Pairwise permutations (factor %d):', f));
            end
        end

        % Upper-tail p-value with +1/(p) correction
        pvals(k) = (sum(stats >= stat_obs) + 1) / p;
        stats_perm{k} = stats;
    end
    close(wb)

    % Package output
    out.factor_index = f;
    out.levels       = levels;
    out.pairs        = pairs;
    out.pvals        = pvals;
    out.stat_obs     = stat_obs;
    out.stats_perm   = stats_perm;
    out.note = ['Statistic = |u_f(l1)-u_f(l2)| using mode-f loadings from GEMANOVA. ', ...
                'H0 (per pair): no difference between the two levels.'];
    out.table = table(pairs(:,1), pairs(:,2), pvals, ...
                      'VariableNames', {'Level1','Level2','pValue'});

    % Optional: level x level grid plots
    if plots == 1
        % Create a K-index map from (l1,l2) -> k
        pairIndex = containers.Map('KeyType','char','ValueType','double');
        for k = 1:K
            key = sprintf('%d_%d', pairs(k,1), pairs(k,2));
            pairIndex(key) = k;
        end

        figure('Name', sprintf('Pairwise permutation tests – Factor %d', f));
        tl = tiledlayout(Lf, Lf, 'Padding','compact', 'TileSpacing','compact');

        for r = 1:Lf
            for c = 1:Lf
                nexttile(sub2ind([Lf Lf], r, c));
                axis off; box on;

                if r == c
                    % Diagonal: label the level
                    text(0.5, 0.5, sprintf('Level %d', r), 'HorizontalAlignment','center', 'FontWeight','bold');
                elseif r < c
                    % Upper triangle: plot histogram for (r,c) if exists
                    key = sprintf('%d_%d', r, c);
                    if isKey(pairIndex, key)
                        k = pairIndex(key);
                        % Histogram
                        histogram(stats_perm{k}*100, 25); hold on;
                        xline(stat_obs*100, 'r', 'LineWidth', 1.5);
                        axis tight; box on; axis on;
                        title(sprintf('(%d,%d)  p=%.3g', r, c, pvals(k)), 'FontSize', 8);
                        xlabel('Explained variance (%)','FontSize',7); ylabel('Count','FontSize',7);
                    end
                else
                    % Lower triangle: keep empty or mirror info if desired
                    % Here: faint grid
                    text(0,0,sprintf(''));
                end
            end
        end
        title(tl, sprintf('Pairwise permutation tests – Factor %d', f), 'FontWeight','bold');
    end
end