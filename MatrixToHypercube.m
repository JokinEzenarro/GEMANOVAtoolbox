function [hypercube, idx_cell, scl, mult] = MatrixToHypercube(X, F)
% Fold data matrix X into an N-dimensional hypercube based on F factor combinations
% Inputs:
%   - X: [n x p] data matrix (n samples, p variables)
%   - F: [n x k] matrix/table with k factors (numeric or categorical)
%   describing the DoE
% Outputs:
%   - hypercube: N-dimensional array with size = [levels_factor1, ..., levels_factork, p]
%   - idx_cell: linear index of each observation in the hypercube
%   - scl: levels in each factor
%   - mult: is the response multivariate?
%
% Author: Jokin Ezenarro (jokin@food.ku.dk)
% Modified: 19/07/2025

    % Convert table to matrix if needed
    if istable(F)
        F = table2array(F);
    end

    % Ensure F is double for indexing
    F = double(F);

    % Get levels per factor
    [F_unique, ~, idx] = unique(F, 'rows');

    % Determine size of the hypercube
    numFactors = size(F, 2);
    size_cube = zeros(1, numFactors);
    for i = 1:numFactors
        size_cube(i) = length(unique(F(:, i)));
    end

    % Initialize hypercube
    p = size(X, 2);
    hypercube = NaN([size_cube, p]);  % N-dimensional + variable dimension
    countCube = zeros(size_cube);     % To count replicates for averaging

    % Map each observation to its position in the hypercube
    idx_cell = num2cell(zeros(size(X,1), 1));
    for i = 1:size(X,1)
        % Find index in each dimension
        ind = zeros(1, numFactors);
        for j = 1:numFactors
            levels = unique(F(:, j), 'stable'); % assume ordered factors
            [~, ind(j)] = ismember(F(i,j), levels);
        end    

        % Linear index (not used here, but can be returned)
        idx_cell{i} = sub2ind(size_cube, ind(:));

        % Accumulate values and counts for averaging if needed
        cube_idx = num2cell(ind);  % convert to cell for indexing
        if isnan(hypercube(cube_idx{:},1))
            hypercube(cube_idx{:},:) = X(i,:);
            countCube(cube_idx{:}) = 1;
        else
            hypercube(cube_idx{:},:) = hypercube(cube_idx{:},:) + reshape(X(i,:), [ones(1, numFactors), size(X,2)]);
            countCube(cube_idx{:}) = countCube(cube_idx{:}) + 1;
        end
    end

    % Average over replicates
    for i = 1:numel(countCube)
        if countCube(i) > 1
            [subs{1:numFactors}] = ind2sub(size_cube, i);
            hypercube(subs{:},:) = hypercube(subs{:},:) ./ countCube(subs{:});
        end
    end

    % Create scl (levels in each factor)
    for i=1:size(F,2) % number of factors
        scl{i}=[1:max(F(:,i))];
    end
    
    mult = 0;
    if size(X,2) > 1 %multivariate signal
        mult = 1;
        scl{numFactors+1}=[1:size(X,2)];
    end
end