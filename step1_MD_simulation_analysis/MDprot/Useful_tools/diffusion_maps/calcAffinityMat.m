function [A, D, I] = calcAffinityMat(X, options)
    % Calculate affinity matrix A from an input dataset X using k nearest
    % neighbors.
    % 
    %% Usage:
    % [A] = calcAffinityMat(X, options)
    %
    %% Description:
    % * X input dataset with size N x M, where N is the number of
    % observations and M is the number of dimensions
    % * A affinity matrix 
    % * D distances between observation i and its k nearest neighbors (k x
    % N)
    % * I indices of observation i and its k nearest neighbors (k x N)
    %
    % * options:
    %       kNN: number of nearest neighbors, play with this number to
    %       capture the relationship between datapoints, defaults to 5
    %
    %       distanceMetric: distance metric used to check for nearest
    %       neighbors, defaults to 'euclidean'. Any options available for
    %       pdist2 can be used
    %
    %       selfTuning: Apply local scaling as in: Zelink-Manor and Perona, 
    %       Self-Tuning Spectral Clustering, 2004. Defaults to 0. The input
    %       number is the nearest neighbor for which the local scaling will
    %       be computed, a recommended value is 7
    %
    % see also calcDiffusionMap
    arguments
        X
        options.kNN = 5;
        options.distanceMetric = 'euclidean'
        options.selfTuning = 0
    end

    [M,~] = size(X);
    
    % Get nearest neighbor list
    [D,I] = pdist2(X,X,options.distanceMetric,'Smallest',options.kNN);
    

    % Get the sparse row and column indices
    rowInds = repmat((1:M),options.kNN,1);
    rowInds = rowInds(:);
    colInds = double(I(:));
    vals    = D(:);
    
    if options.selfTuning % Self-Tuning? 
        % Local scale is simply σi = d(si, sK), where K is the Kth nearest
        % neighbor, Zelink-Manor and Perona used K = 7, which provided good
        % results
        counter = 1;
        sigij = zeros(options.kNN * M,1);
        nnAutotune = min(options.selfTuning,size(D,1));
        sigmaKvec = (D(nnAutotune,:)); % distance from i to K
        for i = 1:M
            sigij(counter : counter + options.kNN - 1) = sigmaKvec(i) * sigmaKvec(I(:,i));
            counter = counter + options.kNN;
        end
    end

if  options.selfTuning % Self-Tuning? 
    kernel = exp(-vals.^2./(sigij+eps));
    K = sparse(rowInds, colInds, kernel, M, M);
else
    sig = median(vals);
    kernel = exp(-vals.^2/sig^2);
    K = sparse(rowInds, colInds, kernel, M, M);
end
    A = (K + K')/2;

end
