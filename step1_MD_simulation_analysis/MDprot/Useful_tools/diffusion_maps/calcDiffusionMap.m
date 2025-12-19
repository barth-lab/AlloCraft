function [diffMapt, Lambda, D2] = calcDiffusionMap(A,options)
    % calcDiffusionMap Calculates diffusion maps from a matrix A. A could
    % be the adjacency matrix of a graph or an affinity matrix of a set of
    % data
    %% Usage:
    % diffMapt = calcDiffusionMap(K)
    % [diffMapt, Lambda, D2] = calcDiffusionMap(K,options)
    %
    %% Description:
    % * diffMapt is the dim-dimensional embedding of the graph if A is an
    % adjacency matrix or of the dataset if A is an affinity matrix
    % * Lambda are the eigenvalues of the matrix M, which represents a Laplacian
    % normalized kernel
    % * D2 squared diffusion distances
    % * A input square matrix that could be adjacency matrix of a graph 
    % or an affinity matrix of the dataset created with calcAffinityMat
    % * options: 
    %       alpha: parameter that tunes the influence of data point density
    %       on the diffusion process. For alpha = 1, the diffusion operator
    %       approximates a Laplace-Beltrami opertator. For alpha = 0.5, it
    %       approximates Fokker-Planck diffusion. For alpha = 0, it reduces
    %       to graph Laplacian normalization
    %       
    %       dim: number of dimensions of the embedding
    %
    %       t: scaling parameter, equivalent to moving the diffusion 
    %       process forward in time. Allows Discovery of different scale of
    %       underlying patterns in the data
    % 
    % see also calcAffinityMat
    %
    % Mahdi Hijazi (2023)
    arguments   
       A
       options.alpha = 0.5; % Defaults to Fokker-Planck diffusion
       options.dim = 4; % Number of dimensions 
       options.t = 1; % Scaling parameter, equivalent to moving the diffusion process forward in time
       options.g = @(x)(x.^(2)); % Parametrization function

    end

% My attempt
 D = diag(sum(A,2)+eps); % Diagonal matrix of summed elements over K
 Lalpha = sparse(D^-options.alpha*A*D^-options.alpha); 

 % Apply graph kernel normalization
 Dalpha = sparse(diag(sum(Lalpha,2)+eps));
 M = Dalpha\Lalpha;
 Mt = M;

options.isreal = true;
options.issym = true;
% Calculate options.dim + 1 eigenvalues and vectors
% Remember that the powers of the eigenvalues are the eigenvalues of M^t
if issparse(Mt)
    [v,lambda] = eigs(Mt,options.dim + 1,'largestabs',options);
    Lambda     = diag(lambda);
else
    [v,lambda] = eig(Mt);
    [lambda,I] = sort(diag(lambda),'descend'); % does eig  return the values sorted?
    Lambda = lambda(1:options.dim+1);
    v    = v(:,I(1:options.dim+1));
end

% Normalize v, Why?
% v        = v./ repmat(sqrt(sum(v.^2)),size(v,1),1);

Psi        = D^-options.alpha * v; % maybe?
inds       = 2:length(Lambda); % disregard the first trivial eigenvalue

diffMapt = v.*  repmat(Lambda'.^ options.t,size(Psi,1),1);
diffMapt = diffMapt(:,inds)';


% Calculate diffusion distances:
D2 = zeros(length(A));
for i = 1:length(A)-1
    for j = i+1:length(A)
        D2(i,j) = sum(Lambda'.^(2*options.t).* (v(i,:) - v(j,:)).^2);
%         D2(i,j) = sum( g(options.t*Lambda') .* (v(i,:) - v(j,:)).^2);

%         D2(i,j) = sum( arrayfun(@(x)g(x), options.t*Lambda') .* (v(i,:) - v(j,:)).^2);
    end
end

end