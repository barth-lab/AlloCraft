function [indexOfCluster_pca_sorted, centroid_pca_sorted, p, ind_centers] = cluster_traj_pca(trj, kPrinComp, kclusters, find_mean, kmax)
%cluster_traj_pca Cluster trajectories by performing kmeans on the PCA
% The output clusters are sorted by size in descending order
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
%
%   Usage:
%   [indexOfCluster_pca] = cluster_traj_pca(trj, kPrinComp)
%   [indexOfCluster_pca, centroid_pca] = cluster_traj_pca(trj, kPrinComp)
%   [indexOfCluster_pca, centroid_pca] = cluster_traj_pca(trj, kPrinComp, kclusters)
%   [indexOfCluster_pca, centroid_pca, p, ind_centers] = cluster_traj_pca(trj, kPrinComp, kclusters, find_mean, kmax)
%
% * trj: is the input trajectory: [nframes x 3natoms]
% * kPrinComp: number of principal components to be used in clustering.
% The function will visualize up to 4 principal components (4th dimension
% being size of the scatter points), but clustering can be performed on
% higher kPrinComp
% * kclusters: number of clusters, if not provided, the default will be
% evaluated using "evalclusters" function, note that this will make the
% function slower
% * find_mean: if this variable equals to 1, mean structure will be determined
% * kmax: maximum k for the algorithm to consider, defaults to square root
% of the number of frames in the trajectory
% * indexOfCluster_pca: indices of the cluster: [nframes x 1]
% * centroid_pca: the centroids of the clusters, kclusters x kPrinComp
% * p (projection):  principal components (projection of the trajectory on to principal modes) [nframes x 3natoms]
% * ind_centers: the index of the frame closest to the centroid of each cluster

if exist('find_mean', 'var') 
    if find_mean == 1
        % Find the mean structure of the trajectory
        [~, trj] = meanstructure(trj);
    end
end

% Calc the PCA and do the clustering:
[p, ~, ~] = calcpca(trj);
[~, ~, ~, ~, explained] = pca(trj); % get %ge explained by every PC

% Find the best number of clusters if the user did not specify
if ~exist('kclusters', 'var') || isempty(kclusters)
    
  nframes = size(trj,1);
  if ~exist('kmax', 'var')  % Default kmax to square root of total frames
      kmax = round(sqrt(nframes));
  elseif isempty(kmax)
      kmax = round(sqrt(nframes));
  end
  % Matlab's function for evaluating optimal clustering
  % Evaluate from 1:sqrt(nframes) clusters
%   E = evalclusters((trj),'kmeans','CalinskiHarabasz','klist',1:round(sqrt(nframes)));
  E = evalclusters(real(p(:,1:kPrinComp)),'kmeans','CalinskiHarabasz','klist',1:kmax);
  kclusters = E.OptimalK;
end

 [indexOfCluster_pca, centroid_pca] = kmeans(p(:,1:kPrinComp), kclusters,'Replicates',5,'Display','final');
 
  cluster_sizes = zeros(kclusters,1);
 % Find sizes of clusters
 for k=1:kclusters 
   cluster_sizes(k) = length(find(indexOfCluster_pca==k));
 end
[cluster_sizes_sorted,cluster_size_ndx] = sort(cluster_sizes,'descend');
% Now sort indexOfCluster_pca and centroid_pca according to the criterion
indexOfCluster_pca_sorted = zeros(length(indexOfCluster_pca),1);
centroid_pca_sorted =  centroid_pca(cluster_size_ndx,:);
% Fill the sorted clusters in the index file
 for k=1:kclusters
      ndx = find(indexOfCluster_pca==cluster_size_ndx(k)); % this will find the kth cluster by size
      indexOfCluster_pca_sorted(ndx) = k;
 end
 
 % Find the centers  (frames that are closest to the centroids):
 ind_centers = zeros(kclusters,1);
 
for k=1:kclusters % For every centroid
 [~,ind_centers(k)] = min(vecnorm(p(:,1:kPrinComp)-centroid_pca_sorted(k,:),2,2));  
end

pcaRange = range(p,'all');

% Plot the data
  figure
  if kPrinComp == 2 % 2D plot for 2 principal components
    scatter(p(:, 1), p(:, 2), 15, indexOfCluster_pca_sorted, 'filled');
  elseif kPrinComp == 3 % 3D plot for 3 principal components
    scatter3(p(:, 1), p(:, 2), p(:, 3),10, indexOfCluster_pca_sorted, 'filled', ...
        'MarkerFaceAlpha',0.5);
    zlabel(['PC 3, ' num2str(explained(3),'%.2f') ' explained'], 'fontsize', 25);
  else % 3D plot for 3+ principal components, where 4th component is size of the scatter points
     % Scale size from 5 to 50
     sizeFcn = 45/range(p(:, 4))*(p(:, 4) + abs(min(p(:, 4)))) + 5;
     scatter3(p(:, 1), p(:, 2), p(:, 3),sizeFcn, indexOfCluster_pca_sorted, 'filled', ...
        'MarkerFaceAlpha',0.5);
     zlabel(['PC 3, ' num2str(explained(3),'%.2f') ' explained'], 'fontsize', 25); 
  end
  xlabel(['PC 1, ' num2str(explained(1),'%.2f') ' explained'], 'fontsize', 25);
  ylabel(['PC 2, ' num2str(explained(2),'%.2f') ' explained'], 'fontsize', 25);
  title(['Clustering with ' num2str(kPrinComp) ' PCs, ' num2str(sum(explained(1:kPrinComp)),'%.2f') ' explained total'], 'fontsize', 25)
  
  hold on
  % Plot centroids:
  if kPrinComp == 2
    scatter(centroid_pca_sorted(:,1),centroid_pca_sorted(:,2),60,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
    for i = 1:size(centroid_pca_sorted,1)
     % Add text label for runs with proper offset
    text(centroid_pca_sorted(i,1)+pcaRange/50,centroid_pca_sorted(i,2)+pcaRange/50,['C' num2str(i)],'FontSize',14) 
    end
  else
    scatter3(centroid_pca_sorted(:,1),centroid_pca_sorted(:,2),centroid_pca_sorted(:,3),...
              60,'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
    for i = 1:size(centroid_pca_sorted,1)
     % Add text label for runs with proper offset
    text(centroid_pca_sorted(i,1)+pcaRange/50,centroid_pca_sorted(i,2)+pcaRange/50, ...
        centroid_pca_sorted(i,3)+pcaRange/50,['C' num2str(i)],'FontSize',14) 
    end
  end
  legend('Data','Centroids')
  legend boxoff
  
  
end

