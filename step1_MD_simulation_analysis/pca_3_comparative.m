% Plot scatter using PC1 from ligand and PC1 from receptor:

pc1ProtLig = [pProt(:,1) p(:,1)]; % Contains first PC for receptor and ligand
figure; scatter(pProt(:,1),p(:,1),5,C(:,1),'filled')
hold on
for i = 1:numRuns
    % Find centroids of the scatter of each run and label it
    if i==1
       lower_lim = 1;
    else
       lower_lim = sum(nFramesEff(1:i-1))+1;
    end
    upper_lim = sum(nFramesEff(1:i));
    p1 = mean(pProt(lower_lim:upper_lim,1));
    p2 = mean(p(lower_lim:upper_lim,1));

    scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',1.5)
              hold on
    % Add text label for runs with proper offset
    text(p1+range(pProt,'all')/50,p2+range(p,'all')/50,['Run ' num2str(i)],'FontSize',14)
end
  xlabel('Receptor PC 1', 'fontsize', 25);
  ylabel('Ligand PC 1', 'fontsize', 25);
%   title(['Scatter plot of ' num2str(kPrinComp) ' receptor PCs, colored by run'])
    legend('Data colored by run','Centroids of runs')
  legend boxoff

% Save the plot
savefig([pcadir 'protlig_pca_run_scatter_' name]);
print2pdf([pcadir 'protlig_pca_run_scatter_' name]);


%% Find the best number of clusters if the user did not specify
if isempty(kClusters)
  if ~exist('kmax', 'var')  % Default kmax to square root of total frames
      kmax = round(sqrt(nframes));
  elseif isempty(kmax)
      kmax = round(sqrt(nframes));
  end
  % Matlab's function for evaluating optimal clustering
  % Evaluate from 1:sqrt(nframes) clusters
%   E = evalclusters((trj),'kmeans','CalinskiHarabasz','klist',1:round(sqrt(nframes)));
  E = evalclusters(real(p(:,1:kPrinComp)),'kmeans','CalinskiHarabasz','klist',1:kmax);
  kClusters = E.OptimalK;
end

 [indexOfCluster_pcaProtLig, centroid_pcaProtLig] = kmeans(pc1ProtLig, kClusters);
 
%  % Find the centers  (frames that are closest to the centroids):
%  ind_centers = zeros(kClusters,1);
%  
% for k=1:kClusters % For every centroid
%  [~,ind_centers(k)] = min(vecnorm(p(:,1:kPrinComp)-centroid_pca(k,:),2,2));  
% end

% Plot the data
  figure
  scatter(pc1ProtLig(:, 1), pc1ProtLig(:, 2), 50, indexOfCluster_pcaProtLig, 'filled');
  xlabel('Receptor PC 1', 'fontsize', 25);
  ylabel('Ligand PC 1', 'fontsize', 25);
  title(['Clustering with ' num2str(kPrinComp) ' PCs'])
  
  hold on
  % Plot centroids:
  scatter(centroid_pcaProtLig(:,1),centroid_pcaProtLig(:,2),60,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
  legend('Data','Centroids')
  legend boxoff
  
  for i = 1:size(centroid_pcaProtLig,1)
     % Add text label for runs with proper offset
    text(centroid_pcaProtLig(i,1)+range(pProt,'all')/50,centroid_pcaProtLig(i,2)+range(p,'all')/50,['C' num2str(i)],'FontSize',14) 
  end