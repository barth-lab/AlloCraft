%% Do PCA for the receptor coordinates
% Or distances
% Or dihedrals
% To do: add highest density structures and plots similar to ligand pca


% Trajectory type:
if strcmp(trjType,'CA')
    atomIndices = mainChain.getAtoms();
    traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
elseif strcmp(trjType,'StructuredCA')
    if ~isempty(settings.helices)
    helixRes = [];
    for i = 1:size(settings.helices,1)
        helixRes = [helixRes settings.helices(i,1):settings.helices(i,2)];
    end
    
    end
    atomIndices = mainEntry.getAtoms('Chain', mainChain.index, 'Name', 'CA','Residues',helixRes);
    traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
elseif strcmp(trjType,'DihAll')
    traj = mainSim.dihedralsMat;
elseif strcmp(trjType,'DihBB')
    traj = mainSim.dihedralsMat(:,(mainSim.reSort(:,2)==1 | mainSim.reSort(:,2)==2));
elseif strcmp(trjType,'DihSC')
    traj = mainSim.dihedralsMat(:,(mainSim.reSort(:,2)==0));
elseif strcmp(trjType,'Distances') % a bit more complicated

     atomIndices = mainChain.getAtoms();
     traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
     resDisCalc = 1:CGPCA:size(traj,2)/3; % Take every Nth residue
     distMat = zeros(length(resDisCalc),length(resDisCalc),size(traj,1));
    
     for frameHere = 1:size(traj,1)
        distMat(:,:,frameHere) = calcdistancematrix(traj(frameHere,to3(resDisCalc)));
     end
     temp = reshape(distMat,length(resDisCalc)*length(resDisCalc),size(traj,1));
    stdPairs = std(temp');
    [a,b] = sort(stdPairs,'descend');   % Take M distances with highest std
    traj = temp(b(1:min(NdisPCA,length(b))),:)'; % Number of distances depends on protein and application
    clear distMat temp;
end

if ~exist('kPrinComp','var')
     kPrinComp = 2; % Number of principal components to take into consideration
end
[indexOfCluster_pcaProt, centroid_pcaProt, pProt, ind_centersProt] = cluster_traj_pca(traj, kPrinComp, settings.kClusters, [], settings.kmax);
clear traj;

figPath = fullfile(pcadir, "prot_pca_scatter_" + trjType + "_" +  mainEntry.name+ ".fig");
savefig(figPath);
print2pdf(figPath);

add2log(md2pathdir, "PCA clustering of receptor performed!");


%% Plot PCA1/2 colored by run and by cluster

figure;
pcaRange = range(pProt, 'all');

scatter(pProt(:,1),pProt(:,2),5,C(:,1),'filled');
hold on;

for i = 1:mainSim.runCount
    % Find centroids of the scatter of each run and label it
    if i==1
       lower_lim = 1;
    else
       lower_lim = sum(nFramesEff(1:i-1))+1;
    end
    upper_lim = sum(nFramesEff(1:i));
    p1 = mean(pProt(lower_lim:upper_lim,1));
    p2 = mean(pProt(lower_lim:upper_lim,2));

    scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',1.5)
              hold on
    % Add text label for runs with proper offset
    text(p1+pcaRange/50,p2+pcaRange/50,['Run ' num2str(i)],'FontSize',14)
end
  xlabel('PC 1', 'fontsize', 25);
  ylabel('PC 2', 'fontsize', 25);
  title([num2str(kPrinComp) ' receptor PCs from ' trjType])
    legend('Data colored by run','Centroids of runs')
  legend boxoff

% Save the plot
figPath = fullfile(pcadir, "prot_pca_run_scatter_"+ trjType + "_" +  mainEntry.name+ ".fig");
savefig(figPath);
print2pdf(figPath);

%% Now get the representative structures:
% ind_centersProt contains the frame numbers of the cluster centers

relFrameIndices = zeros(length(indexOfCluster_pcaProt), 1);

for Ci = 1:length(ind_centersProt)
    filter = (indexOfCluster_pcaProt == Ci);
    relFrameIndices(filter) = 1:sum(filter);
end

bfactor = zeros(length(ind_centersProt), mainEntry.atomCount);
pca_frame_centers = zeros(length(ind_centersProt),3);
protpcaTraj = zeros(length(ind_centersProt),mainEntry.atomCount * 3);
highestDenTrajProt = zeros(length(ind_centersProt),mainEntry.atomCount * 3);
for thisCenter = 1:length(ind_centersProt)
    % Grab the run that the center belongs to:
    runCenter = C(ind_centersProt(thisCenter)); %Misuse color vector to grab the run
    % Make a workaround for finding lower limit with trajectories of
    % unequal length
    if runCenter==1
       lower_lim_center = 1;
    else
       lower_lim_center = sum(nFramesEff(1:runCenter-1))+1;
    end
    pca_frame_centers(thisCenter,1) = ind_centersProt(thisCenter) - lower_lim_center + (settings.frames2skip + 1);
    pca_frame_centers(thisCenter,2) = runCenter;
    pca_frame_centers(thisCenter,3) = length(find(indexOfCluster_pcaProt==thisCenter));
    protpcaTraj(thisCenter,:) = mainSim.traj{runCenter}(pca_frame_centers(thisCenter,1),:);

%     sim = mainSim.createSubset(C(indexOfCluster_pcaProt == thisCenter, :));
%     ref = sim.traj{1}(relFrameIndices(ind_centersProt(thisCenter)), :);
%     atoms = mainChain.getAtoms();
%     
%      % Add mean displacement from cluster center in place of B-factor
%     diff = reshape(sim.traj{1} - ref, [], mainEntry.atomCount, 3);
%     meanDist = mean(sqrt(sum(diff.^2, 3)), 1);
% 
%     bfactor(thisCenter, atoms) = meanDist(atoms);
%     
%     if mainEntry.hasChain(Chains.ligand)
%         ligandAtoms = ligandChain.getLigandAtoms();
%         bfactor(thisCenter, ligandAtoms) = meanDist(ligandAtoms);
%     end
%     % Write pdb of this center with displacement
%     writepdbndx(fullfile(pcadir, "ligandReceptorpca_"  + name + "_C" + num2str(thisCenter) + ".pdb"), mainEntry.pdb, [], 0, protpcaTraj(thisCenter,:), 'default', bfactor(thisCenter,:));

    % Calculate highest density point and find structure closest to it:
    [valuesCluster, centersCluster] = hist3([pProt(indexOfCluster_pcaProt==thisCenter, 1) pProt(indexOfCluster_pcaProt==thisCenter, 2)],[51 51]); % bins may need modification
    [a,b] = max(valuesCluster,[],'all'); % get highest density point in linear space
    [temp1, temp2] = ind2sub([51,51],b);  % get i and j indices
    
    highestDen = [centersCluster{1}(temp1) centersCluster{2}(temp2)]; % that's the highest density bin in PC space
    
     pTemp = pProt(indexOfCluster_pcaProt==thisCenter,1:kPrinComp);
        % The closest point to the highest density point
       [~,indHighestDen] = min(vecnorm(pTemp-highestDen,2,2));  
       
      findpAll = find(pProt(:,1:kPrinComp)==pTemp(indHighestDen,:));
      findpAll = findpAll(1);
    
      runFrameHighestDen = C(findpAll,:); % run and frame of highest density point!
      highestDenTrajProt(thisCenter,:) = mainSim.traj{runFrameHighestDen(1)}(runFrameHighestDen(2),:);

      clear sim;
    
end

% Save frame center data and write out pdbs:
writematrix(pca_frame_centers, fullfile(pcadir, "prot_pcaCenterFrame_run_Nelements.txt"), 'Delimiter','space');
writepdbndx(fullfile(pcadir, sprintf("prot_pcaLigandReceptor_%s_%s.pdb",trjType ,mainEntry.name)), mainEntry.pdb, [], 0, protpcaTraj);
writepdbndx(fullfile(pcadir, sprintf("prot_highestDenpca_%s_%s.pdb",trjType ,mainEntry.name)), mainEntry.pdb, [], 0, highestDenTrajProt);

save(saveVarName,'highestDenTrajProt','-append'); % Required for pathway calculations later

add2log(md2pathdir,{'frame_centers saved to prot_pcaCenterFrame_run.txt','pdbs saved to prot_pcaLigandReceptor.pdb'});


%% Calculate RMSD matrix between centers

rmsdLigandpcaCenter = zeros(length(ind_centersProt));
seleRMSD = mainEntry.getAtoms('Chain', mainChain.index, 'NoName', 'H*');

for i =1:length(ind_centersProt)-1
    for j = i+1:length(ind_centersProt)
        rmsdLigandpcaCenter(i,j) = calcrmsd(protpcaTraj(i,to3(seleRMSD)),protpcaTraj(j,to3(seleRMSD)));
    end
end

figure;
hm = heatmap(rmsdLigandpcaCenter,'Colormap',parula);
hm.Title = ['Receptor RMSD [' char(197) '] matrix PCA centers '];
hm.XLabel = 'PCA Cluster';
hm.YLabel = 'PCA Cluster';
figPath = fullfile(pcadir, "prot_pcarmsdMatrix_"+ trjType + "_" +  mainEntry.name);
if ~endsWith(figPath, '.fig')
    figPath = figPath+'.fig';
end
savefig(figPath);
print2pdf(figPath);

%% Find representative points that have 1- highest density AND 2- contributions from the most runs
% Find the density by histogramming the data into 2D bins
figure
[values, centers] = hist3([pProt(:, 1) pProt(:, 2)],[51 51]); % bins may need modification
imagesc(centers{:},values.')
colorbar
axis xy
xlabel('PCA 1', 'fontsize', 25);
ylabel('PCA 2', 'fontsize', 25);
title(['Density plot of the first ' num2str(kPrinComp) ' PCs'])

hold on
% Plot centroids with text and legend:
% Don't forget that you're a LEGEND:
legend_text = cell(size(centroid_pcaProt,1),1);
  for i = 1:size(centroid_pcaProt,1)
    scatter(centroid_pcaProt(:,1),centroid_pcaProt(:,2),60,'MarkerEdgeColor',[0 .5 .5],...
          'MarkerFaceColor',[0.8500 0.3250 0.0980],...
          'LineWidth',1.5)
     % Add text label for runs with proper offset
    text(centroid_pcaProt(i,1)+pcaRange/50,centroid_pcaProt(i,2)+pcaRange/50,['C' num2str(i)],'FontSize',14)
    legend_text{i} = ['C' num2str(i) ', size: ' num2str( pca_frame_centers(i,3))];
  end
legend(legend_text','location','EastOutside','fontsize', 16)
set(gcf, 'Position',  [100, 100, 700, 400]) % Resize window

figPath = fullfile(pcadir, "prot_pcaDensity_"+ trjType + "_" +  mainEntry.name);
if ~endsWith(figPath, '.fig')
    figPath = figPath + ".fig";
end

savefig(figPath);
print2pdf(figPath);