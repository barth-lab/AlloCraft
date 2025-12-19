function [indexOfCluster_pca, centroid_pca, p, ind_centers, pca_frame_centers, rmsdpcaCenter] = calcPlotpca(mainSim, C, settings, options)
  arguments
    mainSim
    C
    settings
    options.SaveName = ""
    options.SavePath
    options.ChainNdx = 1; % chooses first chain if none are entered
    options.isLigand = 0; % Specific options for ligands
    options.receptorLigandResIds = [] % Important if settings.trjTypeLigand = 'DistancesCoord' or 'VectorsCoord'
    options.saveVarName % save output to variable
    options.customAtomNdx % Custom atom indices to & with other choices
    options.ColorMap %
  end

% This module plots the PCA colored by run as well as a density plot


% Input preparation:
% C is now initialized in md2pathMain
mainEntry = mainSim.entry;
receptorLigandResIds = options.receptorLigandResIds ;
trjTypeLigand = settings.trjTypeLigand;
trjType = settings.trjType;
%% Plot PCA1/2/3/4 colored by run and by cluster


if options.isLigand % Ligand
trjType = trjTypeLigand;
atomIndices = mainEntry.chains{options.ChainNdx}.getLigandAtoms();

if isfield(options, 'customAtomNdx') 
    atomIndices = intersect(options.customAtomNdx ,atomIndices);
end

traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);

if strcmp(trjType, 'DistancesCoord') || strcmp(trjType, 'VectorsCoord')
    if isempty(options.receptorLigandResIds)
        warning(['settings.trjTypeLigand is set to ' trjType  ... 
            ', but receptorLigandResIds is empty!!! Will not output meaningful results']);
    end
    % Ligand atom choice
    ndxLigand = ~selectname(mainEntry.pdb.name,'H*') & selectname(mainEntry.pdb.chainid,settings.chains(options.ChainNdx));
    % Which receptor residues to consider?
    % receptorLigandResIds
    
    % Calculate distances
    % Calculate minimum distance between any residue pair and ligand
    concatTraj = mainSim.concatRuns('StartFrame', settings.frames2skip + 1);
    nFrames = size(concatTraj,1);
    minDis = zeros(nFrames,length(receptorLigandResIds));
    minDisVector = zeros(nFrames,3*length(receptorLigandResIds));
    
    for i = 1:length(receptorLigandResIds)
        resHere = receptorLigandResIds(i);
        ndx1 = selectid(mainEntry.pdb.resseq,resHere) & ~selectname(mainEntry.pdb.name,'H*') & selectname(mainEntry.pdb.chainid,settings.chains(1));
    
        % calcMinDis function seems to be working properly, I tested it
        % against VMD via inspection of selected BB-SC and SC-SC
        % interactions
        [minDis(:,i), minDisVector(:,(3*i-2:3*i))] = calcMinDis(concatTraj,ndxLigand,ndx1); 
    end
if strcmp(trjType, 'DistancesCoord')
    traj = [traj minDis]; % Add minimum distance between ligand binding residues and ligand
elseif strcmp(trjType, 'VectorsCoord')
    traj = [traj minDisVector]; % Add min. dis. vectors
end
end

else % Other trajectory type

if strcmp(trjType,'CA')
    atomIndices = mainEntry.chains{options.ChainNdx}.getAtoms();

    if isfield(options, 'customAtomNdx') 
    atomIndices = intersect(options.customAtomNdx ,atomIndices);
    end

    traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
elseif strcmp(trjType,'StructuredCA')
    if ~isempty(settings.helices)
    helixRes = [];
    for i = 1:size(settings.helices,1)
        helixRes = [helixRes settings.helices(i,1):settings.helices(i,2)];
    end
    
    end
    atomIndices = mainEntry.getAtoms('Chain', options.ChainNdx, 'Name', 'CA','Residues',helixRes);
    if isfield(options, 'customAtomNdx') 
    atomIndices = intersect(options.customAtomNdx ,atomIndices);
    end

    traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
elseif strcmp(trjType,'DihAll')
    traj = mainSim.dihedralsMat;
elseif strcmp(trjType,'DihBB')
    traj = mainSim.dihedralsMat(:,(mainSim.reSort(:,2)==1 | mainSim.reSort(:,2)==2));
elseif strcmp(trjType,'DihSC')
    traj = mainSim.dihedralsMat(:,(mainSim.reSort(:,2)==0));
elseif strcmp(trjType,'Distances') % a bit more complicated
     NdisPCA = settings.NdisPCA;
     CGPCA = settings.CGPCA;

     atomIndices =  mainEntry.chains{options.ChainNdx}.getAtoms();
     
     if isfield(options, 'customAtomNdx') 
     atomIndices = intersect(options.customAtomNdx ,atomIndices);
     end

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
end
%%
kPrinComp = settings.kPrinComp;
[indexOfCluster_pca, centroid_pca, p, ind_centers] = cluster_traj_pca(traj, kPrinComp, settings.kClusters, [], settings.kmax);
if isfield(options, 'saveVarName')
save(options.saveVarName,'indexOfCluster_pca','centroid_pca','p','ind_centers','-append');
end
clear traj;

if isfield(options, 'SavePath')
    figPath = fullfile(options.SavePath, options.SaveName + "pca_scatter_" + trjType + "_" + mainEntry.name);
    savefig( figPath + ".fig");
    print2pdf(figPath);
end
% add2log(md2pathdir, "PCA clustering of ligand pose performed!");

%%
nFramesEff = zeros(mainSim.runCount,1);
for i = 1:mainSim.runCount
    % Take into consideration runs with different number of frames
    nFramesEff(i) = (size(mainSim.traj{i},1) - settings.frames2skip - 1)+1; % Frames used in calculations
end

pcaRange = range(p, 'all');

figure
scatter(p(:,1),p(:,2),5,C(:,1),'filled')
hold on
for i = 1:mainSim.runCount
    % Find centroids of the scatter of each run and label it
    if i==1
       lower_lim = 1;
    else
       lower_lim = sum(nFramesEff(1:i-1))+1;
    end
    upper_lim = sum(nFramesEff(1:i));
    p1 = mean(p(lower_lim:upper_lim,1));
    p2 = mean(p(lower_lim:upper_lim,2));

    scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',1.5)
              hold on
    % Add text label for runs with proper offset
    text(p1+pcaRange/50,p2+pcaRange/50,['Run ' num2str(i)],'FontSize',14)
end
  xlabel('PC 1', 'fontsize', 25);
  ylabel('PC 2', 'fontsize', 25);
  title(['Scatter plot of ' num2str(kPrinComp) ' PCs, colored by run'])
    legend('Data colored by run','Centroids of runs')
  legend boxoff

% Save the plot
if isfield(options, 'SavePath')
figPath = fullfile(options.SavePath,  options.SaveName + "pca_run_scatter_"  + trjType + "_" +  mainEntry.name);
savefig(figPath + ".fig");
print2pdf(figPath);
end

%% Now get the representative structures:
% ind_centers contains the frame numbers of the cluster centers

relFrameIndices = zeros(length(indexOfCluster_pca), 1);

for Ci = 1:length(ind_centers)
    filter = (indexOfCluster_pca == Ci);
    relFrameIndices(filter) = 1:sum(filter);
end

bfactor = zeros(length(ind_centers), mainEntry.atomCount);
pca_frame_centers = zeros(length(ind_centers),3);
pcaTraj = zeros(length(ind_centers),mainEntry.atomCount * 3);
highestDenTraj = zeros(length(ind_centers),mainEntry.atomCount * 3);
for thisCenter = 1:length(ind_centers)
    % Grab the run that the center belongs to:
    runCenter = C(ind_centers(thisCenter)); %Misuse color vector to grab the run
    % Make a workaround for finding lower limit with trajectories of
    % unequal length
    if runCenter==1
       lower_lim_center = 1;
    else
       lower_lim_center = sum(nFramesEff(1:runCenter-1))+1;
    end
    pca_frame_centers(thisCenter,1) = ind_centers(thisCenter) - lower_lim_center + (settings.frames2skip + 1);
    pca_frame_centers(thisCenter,2) = runCenter;
    pca_frame_centers(thisCenter,3) = length(find(indexOfCluster_pca==thisCenter));
    pcaTraj(thisCenter,:) = mainSim.traj{runCenter}(pca_frame_centers(thisCenter,1),:);

%     sim = mainSim.createSubset(C(indexOfCluster_pca == thisCenter, :));
%     ref = sim.traj{1}(relFrameIndices(ind_centers(thisCenter)), :);
%     atoms = mainChain.getAtoms();
%     ligandAtoms = ligandChain.getLigandAtoms();
%     
%      % Add mean displacement from cluster center in place of B-factor
%     diff = reshape(sim.traj{1} - ref, [], mainEntry.atomCount, 3);
%     meanDist = mean(sqrt(sum(diff.^2, 3)), 1);
% 
%     bfactor(thisCenter, atoms) = meanDist(atoms);
%     bfactor(thisCenter, ligandAtoms) = meanDist(ligandAtoms);
%     % Write pdb of this center with displacement
%     writepdbndx(fullfile(options.SavePath, "ligandReceptorpca_"  + name + "_C" + num2str(thisCenter) + ".pdb"), mainEntry.pdb, [], 0, ligandpcaTraj(thisCenter,:), 'default', bfactor(thisCenter,:));

    % Calculate highest density point and find structure closest to it:
    [valuesCluster, centersCluster] = hist3([p(indexOfCluster_pca==thisCenter, 1) p(indexOfCluster_pca==thisCenter, 2)],[51 51]); % bins may need modification
    [a,b] = max(valuesCluster,[],'all'); % get highest density point in linear space
    [temp1, temp2] = ind2sub([51,51],b);  % get i and j indices
    
    highestDen = [centersCluster{1}(temp1) centersCluster{2}(temp2)]; % that's the highest density bin in PC space
    
     pTemp = p(indexOfCluster_pca==thisCenter,1:kPrinComp);
        % The closest point to the highest density point
       [~,indHighestDen] = min(vecnorm(pTemp-highestDen,2,2));  
       
      findpAll = find(p(:,1:kPrinComp)==pTemp(indHighestDen,:));
      findpAll = findpAll(1);
    
      runFrameHighestDen = C(findpAll,:); % run and frame of highest density point!
      highestDenTraj(thisCenter,:) = mainSim.traj{runFrameHighestDen(1)}(runFrameHighestDen(2),:);

      clear sim;
    
end

% Save frame center data and write out pdbs:
if isfield(options, 'SavePath')
writematrix(pca_frame_centers,fullfile(options.SavePath, options.SaveName  + trjType + "_" + "pcaCenterFrame_run_Nelements.txt"),'Delimiter','space')
% writepdbndx(fullfile(options.SavePath, "ligandpca_" + name + ".pdb"), mainEntry.pdb, ligandChain.index, 0, pcaTraj);
writepdbndx(fullfile(options.SavePath, sprintf("%spcaLigandReceptor_%s_%s.pdb",options.SaveName,trjType ,mainEntry.name)), mainEntry.pdb, [], 0, pcaTraj,'default', bfactor);
writepdbndx(fullfile(options.SavePath, sprintf("%shighestDenpca_%s_%s.pdb",options.SaveName,trjType ,mainEntry.name)), mainEntry.pdb, [], 0, highestDenTraj);
end
% add2log(md2pathdir,{'frame_centers saved to pcaCenterFrame_runxt','pdbs saved to ligandpca.pdb and ligandReceptorpca.pdb'});

if isfield(options, 'saveVarName')
save(options.saveVarName,'highestDenTraj','-append'); % Required for pathway calculations later
end
%% Calculate RMSD matrix between centers:

rmsdpcaCenter = zeros(length(ind_centers));
seleRMSD = mainEntry.getAtoms('Chain', options.ChainNdx, 'NoName', 'H*');

for i =1:length(ind_centers)-1
    for j = i+1:length(ind_centers)
        rmsdpcaCenter(i,j) = calcrmsd(pcaTraj(i,to3(seleRMSD)),pcaTraj(j,to3(seleRMSD)));
    end
end

figure
hm = heatmap(rmsdpcaCenter,'Colormap',parula);
hm.Title = "RMSD [" +  char(197) + "] matrix " + options.SaveName  + trjType + " " +  "PCA centers ";
hm.XLabel = 'PCA Cluster';
hm.YLabel = 'PCA Cluster';

if isfield(options, 'SavePath')
figPath = fullfile(options.SavePath, options.SaveName  + trjType + "_" +  "pcaLigandrmsdMatrix");
savefig(figPath);
print2pdf(figPath);
end


%% Plot the clustered data
%   figure
%   scatter(p(:, 1), p(:, 2), 10, indexOfCluster_pca, 'filled');
%   xlabel('PCA 1', 'fontsize', 25);
%   ylabel('PCA 2', 'fontsize', 25);
%   title(['Clustering with ' num2str(kPrinComp) ' PCs'])
%
%   hold on
%   % Plot centroids:
%   scatter(centroid_pca(:,1),centroid_pca(:,2),60,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1.5)
%   legend('Data','Centroids')
%   legend boxoff

%% Find representative points that have 1- highest density AND 2- contributions from the most runs
% Find the density by histogramming the data into 2D bins
figure
[values, centers] = hist3([p(:, 1) p(:, 2)],[51 51]); % bins may need modification
imagesc(centers{:},values.')
colorbar
axis xy
xlabel('PCA 1', 'fontsize', 25);
ylabel('PCA 2', 'fontsize', 25);
title("Density plot of first " + num2str(kPrinComp) + " " + options.SaveName + " PCs")

hold on
% Plot centroids with text and legend:
% Don't forget that you're a LEGEND:
legend_text = cell(size(centroid_pca,1),1);
  for i = 1:size(centroid_pca,1)
    scatter(centroid_pca(:,1),centroid_pca(:,2),60,'MarkerEdgeColor',[0 .5 .5],...
          'MarkerFaceColor',[0.8500 0.3250 0.0980],...
          'LineWidth',1.5)
     % Add text label for runs with proper offset
    text(centroid_pca(i,1)+pcaRange/50,centroid_pca(i,2)+pcaRange/50,['C' num2str(i)],'FontSize',14)
    legend_text{i} = ['C' num2str(i) ', size: ' num2str( pca_frame_centers(i,3))];
  end
legend(legend_text','location','EastOutside','fontsize', 16)
set(gcf, 'Position',  [100, 100, 700, 400]) % Resize window

if isfield(options, 'SavePath')
figPath = fullfile(options.SavePath, options.SaveName  + trjType + "_" +  "pcaLigandDensity");
savefig(figPath);
print2pdf(figPath);
end

%% Plot density as a landscape:


% 
%     xl = minmaxAll(p(:, 1));
%     yl = minmaxAll(p(:, 2));
% 
%     xi = xl(1):0.05:xl(2);
%     yi = yl(1):0.02:yl(2);
% 
%     gaussianWidth = [0.1 0.1];
%     pmf = calcpmf2d([p(:, 1), p(:, 2)], xi, yi, gaussianWidth);
%     % contourLevels = linspace(min(pmf, [], 'all'), max(pmf, [], 'all'), 10);
%     contourLevels = 0:0.5:3.5;
% figure
%     landscape(xi, yi, pmf, contourLevels);
% colorbar
% axis xy
% xlabel('PCA 1', 'fontsize', 25);
% ylabel('PCA 2', 'fontsize', 25);
% title(['Density plot of the first ' num2str(kPrinComp) ' PCs'])
% 
% hold on
% % Plot centroids with text and legend:
% % Don't forget that you're a LEGEND:
% legend_text = cell(size(centroid_pca,1),1);
%   for i = 1:size(centroid_pca,1)
%     scatter(centroid_pca(:,1),centroid_pca(:,2),60,'MarkerEdgeColor',[0 .5 .5],...
%           'MarkerFaceColor',[0.8500 0.3250 0.0980],...
%           'LineWidth',1.5)
%      % Add text label for runs with proper offset
%     text(centroid_pca(i,1)+pcaRange/50,centroid_pca(i,2)+pcaRange/50,['C' num2str(i)],'FontSize',14)
%     legend_text{i} = ['C' num2str(i) ', size: ' num2str( pca_frame_centers(i,3))];
%   end
% 
%   h = zeros(2,1);
%   h(1) = scatter(nan,nan,'MarkerEdgeColor',[0 .5 .5],...
%           'MarkerFaceColor',[0.8500 0.3250 0.0980],...
%           'LineWidth',1.5);
%   h(2) = scatter(nan,nan,'MarkerEdgeColor',[0 .5 .5],...
%           'MarkerFaceColor',[0.8500 0.3250 0.0980],...
%           'LineWidth',1.5);
% legend( h,legend_text','location','EastOutside','fontsize', 16)
% set(gcf, 'Position',  [100, 100, 700, 400]) % Resize window
% 
% 
% %% Functions
% 
% function output = minmaxAll(arr)
%   output = [min(arr, [], 'all'), max(arr, [], 'all')];
% end
end