% Load trajectories of WT and mutants, and then perform a common PCA to all
% of the common trajs

% Father directory for all the runs
metadir = '/scratch/1-D1_D2';
foldersToStudy = {'a-D1R_dpa_gp', 'b-D2R_WT'};
numRunsSys = [5 5]; % number of runs for each system

% foldersToStudy = {'3-d2_bromo_WT','4-d2_bromo_T174M-C220L','9-d2_bromo_V130F', ...
%     '11-d2_brc_L214M','13-d2_brc_L92G','15-d2_brc_F217M'};
% numRunsSys = [5 5 4 5 5 5]; % number of runs for each system

% thisSysLabel = 'abcdefghijklmnop'; % Can you sing the abc?
thisSysLabel = {'WT','T174M-C220L','V130F','L214M','L92G','F217M'};
%thisSysLabel = {'WT','T5.54M-C6.47L','L6.41M','L3.41G','F6.44M'};
frames2skip = 500; % remove first 50 ns of the simulation
nCentersLig = 10; % Number of points closest to the mean point of a system
% consider for analysis and pdb saving
name = 'dd2_dpaWT_F217M';
pdbName = 'prot.pdb';
%%

% Father directory for all the runs
metadir = 'D:\Telework_library\Chemokine_project\6-CXCR4_active_models_mutants';
foldersToStudy = {'a-WT_CXCR4_active_model','b-CXCR4_active_LSC_V3Y','c-CXCR4_active_Fusion_V3Y_Y7L','d-CXCR4_active_2ndGen_Y7L'};
numRunsSys = [7 5 5 5]; % number of runs for each system
% thisSysLabel = 'abcdefghijklmnop'; % Can you sing the abc?
thisSysLabel = {'WT','LSC','Fusion','2ndGen'};
frames2skip = 500; % remove first 50 ns of the simulation
nCentersLig = 10; % Number of points closest to the mean point of a system
% consider for analysis and pdb saving
name = 'CXCR4_active_models';
pdbName = 'prot.pdb';
%% Make the meta-direcory
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end
md2pathdirMeta = [metadir slash 'md2pathMeta' slash];
md2pathdir = [metadir slash 'md2pathMeta' slash name slash];

% Create meta directory
if ~exist(md2pathdirMeta, 'dir')
    mkdir(md2pathdirMeta)
    add2log(md2pathdirMeta, ['Creating md2path directory in ' metadir]);
end

% Create subdirectories for every bundle of systems studied
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir)
    add2log(md2pathdir, ['Creating md2path directory in ' md2pathdirMeta]);
end
%% Do PCA on more than 1 system

% Grab union of all ligand binding residues
thisList = [];
for thisSys = 1:length(foldersToStudy)
    L = importdata([metadir slash foldersToStudy{thisSys} slash md2pathName slash 'BS_residues.txt']);
    thisList = [thisList L];

end
liRes = unique(thisList); % Make a list of common ligand binding residues from all systems

tempProtTraj = [];
tempLigTraj = [];
RMSDCell =  cell(length(foldersToStudy),1);
pdbModel = cell(length(foldersToStudy),1);
crdModel = cell(length(foldersToStudy),1);
nFramesEff = cell(length(foldersToStudy),1); % I will use this to extract
% the pdb of the center
trajCeption = cell(length(thisSysLabel),1);

CSys = []; % Contains labels as: [ system run frameNdx ]
for thisSys = 1:length(foldersToStudy)

    mydir = ([metadir slash foldersToStudy{thisSys}]);
    numRuns = numRunsSys(thisSys);
    % load the trajetories:
    import_align_traj;

    % Load input model and align it: (if it exists)
    if exist([mydir slash 'step1_pdbreader.pdb'],'file')
        [pdbModel{thisSys}, temp] = readpdb([mydir slash 'step1_pdbreader.pdb']);
        [~, crdModel{thisSys}] = superimpose(crd, temp, find(index_CA_all));

    end

    index_CB_all = selectname(pdb.name, 'CB');
    index_CB = index_chainR & index_CB_all; % CBs of receptor only

    % Select ligand binding residues CA and CB only
    ligBindSele = selectid(pdb.resseq,liRes) & (index_CA_all | index_CB_all);

%     sele = index_CB | index_CA; % Do PCA on CA and CB
    sele = index_CA;
    if isLigSmall
        seleLig = index_L; % Non-H atoms
    else
        seleLig =  index_chainL & (index_CA_all | index_CB_all) ; % Same for ligand
    end
    % sele = index_chainR & noH;
    % sele = index_CA; % Selection for doing PCA
    % seleRMSD = index_chainR & noH;% Selection for doing RMSD calcs between centers


    nFramesEff{thisSys} = zeros(numRuns,1);
    RMSDCell{thisSys} = zeros(minTraj,numRuns);
%     startingFrame = 0;
    for runi = 1:numRuns
%        tempProtTraj = [tempProtTraj ; traj{runi}(frames2skip:end,to3(sele))]; % Make one fat trajectory
       tempLigTraj = [tempLigTraj ;  traj{runi}(frames2skip:end,to3(seleLig))]; % Make one fat trajectory
       % Calc RMSD along the way
       RMSDCell{thisSys}(:,runi) = calcrmsd(traj{runi}(1,:), traj{runi}(1:minTraj,:),ligBindSele);
%        if runi ~= 1 % Keep startingFrame = 0 for first run,
%          startingFrame = startingFrame + size(traj{runi-1},1);
%        end
       framesNdx = (frames2skip:size(traj{runi},1))';

       nFramesEff{thisSys}(runi) = size(traj{runi}(frames2skip:end,to3(sele)),1);
       Ctemp = [thisSys*ones(nFramesEff{thisSys}(runi),1)  ...
           runi*ones(nFramesEff{thisSys}(runi),1) framesNdx];
       CSys=[CSys ;Ctemp]; % Used for coloring and for labeling: [ system run frameNdx]
    end

    % Add input models if available:
    if exist([mydir '\step1_pdbreader.pdb'],'file')
%         tempProtTraj = [tempProtTraj ; crdModel{thisSys}(to3(sele))];
        tempLigTraj = [tempLigTraj ;  crdModel{thisSys}(to3(seleLig))];
        Ctemp = [thisSys 0 0]; % 0 for input models
        CSys=[CSys ;Ctemp]; % Used for coloring and for labeling: [ system run frameNdx]
    end
end

[modelNdx,~] = find(CSys==0);

%% Side quest: get RMSD of ligand binding residues:

% Plot RMSD
figure
errStrideRMSD = ceil(nFrames/60);
for thisSys = 1:length(foldersToStudy)

    RMSD_mean= mean(RMSDCell{thisSys},2);
    RMSD_std= std(RMSDCell{thisSys},0,2);


    plot(RMSD_mean,'LineWidth',1);
    hold on
    errorbar(1:errStrideRMSD:length(RMSD_mean),RMSD_mean(1:errStrideRMSD:length(RMSD_mean)), ...
        RMSD_std(1:errStrideRMSD:length(RMSD_mean))/sqrt(numRuns),'o','MarkerSize',5)
    hold on
    legend_entries{2*thisSys-1} = thisSysLabel{thisSys};
    legend_entries{2*thisSys} = '';
end

xlabel('Frame'); ylabel('RMSD [Angstrom]')
legend(legend_entries,'location','best');
legend boxoff
title('RMSD of ligand binding res, C\alpha and C\beta')
set(gca,'FontSize',16)

savefig([md2pathdir 'RMSD_ligand_binding']);
print2pdf([md2pathdir 'RMSD_ligand_binding']);
%% Do the actual PCA now:

% Protein

kClusters = [];
kmax = 15;
kPrinComp = 2; % Number of principal components to take into consideration

% [indexOfCluster_pcaProt, centroid_pcaProt, pProt, ind_centersProt] = cluster_traj_pca( tempProtTraj,kPrinComp, kClusters,[],kmax);
% 
% savefig([md2pathdir 'pcaProt_cluster_' name]);
% print2pdf([md2pathdir 'pcaProt_cluster_' name]);
% title(['Clustering with ' num2str(kPrinComp) ' receptor PCs']);
% 
%     clear tempProtTraj

% Ligand:

[indexOfCluster_pca, centroid_pca, p, ind_centers] = cluster_traj_pca( tempLigTraj,kPrinComp, kClusters,[],kmax);

title(['Clustering with ' num2str(kPrinComp) ' ligand PCs']);
savefig([md2pathdir 'pcaLig_cluster_' name]);
print2pdf([md2pathdir 'pcaLig_cluster_' name]);
%% Plot with runs:
pcaRange = range(pProt,'all');
pCentersProtLig = zeros(length(foldersToStudy),2);
pInputModelsProtLig = zeros(length(foldersToStudy),2);

figure
% subplot(3,1,1)
% scatter(pProt(:,1),pProt(:,2),5,CSys(:,1),'filled')
% hold on
% subplot(3,1,2)
% scatter(p(:,1),p(:,2),5,CSys(:,1),'filled')
% hold on
% subplot(3,1,3)
scatter(pProt(:,1),p(:,1),5,CSys(:,1),'filled')
hold on
  xlabel('Receptor PC 1', 'fontsize', 25);
  ylabel('Ligand PC 1', 'fontsize', 25);
%
for thisSys = 1:length(foldersToStudy)
    numRuns = numRunsSys(thisSys);

%     (find(CSys(:,1)==thisSys))

    p1 = mean(pProt(find(CSys(:,1)==thisSys),1));
    p2 = mean(p(find(CSys(:,1)==thisSys),1));
    pCentersProtLig(thisSys,:) = [p1 p2];
%     scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
%                   'MarkerFaceColor',[0 .7 .7],...
%                   'LineWidth',1.5)
    colorHere = colors(max(round((size(colors,1)/(length(thisSysLabel)-1)*(thisSys-1))),1),:);

    scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',colorHere,...
                  'LineWidth',1.5)
    hold on
    text(p1+pcaRange/50,p2+pcaRange/50,thisSysLabel{thisSys},'FontSize',14)

    % If input model is available, plot it!
    if ~isempty(crdModel{thisSys})
        p1 = pProt(modelNdx(thisSys),1);
        p2 = p(modelNdx(thisSys),1);
        scatter(p1,p2,60,'MarkerEdgeColor',[0.5 .5 0],...
                      'MarkerFaceColor',[0.7 .7 0],...
                      'LineWidth',1.5)
        hold on
        text(p1+pcaRange/50,p2+pcaRange/50,[thisSysLabel{thisSys} 'M'],'FontSize',14)
        pInputModelsProtLig(thisSys,:) = [p1 p2];
    end


%     for runi = 1:numRuns
%         % Find centroids of the scatter of each run and label it
%         if runi==1
%            lower_lim = 1;
%         else
%            lower_lim = sum(nFramesEff(1:runi-1))+1;
%         end
%         upper_lim = sum(nFramesEff(1:runi));
%         p1 = mean(pProt(lower_lim:upper_lim,1));
%         p2 = mean(pProt(lower_lim:upper_lim,2));
%
%         scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
%                       'MarkerFaceColor',[0 .7 .7],...
%                       'LineWidth',1.5)
%                   hold on
%         % Add text label for runs with proper offset
%         text(p1+pcaRange/50,p2+pcaRange/50,['Run ' num2str(runi)],'FontSize',14)
%     end
end

%   xlabel('PC 1', 'fontsize', 25);
%   ylabel('PC 2', 'fontsize', 25);
%   title(['Scatter plot of ' num2str(kPrinComp) ' receptor PCs, colored by run'])
    legend('Data colored by run','Centroids of runs','Input models')
  legend boxoff

  savefig([md2pathdir 'pcaProtLig_system_scatter_' name]);
print2pdf([md2pathdir 'pcaProtLig_system_scatter_' name]);
writematrix(pCentersProtLig,[md2pathdir 'pCentersProtLigCoordPC1.txt'],'Delimiter','space')

  %% Ligand only PCAs:
  pcaRange = range(p,'all');
  pCentersLig = zeros(length(foldersToStudy),2);
  pInputModelsLig = zeros(length(foldersToStudy),2);
  ind_centersRunsLig = zeros(length(foldersToStudy),nCentersLig);
  run_FrameNdx = cell(length(foldersToStudy),1);

  figure
colormapName = 'turbo';
colors = colormap(colormapName);

% Plot the PCs first
for thisSys = 1:length(foldersToStudy)
    colorHere = colors(max(round((size(colors,1)/(length(thisSysLabel)-1)*(thisSys-1))),1),:);
    scatter(p(CSys(:,1)==thisSys,1),p(CSys(:,1)==thisSys,2),5,colorHere,'filled');
    hold on
end

% Now plot centers and add text labels
for thisSys = 1:length(foldersToStudy)
    numRuns = numRunsSys(thisSys);

%     (find(CSys(:,1)==thisSys))

    p1 = mean(p(find(CSys(:,1)==thisSys),1));
    p2 = mean(p(find(CSys(:,1)==thisSys),2));
    pCentersLig(thisSys,:) = [p1 p2];

    colorHere = colors(max(round((size(colors,1)/(length(thisSysLabel)-1)*(thisSys-1))),1),:);

    scatter(p1,p2,60,'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',colorHere,...
                  'LineWidth',1.5)
    hold on
    text(p1+pcaRange/50,p2+pcaRange/50,thisSysLabel{thisSys},'FontSize',14)

    % If input model is available, plot it!
    if ~isempty(crdModel{thisSys})
        p1 = p(modelNdx(thisSys),1);
        p2 = p(modelNdx(thisSys),2);
        scatter(p1,p2,60,'MarkerEdgeColor',[0.5 .5 0],...
                      'MarkerFaceColor',[0.7 .7 0],...
                      'LineWidth',1.5)
        hold on
        text(p1+pcaRange/50,p2+pcaRange/50,[thisSysLabel{thisSys} 'M'],'FontSize',14)
        pInputModelsLig(thisSys,:) = [p1 p2];
    end

    % Save the pdbs of the cluster centers:
    % Only consider pcs from this system:
    pTemp = p(CSys(:,1)==thisSys,1:kPrinComp);
    % The closest point to the mean point
   [~,ind_centersRunsLig(thisSys)] = min(vecnorm(pTemp-pCentersLig(thisSys,:),2,2));

   % The X closest points to the mean point
   [B,I] = sort(vecnorm(pTemp-pCentersLig(thisSys,:),2,2));
   ind_centersRunsLig(thisSys,:) = I(1:nCentersLig);

   % ind_centersRunsLig is the index of the center of all the runs of a
   % system, the number is the frame number of given sys with @frames2skip
   % removed

    % Now use frame index to grab the real frame and extract it from traj:
    CSysTemp = CSys(CSys(:,1)==thisSys,:); % Csys for this system only

    run_FrameNdx{thisSys} = CSysTemp(ind_centersRunsLig(thisSys,:),2:3); % Got 'em!

end
  xlabel('PC 1', 'fontsize', 25);
  ylabel('PC 2', 'fontsize', 25);
%     legend('Data colored by run','Centroids of runs','Input models')
  legend(thisSysLabel)
  legend boxoff
  set(gca,'FontSize',16)

   title(['Scatter plot of ' num2str(kPrinComp) ' ligand PCs, colored by run'])

   % Save the plot
savefig([md2pathdir 'pcaLig_system_scatter_' name]);
print2pdf([md2pathdir 'pcaLig_system_scatter_' name]);
% Write the PCA values of the centers and the models
writematrix(pCentersLig,[md2pathdir 'pCentersLigCoordPC1_2.txt'],'Delimiter','space')
writematrix(p(CSys(:,3)==0,1:2),[md2pathdir 'pModelCentersLigCoordPC1_2.txt'],'Delimiter','space')

myLimits = axis; % Grab limits for density plots (next section)
%% Caclulate density plots of every system and save highest density cluster
% centers

  pHighestDen = zeros(length(foldersToStudy),2);
  ind_HighestDenLig = zeros(length(foldersToStudy),nCentersLig);
  run_FrameHighestDenNdx = cell(length(foldersToStudy),1);
  figure
  nPlots = ceil(sqrt(length(foldersToStudy)));

for thisSys = 1:length(foldersToStudy)
    subplot(nPlots,nPlots,thisSys);
     % Only consider pcs from this system:
    pTemp = p(CSys(:,1)==thisSys,1:kPrinComp);

    [values, centers] = hist3([pTemp(:, 1) pTemp(:, 2)],[51 51]); % bins may need modification
    imagesc(centers{:},values.')
    axis(myLimits);
    colorbar
    axis xy
    xlabel('PCA 1', 'fontsize', 16);
    ylabel('PCA 2', 'fontsize', 16);
    title( thisSysLabel{thisSys})
    sgtitle( ['Density plot of the first ' num2str(kPrinComp) ' PCs'])
    % Find highest density point/s

    % How to decide how many maxes to take?
    maxDen = max(values,[],'all');
    [a,b] = find(values==maxDen);
    pHighestDen(thisSys,:) = [centers{1}(a) centers{2}(b)]; % Highest density PC

    % Find closest PC point to our highest density center
   [~,ind_HighestDenLig(thisSys)] = min(vecnorm(pTemp-pHighestDen(thisSys,:),2,2));
   % The X closest points to the mean point
   [B,I] = sort(vecnorm(pTemp-pHighestDen(thisSys,:),2,2));
   ind_HighestDenLig(thisSys,:) = I(1:nCentersLig);

   CSysTemp = CSys(CSys(:,1)==thisSys,:); % Csys for this system only

   run_FrameHighestDenNdx{thisSys} = CSysTemp(ind_HighestDenLig(thisSys,:),2:3); % Got 'em!
    % extract structure et voila!

end
savefig([md2pathdir 'pcaLig_density_' name]);
print2pdf([md2pathdir 'pcaLig_density_' name]);
% Write the PCA values of the highest density points
writematrix(pHighestDen,[md2pathdir 'pHighestDenLigCoordPC1_2.txt'],'Delimiter','space')

%% Finally save the pdbs of the centers of the clusters of separate systems:

% CSys and PCs of the ligand cluster centers
CSysClusterCenters = [CSys(ind_centers,:)  (1:length(ind_centers))'];
pClusterCenters = p(ind_centers,1:kPrinComp);
pdbCell = cell(length(foldersToStudy),1);
indexLigCell = cell(length(foldersToStudy),1);

 pcaTraj = cell(length(foldersToStudy),1);
 pcaTrajHighestDen = cell(length(foldersToStudy),1);
 pcaClusterTraj = cell(length(ind_centers),1);
  for thisSys = 1:length(foldersToStudy)

    mydir = ([metadir slash foldersToStudy{thisSys}]);
    numRuns = numRunsSys(thisSys);
    % load the trajetories:
    import_align_traj;
    pdbCell{thisSys} = pdb;
    indexLigCell{thisSys} = index_chainL;
    % Save pdbs of centers of systems
    for modelN = 1:nCentersLig
        pcaTraj{thisSys} = [ pcaTraj{thisSys} ;traj{run_FrameNdx{thisSys}(modelN,1)}(run_FrameNdx{thisSys}(modelN,2),:)];
        pcaTrajHighestDen{thisSys} = [ pcaTrajHighestDen{thisSys} ;traj{run_FrameHighestDenNdx{thisSys}(modelN,1)}(run_FrameHighestDenNdx{thisSys}(modelN,2),:)];
    end
    writepdbndx([md2pathdir 'pcaLigand_' thisSysLabel{thisSys} '.pdb'], pdb, [], 0, pcaTraj{thisSys});

    % Save pdbs of highest density points of every system
    writepdbndx([md2pathdir 'pcaLigandHighestDen_' thisSysLabel{thisSys} '.pdb'], pdb, [], 0, pcaTrajHighestDen{thisSys});

    % Save pdbs of PCA cluster centers:

    for thisCenter = 1:length(ind_centers)
        if CSysClusterCenters(thisCenter,1) == thisSys % Cluster is in this system, save PDB!
            pcaClusterTraj{thisCenter} = traj{CSysClusterCenters(thisCenter,2)}(CSysClusterCenters(thisCenter,3),:);
            writepdbndx([md2pathdir 'pcaClusterLigand_' thisSysLabel{thisSys} '_C' num2str(thisCenter) '.pdb'], pdb, [], 0, pcaClusterTraj{thisCenter});
        end
    end

  end


  %% Calculate RMSD matrix between centers: (Rob specific)


for k = 1:2 % First run for system centers, second for
    indexLigCell = cell(length(foldersToStudy),1);
    commonpcaTraj = [];

    for thisSys = 1:length(foldersToStudy) % Prepare the index for all the common atoms of the ligand
         lig37Common =  selectid(pdbCell{thisSys}.resseq,ligandRes([3 7])) & selectname(pdbCell{thisSys}.name, 'CA', 'C', 'N', 'O','CB','CG','CG1');
         indexLigCell{thisSys} = selectid(pdbCell{thisSys}.resseq,ligandRes([1:2 4:6 8:end])) | lig37Common;
         % Make the trajectory of common atoms only:
         if k==1
            commonpcaTraj = [commonpcaTraj ; pcaTraj{thisSys}(:,to3(indexLigCell{thisSys}))];
         elseif k==2
            commonpcaTraj = [commonpcaTraj ; pcaTrajHighestDen{thisSys}(:,to3(indexLigCell{thisSys}))];
         end
    end

    % Add input models if they exist:
    for thisSys = 1:length(foldersToStudy)
        if ~isempty(crdModel{thisSys})
            commonpcaTraj = [commonpcaTraj ; crdModel{thisSys}(:,to3(indexLigCell{thisSys}))];
        end
    end
    %

    rmsdLigandpcaCenter = zeros(size(commonpcaTraj,1));
    for i =1:size(commonpcaTraj,1)-1
        i
        for j = i+1:size(commonpcaTraj,1)
            rmsdLigandpcaCenter(i,j) = calcrmsd(commonpcaTraj(i,:),commonpcaTraj(j,:));
        end
    end

    figure
    hm = heatmap(rmsdLigandpcaCenter,'Colormap',colors);
    hm.Title = ['Ligand all common-atom RMSD [' char(197) '] matrix PCA top ' num2str(nCentersLig)];
    hm.XLabel = 'PCA Cluster';
    hm.YLabel = 'PCA Cluster';

    if k==1
       figName = 'pcaSysCenterrmsdMatrix';
    elseif k==2
       figName = 'pcaHighestDenrmsdMatrix';
    end
    savefig([md2pathdir figName]);
    print2pdf([md2pathdir figName]);

    % Write out the matrix to an excel file:

    % Write labels first:
    tableLabels = cell(size(commonpcaTraj,1),1);
    for thisCenter =1:nCentersLig
        for thisSys = 1:length(foldersToStudy)
            tableLabels{thisCenter + nCentersLig*(thisSys-1)} = [thisSysLabel{thisSys} '-' num2str(thisCenter)];
        end
    end
    for thisSys = 1:length(foldersToStudy)
        if ~isempty(crdModel{thisSys})
            tableLabels{nCentersLig*length(foldersToStudy) + thisSys} = [thisSysLabel{thisSys} '-Input']  ;
        end
    end

    if k==1
       xlsName = [md2pathdir 'pcaModelsRMSD_' name '.xls',];
    elseif k==2
       xlsName = [md2pathdir 'pcaHighestDenRMSD_' name '.xls',];
    end

    if exist(xlsName,'file')
        add2log(md2pathdir, ['Warning!! ' xlsName ' already exists! Will append to it' ]);
    end

    writecell(tableLabels',xlsName,'Range','B1')
    writecell(tableLabels,xlsName,'Range','A2')
    writematrix(rmsdLigandpcaCenter,xlsName,'Range','B2')
end
%% Calcualte dihedrals of ligand (if peptide) and ligand binding residues of
% top models and input models

% Elements of the trajectory:
% pdbCell, crdModel, pcaTraj
higherOrder = 'all';

if isLigSmall
    dihRes = liRes;
else
    dihRes = [ liRes' ; ligandRes];
end

dihedrals =[];
dihedrals_mat = cell(length(foldersToStudy),1);
reSortCell = cell(length(foldersToStudy),1);
for thisSys = 1:length(foldersToStudy)
    % Add input models if they exist:
    if ~isempty(crdModel{thisSys})
        tempTraj = [pcaTraj{thisSys} ; crdModel{thisSys}];
    else
        tempTraj = pcaTraj{thisSys};
    end

% calculate dihedrals for ligand and ligand binding residues
    [dihedralsTemp, dihIndex,reSortCell{thisSys}] = calcalldihedralsfromtrajs(pdbCell{thisSys} ...
        ,tempTraj,dihRes,1,higherOrder);
%     dihedrals = [ dihedrals dihedralsTemp];
    dihedrals = dihedralsTemp;

    dihedrals_mat{thisSys} = [];
    for dih = 1:length(dihedrals)
        X = dihedrals{dih};
        if dih == 1 && isnan(X(1,1))
            X(:,1) = []; % Remove the NaNs from Nterm column 1
        elseif dih == length(dihedrals) && isnan(X(1,2))
            X(:,2) = []; % Remove the NaNs from Cterm column 2
        end
        dihedrals_mat{thisSys} = [dihedrals_mat{thisSys} X];
    end
    dihedrals_mat{thisSys} = real(dihedrals_mat{thisSys});
end

% Write out the dihedrals to an excel file:
% we want to write them like this:
% [dihRes(reSortCell{thisSys}(:,1)) reSortCell{thisSys}(:,2) dihedrals_mat{thisSys}']

% Write labels first:
dihedralLabels = {'phi','psi','chi1','chi2','chi3','chi4','chi5'};
xlsName = [md2pathdir 'pcaModelsDihedrals_' name '.xls',];
if exist(xlsName,'file')
    add2log(md2pathdir, ['Warning!! ' xlsName ' already exists! Will overwrite/append to it' ]);
end



for thisSys = 1:length(foldersToStudy)

    tableDihLabels = cell(nCentersLig+1,1);
    for thisCenter =1:nCentersLig
        tableDihLabels{thisCenter} = [thisSysLabel{thisSys} '-' num2str(thisCenter)];
    end
    if ~isempty(crdModel{thisSys})
        tableDihLabels{nCentersLig+1} = [thisSysLabel{thisSys} '-Input']  ;
    end


    writecell(tableDihLabels',xlsName,'Sheet', thisSysLabel{thisSys},'Range','C1')
    writematrix([dihRes(reSortCell{thisSys}(:,1)) reSortCell{thisSys}(:,2)],xlsName,'Sheet', thisSysLabel{thisSys},'Range','A2')
    writematrix(dihedrals_mat{thisSys}',xlsName,'Sheet', thisSysLabel{thisSys},'Range','C2')
    writecell({'Residue','Dih Type'},xlsName,'Sheet', thisSysLabel{thisSys},'Range','A1')

end

%% Modernized meta PCA, fit for AlloDy 1.0:

metaLoadAlign;
colormapName = 'parula';
colors = colormap(colormapName);

% This is for receptor chain: ONLY CA PCA FOR NOW
if attempt2AlignTrajs
    simResList =  database.residues{1}{:,1:length(foldersToStudy)}; % list from simulation set
    allAligned = sum(simResList==0,2)==0;
    alignmentVector = simResList(allAligned,:); % Attempt to line different systems
    % alignmentVector works if:
    % 1- residue numbers start from 1
    % 2- residue numbers of same chain are sorted in increasing order (no funny 
    % PDB files where residue 500 is before residue 1 in the same chain)
end

tempProtTraj = [];
tempLigTraj = [];
CSys = []; % Contains labels as: [ system run frameNdx ]
nFramesEff = cell(length(foldersToStudy),1); % I will use this to extract
% the pdb of the center

% Read and concatenate trajs
for thisSys = 1:length(foldersToStudy)

    myDir = fullfile(metadir,foldersToStudy{thisSys});
    entryHere  = database.entries{thisSys};
    chainHere =  entryHere.chains{Chains.receptor};
    simHere = entryHere.addSimulation(fullfile(myDir, "run*"),'align2chain',chains(Chains.receptor));
    
    % Protein stuff
    protAtomIndices = chainHere.getAtoms(); % Grab CA atoms
    traj = simHere.concatRuns('Atoms', protAtomIndices, 'StartFrame', frames2skip + 1);
    if attempt2AlignTrajs
        tempProtTraj = [tempProtTraj ;traj(:,to3(alignmentVector(:,thisSys)))]; % Matlab will complain about this :D :D 
    else
        tempProtTraj = [tempProtTraj ;traj];
    end

    % Ligand stuff
    ligAtomIndices = entryHere.chains{Chains.ligand}.getLigandAtoms; % Assumes same ligand among different simulations
    trajLig = simHere.concatRuns('Atoms', ligAtomIndices, 'StartFrame', frames2skip + 1);
    tempLigTraj = [tempLigTraj ;trajLig];
    % Grab the index Csys
    numRuns = simHere.runCount;
    for runi = 1:numRuns
        framesNdx = ((frames2skip + 1):size(simHere.traj{runi},1))';
        nFramesEff{thisSys}(runi) = size(simHere.traj{runi}((frames2skip + 1):end,to3(protAtomIndices)),1);
        Ctemp = [thisSys*ones(nFramesEff{thisSys}(runi),1)  ...
           runi*ones(nFramesEff{thisSys}(runi),1) framesNdx];
        CSys=[CSys ;Ctemp]; % Used for coloring and for labeling: [ system run frameNdx]
    end
end

% Do the actual PCA now:

%% Protein

kClusters = [];
kmax = 15;
kPrinComp = 2; % Number of principal components to take into consideration

[indexOfCluster_pcaProt, centroid_pcaProt, pProt, ind_centersProt] = cluster_traj_pca( tempProtTraj,kPrinComp, kClusters,[],kmax);

savefig([md2pathdir 'pcaProt_cluster_' name]);
print2pdf([md2pathdir 'pcaProt_cluster_' name]);
title(['Clustering with ' num2str(kPrinComp) ' receptor PCs']);

clear tempProtTraj traj

%% Scatter PCA (protein) colored by system:
[pCenters, run_FrameNdx] = plotPCAScatter(pProt, kPrinComp, CSys, 'pathName', ...
    [md2pathdir 'pcaProt_system_' name], 'thisSysLabel', thisSysLabel);


%% Density of PCA space:

[pHighestDen, run_FrameHighestDenNdx] = plotPCADensity(pProt, kPrinComp, CSys, 'pathName', ...
    [md2pathdir 'pcaProt_density_' name], 'thisSysLabel', thisSysLabel);

%% Finally save the pdbs of the centers of the clusters of separate systems:

% What to input here?
% ind_centers, run_FrameNdx, run_FrameHighestDenNdx, nCenters
% file name: pcaLigand or pcaProt or pcaProtLigand

pcaName = 'pcaProt';
ind_centersHere = ind_centersProt;
run_FrameNdxHere = run_FrameNdx;
run_FrameHighestDenNdxHere = run_FrameHighestDenNdx;

metapca_savePDBs;


%% Cluster every system separately in common PC space?
figure
tiledlayout('flow')
for thisSys = 1:length(foldersToStudy)
    nexttile
    pTemp = pProt(CSys(:,1)==thisSys,1:kPrinComp);
    [indexOfCluster_pca_sorted, centroid_pca_sorted, ind_centers] = cluster_traj_only(pTemp, kPrinComp, kClusters, kmax);
    title( thisSysLabel{thisSys})
end
sgtitle( ['Scatter of first ' num2str(kPrinComp) ' PCs'])

%% Calculate RMSD matrix between centers



%% Ligand 
[indexOfCluster_pca, centroid_pca, p, ind_centers] = cluster_traj_pca( tempLigTraj,kPrinComp, kClusters,[],kmax);

savefig([md2pathdir 'pca_cluster_' name '.fig']);
print2pdf([md2pathdir 'pca_cluster_' name]);
title(['Clustering with ' num2str(kPrinComp) ' receptor PCs']);

% clear tempLigTraj trajLig

%% Scatter PCA colored by system:
[pCenters, run_FrameNdx] = plotPCAScatter(p, kPrinComp, CSys, 'pathName', ...
    [md2pathdir 'pca_system_' name], 'thisSysLabel', thisSysLabel,'colormapName','turbo');


% Density of PCA space:
[pHighestDen, run_FrameHighestDenNdx] = plotPCADensity(p, kPrinComp, CSys, 'pathName', ...
    [md2pathdir 'pca_density_' name], 'thisSysLabel', thisSysLabel);

% Save PDBs of centers:
pcaName = 'pca';
ind_centersHere = ind_centers;
run_FrameNdxHere = run_FrameNdx;
run_FrameHighestDenNdxHere = run_FrameHighestDenNdx;

metapca_savePDBs;

%% Calc RMSD between highest density centers for different systems:
% Elements that we will use:
% pcaTrajHighestDen
% pcaTraj

rmsdHighestDenMean = zeros(length(foldersToStudy));
rmsdLigandpcaCenter = zeros(length(nCenters)*length(foldersToStudy));


for thisSysi = 1:length(foldersToStudy)
    myDir = fullfile(metadir,foldersToStudy{thisSysi});
    entryi  = database.entries{thisSysi};
    chaini =  entryi.chains{Chains.receptor};
    
    protAtomIndicesi = chaini.getAtoms(); % Grab CA atoms
    ligAtomIndicesi = entryi.chains{Chains.ligand}.getLigandAtoms; % Assumes same ligand among different simulations
    for thisSysj = thisSysi:length(foldersToStudy)
    entryj  = database.entries{thisSysj};
    chainj =  entryj.chains{Chains.receptor};
    
    protAtomIndicesj = chainj.getAtoms(); % Grab CA atoms
    ligAtomIndicesj = entryj.chains{Chains.ligand}.getLigandAtoms; % Assumes same ligand among different simulations
    temp = zeros(nCenters);
    for i =1:nCenters%-1
        ndxi = nCenters*(thisSysi-1) + i;
        for j = 1:nCenters
            ndxj = nCenters*(thisSysj-1) + j;
            rmsdLigandpcaCenter(ndxi,ndxj) = calcrmsd(pcaTrajHighestDen{thisSysi}(i,to3(ligAtomIndicesi)), ...
                pcaTrajHighestDen{thisSysj}(j,to3(ligAtomIndicesj)));
            temp(i,j) = rmsdLigandpcaCenter(ndxi,ndxj) ;
            
        end
    end
    rmsdHighestDenMean(thisSysi,thisSysj) = mean(nonzeros(temp),'all');
    end
end

figure
hm = heatmap(rmsdLigandpcaCenter,'Colormap',parula);
hm.Title = ['Ligand RMSD [' char(197) '] matrix PCA highest density '];
hm.XLabel = 'PCA Cluster';
hm.YLabel = 'PCA Cluster';
 figPath = fullfile(md2pathdir, "pcaLigandrmsdMatrixHighestDen");
savefig(figPath);
print2pdf(figPath);

figure
hm2 = heatmap(thisSysLabel,thisSysLabel, rmsdHighestDenMean,'Colormap',turbo);
hm2.Title = ['Ligand RMSD [' char(197) '] matrix highest density centers '];

% 
 figPath = fullfile(md2pathdir, "pcaLigandrmsdMatrixHighestDenMean");
savefig(figPath);
print2pdf(figPath);
%% Support functions

function [pCenters, run_FrameNdx] = plotPCAScatter(p, kPrinComp, CSys, options)
arguments
    p
    kPrinComp
    CSys
    options.colormapName = 'parula';
    options.pathName = []; % For saving 
    options.thisSysLabel = cell(max(CSys(:,1)),1)
    options.nCenters = 10 % How many points close to the center to output
end
    nSys = max(CSys(:,1)); % Number of systems to study
    if nargout > 1
        ind_centersRuns = zeros(nSys,options.nCenters);
        run_FrameNdx = cell(nSys,1);
    end

    pcaRange = range(p,'all');
    pCenters = zeros( nSys,kPrinComp);
%     pInputModelsProt = zeros( max(CSys(:,1)),kPrinComp);
    
    figure
    scatter(p(:,1),p(:,2),5,CSys(:,1),'filled')
    colors = colormap(options.colormapName);
    colormap(options.colormapName);
    hold on
      xlabel('PC 1', 'fontsize', 25);
      ylabel('PC 2', 'fontsize', 25);
    
    for thisSys = 1: nSys
%         numRuns =  max(CSys(CSys(:,1)==thisSys,2));
    
        pCenters(thisSys,:) = mean(p(find(CSys(:,1)==thisSys),1:kPrinComp));
        p1 = pCenters(thisSys,1);
        p2 =  pCenters(thisSys,2);
    
        colorHere = colors(max(round((size(colors,1)/(length(options.thisSysLabel)-1)*(thisSys-1))),1),:);
    
        scatter(p1,p2,60,'MarkerEdgeColor',[0.8 0.8 0.8],... %[0 .5 .5]
                      'MarkerFaceColor',colorHere,...
                      'LineWidth',1.5)
        hold on
        if ~isempty(options.thisSysLabel{thisSys})
            text(p1+pcaRange/50,p2+pcaRange/50,options.thisSysLabel{thisSys},'FontSize',14)
        end
       
        if nargout > 1
           % Save the pdbs of the cluster centers:
           % Only consider pcs from this system:
           pTemp = p(CSys(:,1)==thisSys,1:kPrinComp);
           % The closest point to the mean point
           [~,ind_centersRuns(thisSys)] = min(vecnorm(pTemp-pCenters(thisSys,:),2,2));
        
           % The X closest points to the mean point
           [B,I] = sort(vecnorm(pTemp-pCenters(thisSys,:),2,2));
           ind_centersRuns(thisSys,:) = I(1:options.nCenters);
        
           % ind_centersRuns is the index of the center of all the runs of a
           % system, the number is the frame number of given sys with @frames2skip
           % removed
        
           % Now use frame index to grab the real frame and extract it from traj:
           CSysTemp = CSys(CSys(:,1)==thisSys,:); % Csys for this system only
        
           run_FrameNdx{thisSys} = CSysTemp(ind_centersRuns(thisSys,:),2:3); % Got 'em
        end
    end
    
      title(['Scatter plot of ' num2str(kPrinComp) ' receptor PCs, colored by system'])
        legend('Data colored by run','Centroids of runs','Input models')
      legend boxoff
      
      if ~isempty(options.pathName)
        savefig([options.pathName '.fig']);
        print2pdf(options.pathName);
        writematrix(pCenters,[ options.pathName '.txt'],'Delimiter','space')
      end
end


function [pHighestDen, run_FrameHighestDenNdx] = plotPCADensity(p, kPrinComp, CSys, options)
arguments
    p
    kPrinComp
    CSys
    options.colormapName = 'parula';
    options.pathName = [];
    options.thisSysLabel = cell(max(CSys(:,1)),1)
    options.nCenters = 10 % How many points close to the center to output
end
      nSys = max(CSys(:,1)); % Number of systems to study 
      pHighestDen = zeros(nSys,2);
      ind_HighestDen = zeros(nSys,options.nCenters);
      run_FrameHighestDenNdx = cell(nSys,1);
      figure
      tiledlayout('flow')
    for thisSys = 1:nSys
        nexttile
         % Only consider pcs from this system:
        pTemp = p(CSys(:,1)==thisSys,1:kPrinComp);
    
        [values, centers] = hist3([pTemp(:, 1) pTemp(:, 2)],[51 51]); % bins may need modification
        imagesc(centers{:},values.')
        colormap(options.colormapName);
    %     axis(myLimits);
        colorbar
        axis xy
        xlabel('PCA 1', 'fontsize', 16);
        ylabel('PCA 2', 'fontsize', 16);
        title( options.thisSysLabel{thisSys})
        sgtitle( ['Density plot of the first ' num2str(kPrinComp) ' PCs'])
        % Find highest density point/s
    
        % How to decide how many maxes to take?
        maxDen = max(values,[],'all');
        [a,b] = find(values==maxDen);
        pHighestDen(thisSys,:) = [centers{1}(a(1)) centers{2}(b(1))]; % Highest density PC
        hold on
        scatter(pHighestDen(thisSys,1),pHighestDen(thisSys,2),50,"k", 'marker','x', 'LineWidth',1.5)
        % Find closest PC point to our highest density center
       [~,ind_HighestDen(thisSys)] = min(vecnorm(pTemp-pHighestDen(thisSys,:),2,2));
       % The X closest points to the mean point
       [B,I] = sort(vecnorm(pTemp-pHighestDen(thisSys,:),2,2));
       ind_HighestDen(thisSys,:) = I(1:options.nCenters);
    
       CSysTemp = CSys(CSys(:,1)==thisSys,:); % Csys for this system only
    
       run_FrameHighestDenNdx{thisSys} = CSysTemp(ind_HighestDen(thisSys,:),2:3); % Got 'em!
        % extract structure et voila!
    end
    if ~isempty(options.pathName)
        savefig([ options.pathName '.fig']);
        print2pdf(options.pathName);
        % Write the PCA values of the highest density points
        writematrix(pHighestDen,[ options.pathName '.txt'],'Delimiter','space')
    end
end
