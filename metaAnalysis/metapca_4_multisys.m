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

%% Modernized meta PCA

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
