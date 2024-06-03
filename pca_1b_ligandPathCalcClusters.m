%% What to do with clusters?
% Extract trajectories of every cluster
% Find ligand binding residues
% We need the cluster indices for this section: clusterNdx_filtered_all

% Inputs
%  - C
%  - clusterNdx_filtered_all
%  - receptorGpResIds
%  - pcadir

% Options for subruns: (Now they're in md2pathmain)
%         % Extract cluster trajectories:
%         c2Analyze = [1 2]; % clusters to analyze
%         c2Analyze(c2Analyze > max(clusterNdx_filtered_all)) = []; % Remove clusters that
%         % don't exist to avoid errors later
%         
%         importantCutOff = 0.7; % use more stringent cutoff here than md2pathmain
%         % Dual cutoff contact scheme: two atoms come into contact if within
%         % rcut1 but only disconnect when they're further than rcut2
%         rcut1 = 3.5;
%         rcut2 = 5;
%         pcaContactCut = 0.1; % Would be different than normal contact cut where the
%         % cutoff is "spread" over different runs
%% Create one simulation object for each cluster

% Don't forget to have the dihedrals loaded in mainSim!

clusterSims = arrayfun(@(Ci) mainSim.createSubset(C(clusterNdx_filtered_all == Ci, :), 'DihedralMatIndices', clusterNdx_filtered_all == Ci), c2Analyze);
% clusterSims = arrayfun(@(Ci) mainSim.createSubset(C(indexOfCluster_pca_filtered == Ci, :), 'DihedralMatIndices', indexOfCluster_pca_filtered == Ci), c2Analyze);

%% Calculate ligand contact map

contact_ligCluster = cell(length(c2Analyze), 1);

for Ci = 1:length(c2Analyze)
    [contactsPerRun, ligandAtomIndices] = clusterSims(Ci).computeContacts(mainChain, ligandChain, 'ContactCut', pcaContactCut,'RCut1', rcut1, 'RCut2', rcut2);
    contact_ligCluster{Ci} = contactsPerRun(:, :, 1);
end


%% Plot contact map containing all non-zero elements

% Excel file to save ligand contact data
xlsName = fullfile(pcadir, "ligandContactDataClusters_" + mainEntry.name + ".xls");

if exist(xlsName, 'file')
    add2log(md2pathdir, "Warning!! " + xlsName + " already exists! Will append to it");
end

receptorLigandResIdsPerCluster = cell(length(c2Analyze), 1);

for Ci = 1:length(c2Analyze)
    clusterNumber = c2Analyze(Ci);

    receptorLigandResIdsPerCluster{Ci} = plot_save_contact_map(contact_ligCluster{Ci}, mainChain, ligandChain, ligandAtomIndices, ...
        'ImportantCutoff', importantCutOff, 'RCut', rcut2, ...
        'SaveDir', pcadir, ...
        'SaveName', ("%s_C" + clusterNumber), ...
        'SaveXlsName', "%s", ...
        'SaveXlsSheet', "Cluster C" + clusterNumber);
end



%% Calculate MI!

for Ci = 1:length(c2Analyze)
    miPath = fullfile(pcadir, sprintf("MI_C%d.mat", c2Analyze(Ci)));
    clusterSims(Ci).computeMI('Path', miPath);
end


%% Assess convergence of the entropies

for Ci = 1:length(c2Analyze)
    clusterNumber = c2Analyze(Ci);

    assessEntropyConvergence(clusterSims(Ci), ...
        'Name', mainEntry.name + ", cluster C" + clusterNumber, ...
        'SavePath', fullfile(pcadir, ['Convergence_C'  num2str(clusterNumber)]));
end


%% Now that we have contact maps for specific ligand poses, run path calculation!

for Ci = 1:length(c2Analyze)
    clusterNumber = c2Analyze(Ci);

    pathCalcdir = prepareAlloPathCalc( ...
        clusterSims(Ci), ...
        mainChain, ...
        settings, ...
        'Dir', pcadir, ...
        'Name', ("%s_C" + clusterNumber), ...
        'LogPath', md2pathdir, ...
        'ReceptorLigandResIds', receptorLigandResIdsPerCluster{Ci}, ...
        'ReceptorGpResIds', receptorGpResIds, ...
        'ReceptorResIds', receptorResIds, ...
        'Traj', highestDenTraj(clusterNumber, :) ...
    );

    load(fullfile(pathCalcdir,"workspace.mat"))

    prepMI;
    MIStatsResLevel;
    graphanalysis;
    ClusterMIpathways;
    analyzeClusters;
    writeChannels;
    analyzePathDomains;
    MIanalysisBS2Effector;
    pathwayBiasAnalysis;
    
    [pHere] = visualizeClsGraph(PDB,pathstruc,Gmatmajor,1,'MIFractionCutoff',MIFractionCutoff);
    % If protein is inversed along z, use this command to make it upright again
    % set(gca,'zdir','reverse')
    
    % If the "chirality" is off, you can use this:
    % set(gca,'ydir','reverse')
end


%% Add an analysis layer similar to metaPathCalcAnalysis
