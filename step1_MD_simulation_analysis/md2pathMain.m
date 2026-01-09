% Features to add:
% 
% Inherent double check for pathway calculation (calculate pathways again at
% 75% of simulation length?), in addition to convergence criteria
%
% Protein-protein contact map, in addition to prot-lig and prot-effector
% Protein pose clustering (RMSD?) + extraction of representative protein
% pose from every cluster. This is already done with PCA, but it leaves a
% lot to be desired
% Ligand binding residues to RMSF plot
% Ligand atoms instead of atom numbers on RMSF ligand plot
%% This set of scripts/functions aims to take the output of an MD simulation
% and calculate allosteric pathways, going thru several steps of analysis

% Preparations needed: PDB file with separate chains for receptor, ligand,
% and G-protein: chains can be set up in the input_md2path script
% NOTE: all input residues should be in the same pdb numbering as the input
% PDB
% DO NOT FORGET TO SETUP AND RUN input_md2path before running this script

% Small roadmap:
% 1- transform xtc to dcd
% 2- load trajs
% 3- Calc RMSD and RMSF
% 4- Calc contact map with ligand and/or Gp (if they exist)
% 5- Make a PCA of ligand binding poses and receptor coordinates 
% 6- Calc GPCR order param for activation/state
% 7- Calc dihedrals
% 8- Calc entropies
% 9- Evalutate convergence of entropies
% 10- Prep input for pathway calculation
% 11- Run pathway calculation

global settings
%% Check that the data directory exists

assert(logical(exist(settings.mydir, 'dir')), "Directory does not exist!!! Please point to an existing directory!");

% Initialize variables:
if ~isfield(settings,'pdbCodeExtra') 
    settings.pdbCodeExtra='';
end
if ~isfield(settings,'includeLigPathwayCalc') 
    settings.includeLigPathwayCalc=false;
end
if ~isfield(settings,'diagnosticsOn') 
    settings.diagnosticsOn=false;
end
if ~isfield(settings,'MIWeightPaths') 
    settings.MIWeightPaths=true;
end
if settings.doPCA && ~isfield(settings,'trjType') 
    trjType = 'CA';
elseif settings.doPCA
    trjType = settings.trjType;
    % Special case: distances
    if strcmp(trjType,'Distances') && ~isfield(settings,'NdisPCA') 
        NdisPCA = 1000; % Take highest variance distances for protein PCA
    end
    if strcmp(trjType,'Distances') && ~isfield(settings,'CGPCA')
        CGPCA = 4; % Coarse graining by taking every Nth residue for distance calc
    end
end
if settings.doPCA && ~isfield(settings,'trjTypeLigand') 
    trjTypeLigand = 'Default'; % CA for peptide ligand, heavy atoms for small molecules
elseif settings.doPCA
    trjTypeLigand = settings.trjTypeLigand;
end

%% Load, fetch and align PDB files

database = Database(settings.databasePath);

% Chains: [receptor, G protein, ligand]
disp(settings.mydir)
database.read(fullfile(settings.mydir, "prot.pdb"), settings.chains, "Protein");

if ~isempty(settings.pdbCode)
    database.fetch(settings.pdbCode, settings.pdbChains, 1);

    if ~isempty(settings.pdbCodeInactive)
        database.fetch(settings.pdbCodeInactive, settings.pdbInactiveChains, 2);
    end
    for i = 1:length(settings.pdbCodeExtra)
        database.fetch(settings.pdbCodeExtra{i}, settings.pdbChainsExtra{i});
    end
end

% Align sequences and produce the database.residues table
database.align();

% Align structures to the first database entry
database.alignStructures();

% Store names for each chain
Chains.receptor = 1;
Chains.gprotein = 2;
Chains.ligand = 3;

% Keep references to entries and chains used frequently
mainEntry = database.entries{1};
mainChain = mainEntry.chains{Chains.receptor};

if mainEntry.hasChain(Chains.ligand)
    ligandChain = mainEntry.chains{Chains.ligand};
end
if mainEntry.hasChain(Chains.gprotein)
    effectorChain = mainEntry.chains{Chains.gprotein};
end

% Temp section until Mathworks decides to fix how copyfile works on Linux
% systems
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end

%% Initialize the directory where output and figures will be saved

% TODO: Remove? Only used once

[~,name,~] = fileparts(settings.mydir); % gives path and current folder name
name(name=='_')='-';

if length(name) > 3 && name(2) == '-' % just my directory naming convention
    name = name(3:end);
elseif length(name) > 3 && name(3) == '-'
    name = name(4:end);
end

mainEntry.name = name;


md2pathdir = fullfile(settings.mydir, "md2pathdev");
saveVarName = fullfile(md2pathdir, "workspace.mat"); % Where to save important variables

% Create md2pathdir if it does not exist
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir);
    add2log(md2pathdir, "Creating md2path directory in " + settings.mydir);
end
if ~exist(saveVarName,'file')
    save(saveVarName,'md2pathdir');
else
    save(saveVarName,'md2pathdir','-append');
end

% Possibly copy some of the PDB to md2pathdir
% e.g. copyfile(database.entries{1}.path, fullfile(md2pathdir, "foo.pdb"));


%% Add labels

% Add labels with Ballesteros-Weinstein notation, aligned with the 2nd
% database entry
% Exported from GPCRdb (https://gpcrdb.org/residue/residuetable)
if length(database.entries) > 1 && isfield(settings, 'systemName')
    residueTablePath = fullfile(database.dir, settings.systemName + "_residue_table.xlsx");
    database.label(2, residueTablePath);
end

%% Find mutations compared to a downloaded reference sequence:

if exist(fullfile(settings.databasePath, settings.systemName + ".fasta"),'file') % Quick and dirty if statement, change later
    fastaPath = fullfile(settings.databasePath, settings.systemName + ".fasta");
    mutations = database.findMut(1, Chains.receptor, fastaPath);
    mutations
end
%% Transform .xtc files into .dcd

% Are the DCD files there? If yes, skip this step
% areThereDCDs = dir([settings.mydir '/run*/traj.dcd']);
areThereDCDs = dir(fullfile(settings.mydir, "run*", "traj.dcd"));

if length(areThereDCDs) < settings.numRuns
    % run VMD from command line:  vmd -dispdev text -e
    [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
    pathToScript = fullfile(scriptPath, "load_save.tcl");  % assumes script is in the same directory
    cmdStr       = "vmd -dispdev text -e " + pathToScript + " -args " + settings.mydir + " " + num2str(settings.numRuns) + " " + num2str(settings.stride) + " " + database.entries{1}.path + " " + settings.xtcName;
    system(cmdStr);
    add2log(md2pathdir, "Transformed .xtc files to .dcd files with the following command: """ + cmdStr + """");
end


%% Import trajectories, decenter and align them using CAs to first frame

mainSim = mainEntry.addSimulation(fullfile(settings.mydir, "run*"),'align2chain',settings.chains(Chains.receptor));


%% Calculate RMSD and RMSF of CAs and ligand


% Calculate RMSD and RMSF of the CA atoms over the course of the sims

% Receptor:
receptorAtoms = mainChain.getAtoms();
RMSD = mainSim.computeRmsd(receptorAtoms);
RMSF = mainSim.computeRmsf(receptorAtoms, 'StartFrame', settings.frames2skip + 1);
save(saveVarName,'RMSD','RMSF','-append');

bfactor = zeros(1, mainEntry.atomCount);
bfactor(1, receptorAtoms) = mean(RMSF, 2);

% Ligand
if mainEntry.hasChain(Chains.ligand)
    ligandAtoms = ligandChain.getLigandAtoms();
    RMSD_L = mainSim.computeRmsd(ligandAtoms);
    RMSF_L = mainSim.computeRmsf(ligandAtoms, 'StartFrame', settings.frames2skip + 1);
    bfactor(1, ligandAtoms) = mean(RMSF_L, 2);
    save(saveVarName,'RMSD_L','RMSF_L','-append');
end

% Effector protein: double check structural aligment of traj, if it's
% centered at the reeptor chain, the effector RMSD may be artificially high
if mainEntry.hasChain(Chains.gprotein)
    effectorAtoms = effectorChain.getAtoms();
    RMSD_G = mainSim.computeRmsd(effectorAtoms);
    RMSF_G = mainSim.computeRmsf(effectorAtoms, 'StartFrame', settings.frames2skip + 1);
    bfactor(1, effectorAtoms) = mean(RMSF_G, 2);
    save(saveVarName,'RMSD_G','RMSF_G','-append');
end

% Save the RMSF as the B-factor
writepdb(fullfile(md2pathdir, "prot_rmsf.pdb"), mainEntry.pdb, [], 'default', bfactor);


% Options:
errStrideRMSF = 3; % how far are the error bars on the RMSF plot

plot_save_rmsd_rmsf;
% Plots RMSD/RMSF plots for receptor (and ligand) and saves them to md2pathdir


%% Calculate contact map of receptor with ligand
% If ligand/gp exists in simulation!

% Dual cutoff contact scheme: two atoms come into contact if within
% rcut1 but only disconnect when they're further than rcut2
rcut1 = 3.5;
rcut2 = 5;

if isempty(settings.bsResidues) && mainEntry.hasChain(Chains.ligand)
    add2log(md2pathdir, 'Calculating contact map of receptor with ligand');

    [contactsPerRun, chainBAtomIndices, ~] = mainSim.computeContacts(mainChain, ligandChain, 'StartFrame', settings.frames2skip + 1,'RCut1', rcut1, 'RCut2', rcut2);
    receptorLigandResIds = plot_save_contact_map(contactsPerRun, mainChain, ligandChain, chainBAtomIndices, ...
        'ImportantCutoff', 0.55, 'RCut', rcut2, 'SaveDir', md2pathdir,'SaveXlsSheet',"Full simulation"); %FIX MUMBER SHOWN ON TITLE OF HEATMAP

    save(saveVarName,'contactsPerRun','chainBAtomIndices','receptorLigandResIds','-append');
elseif ~isempty(settings.bsResidues)
    receptorLigandResIds = settings.bsResidues;
end


%% Calculate contact map of receptor with G protein

if isempty(settings.gpiResidues) && (mainEntry.hasChain(Chains.gprotein) || (length(database.entries) > 1))
    add2log(md2pathdir, "Calculating Gp interacting residues from downloaded PDB: " + settings.pdbCode);

    if length(database.entries) > 1
        refEntry = database.entries{2};
        [~, receptorGpResIds] = calc_plot_gp_contact_map(mainEntry, refEntry, Chains, 'SaveDir', md2pathdir, 'StartFrame', settings.frames2skip + 1);
    else
        [receptorGpResIds, ~] = calc_plot_gp_contact_map(mainEntry, [], Chains, 'SaveDir', md2pathdir, 'StartFrame', settings.frames2skip + 1);
    end

    % Write out gp contacting residues to file
    writematrix(receptorGpResIds', fullfile(md2pathdir, "GPI_residues.txt"), 'Delimiter', 'space');

    % Save figure as pdf and as .fig
    % figPath = fullfile(md2path, "contact_map" + name + "_G-protein");
    % savefig(figPath);
    % print2pdf(figPath, 1);
    save(saveVarName,'receptorGpResIds','chainBAtomIndices','receptorLigandResIds','-append');
    add2log(md2pathdir, ["G-protein contacting residues written to GPI_residues.txt", "Contact maps saved to pdf and fig", ""]);
elseif ~isempty(settings.gpiResidues)
    % Set the variable used by the script for Gp interacting residues
    receptorGpResIds = settings.gpiResidues;
end

%% PCA section:

if settings.doPCA
    add2log(md2pathdir, "PCA analysis commenced");
    pcadir = fullfile(md2pathdir, "pca");

    if ~exist(pcadir, 'dir')
        mkdir(pcadir);
        add2log(md2pathdir, "Creating pca directory in " + md2pathdir);
    end
    kPrinComp = 2; % Number of principal components to take into consideration
    rng(1); % For reproducibility
    
    % Initalize color index C
    C=[];
    nFramesEff = zeros(mainSim.runCount,1);
    for i = 1:mainSim.runCount
        % Take into consideration runs with different number of frames
        nFramesEff(i) = (size(mainSim.traj{i},1) - settings.frames2skip - 1)+1; % Frames used in calculations
        framesNdx = ((settings.frames2skip+1):size(mainSim.traj{i},1))';
        Ctemp = [i*ones(nFramesEff(i),1) framesNdx];
        C=[C ;Ctemp];% Used for coloring and for labeling: [ run frameNdx]
    end

    %% Do a PCA analysis on the ligand pose
    if  mainEntry.hasChain(Chains.ligand)
        % Ligand: we do PCA on heavy atoms in case of small ligand and CA in
        % case of peptide

       % Run the ligand PCA module: use trjTypeLigand to change what goes
       % into the PCA 
        pca_1_ligand;
        % output:
        % FILE: pcaCenterFrame_run.txt
        %  ligandReceptorpca.pdb
        %  ligandpca.pdb
        %
        % * pca_frame_centers: frames and runs of PCA cluster centers
        %
        save(saveVarName,'rmsdLigandpcaCenter','-append');
    end

    % Run the receptor PCA module: use trjType to change what goes into the
    % PCA
     pca_2_prot;
    % output: same as pca_ligand but for receptor

end


%% Calculate GPCR order parameters to assess state of the system
% Filter frames according to activation states if needed
if ~isfield(settings,'cullStates') 
settings.cullStates=true;
end

if settings.isGPCR
   [params, ndxXY] = calcPlotOrderParameters(database, database.entries{3}, 'SavePath', ...
       md2pathdir,'frames2skip',settings.frames2skip + 1,'cullFrames', settings.cullStates);
   if settings.doPCA
   [~, ndxXYClusters, clusterNdx_filtered] = calcPlotOrderParametersClusters(database, database.entries{3}, ...
       indexOfCluster_pca , 'SavePath', md2pathdir,'frames2skip',settings.frames2skip + 1 ,'cullFrames', settings.cullStates);
   save(saveVarName,'ndxXYClusters','clusterNdx_filtered','-append')
   end

   save(saveVarName,'params','ndxXY','-append');

end

%% Calculate dihedrals from trajs
if settings.includeLigPathwayCalc
    receptorResIds = mainChain.concatResIds(ligandChain);
else
    receptorResIds = mainChain.resIds;
end

mainSim.computeDihedrals(  ...
    'Path', fullfile(md2pathdir, "dihedrals.mat"), ...
    'ReSortPath', fullfile(md2pathdir, "reSort.txt"), ...
    'ResIds', receptorResIds, ...
    'StartFrame', settings.frames2skip + 1 ...
);


% %% Compare dihedral distributions of selected dihedrals to references
% % Ex: compare chi's of conserved rotamers in class A GPCRs
% 
% if settings.isGPCR
%     % Choose features
%     resImp = ["3.50","4.46","5.54","5.58","6.44","6.47","6.48","7.45"];
% %     resImp = ["3.32","5.42","5.46","7.43"];
%     dihPairs = [1 2; 2 3;3 4];
% 
%     for pair = 1:size(dihPairs,1)
%      [idxCell] = dihedralClusterReferences(database,resImp,'dihij',dihPairs(pair,:));
%      savefig(fullfile(md2pathdir,"dihedralClusters_"+num2str(pair)))
%     end
%     % To do: 
%     % 1D or 3D maps?
% end




%% Run pathway calculation for clustered data rather than full trajectory?
% Includes culling of frames

% 
% if settings.isGPCR
%     if settings.doPCA
%         % Include filtering from both criteria (if applicable)
%         clusterNdx_filtered_all = min(clusterNdx_filtered,[],2); % This works since 
%         % a row is either a zero or cluster number
%         
%         % Options for subruns:
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
%         
%         pca_1b_ligandPathCalcClusters;
%         clear clusterSims % Clear large variables
%     end
% end

%% Run pathway calculation for culled frames without clustering

if settings.isGPCR && settings.cullStates
    ndxXY_all = ((ndxXY{1}) & (ndxXY{2})); % Filter using both criteria
    % Make subset
    culledSims = mainSim.createSubset(C(ndxXY_all, :), 'DihedralMatIndices', ndxXY_all);

    if isempty(settings.bsResidues) && mainEntry.hasChain(Chains.ligand)
        % Contact map
        [contactsPerRunCulled, ligandAtomIndicesCulled] = culledSims.computeContacts(mainChain, ligandChain,'RCut1', rcut1, 'RCut2', rcut2);
        disp(mainChain)
        disp(ligandChain)
        disp(rcut1)
        disp(rcut2)
        %disp(ligandAtomIndicesCulled)
        disp(contactsPerRunCulled)
        % Plot contact map
        receptorLigandResIdsCulled = plot_save_contact_map(contactsPerRunCulled, mainChain, ligandChain, ligandAtomIndicesCulled, ...
            'ImportantCutoff', 0.55, 'RCut', rcut2, ...
            'SaveDir', md2pathdir, ...
            'SaveName', ("%s_culled" ), ...
            'SaveXlsName', "%s", ...
            'SaveXlsSheet', "Culled data");
    elseif ~isempty(settings.bsResidues)
     receptorLigandResIdsCulled = settings.bsResidues;
    end
    % Calculate MI
    miPath = fullfile(md2pathdir, "MI_culled.mat");
    culledSims.computeMI('Path', miPath);
    
    % MI convergence
    assessEntropyConvergence(culledSims, ...
        'Name', mainEntry.name + ", culled data", ...
        'SavePath', fullfile(md2pathdir, "Convergence_culled"));

    % Path calculation
pathCalcdir = prepareAlloPathCalc( ...
        culledSims, ...
    mainChain, ...
    settings, ...
    'Dir', md2pathdir, ...
        'Name', ("%s_Culled_data"), ...
    'LogPath', md2pathdir, ...
        'ReceptorLigandResIds', receptorLigandResIdsCulled, ...
    'ReceptorGpResIds', receptorGpResIds, ...
    'ReceptorResIds', receptorResIds ...
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
    clear culledSims; % Clear large variables
end

%% Calculate MI! (Including ALL frames)

mainSim.computeMI('Path', fullfile(md2pathdir, "MI.mat"));


%% Assess convergence of entropies (Including ALL frames)

assessEntropyConvergence(mainSim, 'SavePath', fullfile(md2pathdir, "convergence"));
add2log(md2pathdir, "S2 convergence calculated and plotted");


%% Prepare input files for path calculation! (Including ALL frames)

pathCalcdir = prepareAlloPathCalc( ...
    mainSim, ...
    mainChain, ...
    settings, ...
    'Dir', md2pathdir, ...
    'LogPath', md2pathdir, ...
    'ReceptorLigandResIds', receptorLigandResIds, ...
    'ReceptorGpResIds', receptorGpResIds, ...
    'ReceptorResIds', receptorResIds ...
);


%% Finally run path calculation!!

%These options are now specified in the input script
% Distance cutoff for considering "allosteric" signal, default is 10 A
disCutOff = 10;

% cutoff for pathways to be considered near, and amount of overlap required
% to cluster them into one pipeline, default to 7.5 and 0.75
nearcutoff = 7.5;
overlapcutoff = 0.75;

%Reinitialize some needed variables for alloPathCalc
if ~exist("md2pathdir",'var')
    md2pathdir = fullfile(settings.mydir, "md2pathdev");
end
if ~exist("pathCalcdir",'var')
    pathCalcdir = fullfile(md2pathdir, sprintf("%s", 'alloPathCalc'));
end

load(fullfile(pathCalcdir,"workspace.mat"))

prepMI;
MIStatsResLevel;
graphanalysis;
ClusterMIpathways;
analyzeClusters;
writeChannels;
analyzePathDomains;
MIanalysisBS2Effector;
if settings.isGPCR
    pathwayBiasAnalysis;
end
[pHere] = visualizeClsGraph(PDB,pathstruc,Gmatmajor,1,'MIFractionCutoff',MIFractionCutoff);
% If protein is inversed along z, use this command to make it upright again
% set(gca,'zdir','reverse')
% 
% If the "chirality" is off, you can use this:
% set(gca,'ydir','reverse')
