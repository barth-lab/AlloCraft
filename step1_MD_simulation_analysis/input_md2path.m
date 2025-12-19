%% Instructions
% Copy this script to your running directory and run it before running
% md2path.m

% Preparations needed: PDB file with separate chains for receptor, ligand,
% and G-protein: 1st chain is receptor, 2nd is ligand, 3rd is Gp
% Input: directory, number of runs to analyze, and stride of trajectory
% NOTE: all input residues should be in the same pdb numbering as the input
% PDB


%% Main input variables:
global settings;

% Get the full path to the currently running script
current_script_path = mfilename('fullpath');
% Get the directory of the currently running script
[current_script_dir, ~, ~] = fileparts(current_script_path);
% Define the project root as the parent directory of the script's directory
PROJECT_ROOT = fileparts(current_script_dir);


% Main data directo
% settings.mydir = '../data/DA_F6.44M_S5.46G';
settings.mydir = fullfile(PROJECT_ROOT, 'data', 'DA_F6.44M_S5.46G');

% Add src directory to path
% addpath('src');
% addpath('MDprot/mdtoolbox');
% addpath('MDprot/Useful_tools');
% addpath('util');
% addpath('MDprot');
% addpath('MDprot/matdcd-1.0');
% addpath('MDprot/Entropy_MI');
% addpath('MDprot/GPCR_tools');
% addpath('alloPathCalc');
% addpath('extra_scripts');
addpath(fullfile(current_script_dir, 'src'));
addpath(fullfile(current_script_dir, 'MDprot', 'mdtoolbox'));
addpath(fullfile(current_script_dir, 'MDprot', 'Useful_tools'));
addpath(fullfile(current_script_dir, 'util'));
addpath(fullfile(current_script_dir, 'MDprot'));
addpath(fullfile(current_script_dir, 'MDprot', 'matdcd-1.0'));
addpath(fullfile(current_script_dir, 'MDprot', 'Entropy_MI'));
addpath(fullfile(current_script_dir, 'MDprot', 'GPCR_tools'));
addpath(fullfile(current_script_dir, 'alloPathCalc'));
addpath(fullfile(current_script_dir, 'extra_scripts'));


% Path for the database used by md2path
% settings.databasePath = '../data/Databases/';
settings.databasePath = fullfile(PROJECT_ROOT, 'data', 'Databases');
settings.systemName = 'DD2R'; % Points to the residue table in the database



%% Define the project root directory
% PROJECT_ROOT = '/home/asengar/Downloads/trial_github/AlloDy_Analysis_Repo'; 

% Main data directory
%settings.mydir = fullfile(PROJECT_ROOT, 'data', 'DA_F6.44M_S5.46G');

% Add directories to path
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'src'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'MDprot', 'mdtoolbox'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'MDprot', 'Useful_tools'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'util'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'MDprot'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'MDprot', 'matdcd-1.0'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'MDprot', 'Entropy_MI'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'MDprot', 'GPCR_tools'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'alloPathCalc'));
%addpath(fullfile(PROJECT_ROOT, 'scripts', 'extra_scripts'));

% Path for the database used by md2path
%settings.databasePath = fullfile(PROJECT_ROOT, 'data', 'Databases', filesep);
%settings.systemName = 'DD2R'; % Points to the residue table in the database



% Name of the trajectory file in each run
%settings.xtcName = 'step7_noPBC_prot.xtc';
settings.xtcName = 'traj.dcd';

% Number of frames to remove at the start of the simulation
settings.frames2skip =1000;

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'ACB';

% Reference PDB code that this simulation is based on
settings.pdbCode = '6VMS';
settings.pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
settings.pdbCodeInactive = '6CM4';
settings.pdbInactiveChains = 'A--';

% Extra PDBs for reference (in conformational landscape plots for example)
% Note: double check the alignment in database.residues{1}: if a
% misalignment exists in the residues used from conf. landscape
% calculation, there will an error!
% settings.pdbCodeExtra = {'6LFO','7F1R'};
% settings.pdbChainsExtra = {'RAD','RA-'};

% Number of run3
settings.numRuns = 3;

% How many frames to skip when saving xtcs to dcds
settings.stride = 1;

% Actiavtes GPCR specfic options (calculating 7 TM helices and
% GPCR allosteric path calculation, as well as planned GPCR order parameters)
settings.isGPCR = true;

% Residues to be ignored in the pathway calculation, usually you want to
% ignore floppy termini or long loops that introduce noise into the system.
% Defaults to first and last 4 residues of the receptor chain. Assumes that
% protein-renum, true to its name, is renumbered starting from 1.
settings.excludedResidues = [];

% The residues belonging to the 7 helical turns of a GPCR, If
% empty, the code will attempt to find them by calling VMD and extracting a
% secondary structure variable.
settings.helices = [ [2 32]
    [37 68]
    [73 108]
    [114 143]
    [156 196]
    [198 235]
    [240 265]
    ];

settings.helices = [ [1 31]
    [36 67]
    [72 107]
    [113 142]
    [155 195]
    [197 234]
    [239 264]
    ];

% If binding site/ligand interacting residues are not defined,
% the code will calculate contact maps of the protein and ligand.
settings.bsResidues = [];

% If intracellular binding partner is present, check the contacts from the
% original PDB structure and NOT from the simulation, since G-protein is
% usually truncated in simulations and some contacts are lost.
settings.gpiResidues = [];

% Same as gpi but for beta arrestin, keep empty unless you want to do a beta
% arrestin bias calculation.
settings.baiResidues = [];

% Do principle component analysis on the ligand binding pose as
% well as the receptor conformations
settings.doPCA = true;

% Number of clusters for PCA analysis of the ligand, if empty,
% the code will calculate it automatically, BUT it will take time!!!
settings.kClusters = [];

% Maximum number of clusters to consider during PCA clustering,
% only give this option if kClusters is empty!!!
settings.kmax = 15;

% Type of input to pca_2_prot to do PCA on the protein trajectory
% Could be: 'CA', 'DihAll', 'DihBB', 'DihSC', and 'Distances'
% 'Distances' requires additional options: settings.NdisPCA and
% settings.CGPCA
settings.trjType = 'CA';

% Type of input to pca_1_ligand to do PCA on ligand trajectory:
% Could be 'Default', 'VectorsCoord', and 'DistancesCoord'
% 'Default' gives CAs if ligand is peptide or not H atoms if ligand is
% small
settings.trjTypeLigand = 'Default'; 

% Remove frames that do not belong to desired states from dihedral/MI
% calculation
% WARNING: requires user input during run at the "Calculate GPCR order
% parameters" section
settings.cullStates = true;

% Include ligand (peptide only for now) in pathway calculations
% NOTE: take extra care when using this option, as the pathway part of
% the code does not deal well with different chains, to avoid problems, do
% the following:
% 1- Make sure there are no residue numbers repeated in the peptide and
% receptor
% 2- Renumber both receptor and peptide in the same order in which they
% appear in the pdb file
settings.includeLigPathwayCalc = false;


% Allosteric pathway run options (Options controlled from md2path script now)

% Distance cutoff for considering "allosteric" signal, default is 10 A
settings.disCutoff = 10;

% Minimal fraction of MI to be conserved when filtering pathways
settings.miFractionCutoff = 0.85;

% Cutoff for pathways to be considered near to cluster them into one pipeline
settings.nearCutoff = 7.5;

% Amount of overlap required for pathways to be considered near
settings.overlapCutoff = 0.75;

% Diagnostic plots on?
settings.diagnosticsOn = true;

% Weight pathways by MI post-clustering? (WARNING: experimental)
settings.MIWeightPaths = true;