%% Instructions
% Copy this script to your running directory and run it before running
% kldivMain.m


%% Main input variables:

% Main data directory for perturbed and reference directories
settings.mydir = '/path/to/AlloCraft/MD_simulation_analysis/demo/2-D2_BRC_WT';
settings.refdir =  '/path/to/AlloCraft/MD_simulation_analysis/demo/1-D2_DA_WT';
settings.mainName = 'BRC-WT';
settings.refName = 'DA-WT';

% Path for the database used by md2path
settings.databasePath = '/path/to/AlloCraft/MD_simulation_analysis/demo/database';
settings.systemName = 'DD2R'; % Points to the residue table in the database

% Name of the trajectory file in each run
settings.xtcName = 'step7_noPBC_prot.xtc';

% Number of frames to remove at the start of the simulation
settings.frames2skip = 0;

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'ACB';

% Reference PDB code that this simulation is based on
settings.pdbCode = '6VMS';
settings.pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
settings.pdbCodeInactive = '6CM4';
settings.pdbInactiveChains = 'A--';
% Number of runs
settings.numRuns = 2;

% How many frames to skip when saving xtcs to dcds
settings.stride = 1;

settings.helices = [   2    30
    37    67
    73   107
   118   143
   156   189
   198   234
   240   264];

% KLDiv options

% number of blocks to divide the simulation for signficance testing,
% recommended number is the number of independent simulations or less.
% MUST BE an EVEN number
settings.nBlocks = 6; 

% Histogram finite size effect correction
settings.Grassberger = false ;

% Significance threshold for KLDiv, recommended: 0.1 for nBlocks = 4
% 0.05 for nBlocks >= 6
settings.st = 0.05;

% Actiavtes GPCR specfic options
settings.isGPCR = true;

% Specific residues to be highlighted during KL1 visualziation and the
% corresponding text
settings.highlightRes =  [36    38    39   101   104   105   108   109   111   112   115   181   185   188   197   198   201   202   205   206   209   210   213   264   265];
settings.highlightText = 'G-p binding';
