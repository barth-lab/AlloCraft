%% Instructions
% Copy this script to your running directory and run it before running
% kldivMain.m


%% Main input variables:

% Main data directory for perturbed and reference directories
settings.mydir = 'D:\Telework_library\Chemokine_project\6-CXCR4_active_models_mutants\c-CXCR4_active_Fusion_V3Y_Y7L';
settings.refdir =  'D:\Telework_library\Chemokine_project\6-CXCR4_active_models_mutants\a-WT_CXCR4_active_model';
settings.mainName = 'Fusion';
settings.refName = 'WT';

% Path for the database used by md2path
settings.databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\';
settings.systemName = 'CXCR4'; % Points to the residue table in the database

% Name of the trajectory file in each run
settings.xtcName = 'step7_noPBC_prot.xtc';

% Number of frames to remove at the start of the simulation
settings.frames2skip = 500;

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'A-B';

% Reference PDB code that this simulation is based on
settings.pdbCode = '6LFM';
settings.pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
settings.pdbCodeInactive = '4RWS';
settings.pdbInactiveChains = 'A--';

% Which pdb should be used as reference for generic GPCR numbering?
% 2 means settings.pdbCode will be taken as reference
% 3 means settings.pdbCodeInactive will be
settings.refPDBNdx = 2; 
% Number of runs
settings.numRuns = 7;

% How many frames to skip when saving xtcs to dcds
settings.stride = 5;

% Numbering based on PDB in mydir, NOT refdir
settings.helices = [5 39
46 74
80 114
127 149
168 202
207 242
248 276];

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

% Calculate 2nd order KL and mutual divergence
% NOTE: VERY slow calculation and mostly redundant with KL1
settings.calc2ndOrder = false;

% Specific residues to be highlighted during KL1 visualziation and the
% corresponding text, numbering based on PDB in mydir, NOT refdir
settings.highlightRes =  [19    46    47    51    52    59    60    79    82    86    91   117   122   129   159   183   187   209   227   243   255];
settings.highlightText = 'Valines';

% Activates GPCR specfic options
settings.isGPCR = true;
settings.includeLigPathwayCalc = true;