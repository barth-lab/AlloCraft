%% Instructions
% Copy this script to your running directory and run it before running
% kldivMain.m


%% Main input variables:
global settings;

% Get the full path to the currently running script
current_script_path = mfilename('fullpath');
% Get the directory of the currently running script
[current_script_dir, ~, ~] = fileparts(current_script_path);
% Define the project root as the parent directory of the script's directory
PROJECT_ROOT = fileparts(current_script_dir);

% Main data directory for perturbed and reference directories
% settings.mydir = '../data/DA_F6.44M_S5.46G';
% settings.refdir =  '../data/1-d2_dop_WT';
settings.mydir = fullfile(PROJECT_ROOT, 'data', 'DA_F6.44M_S5.46G');
settings.refdir =  fullfile(PROJECT_ROOT, 'data', '1-d2_dop_WT');
settings.mainName = ['F6.44_S5.46G'];
settings.refName = 'WT';

% Add necessary directories to path
% addpath('src');
% addpath('MDprot/mdtoolbox');
% addpath('MDprot/Useful_tools');
% addpath('util');
% addpath('MDprot/GPCR_tools');
% addpath('alloPathCalc');
% addpath('MDprot/matdcd-1.0');
% addpath('MDprot');
% addpath('extra_scripts');
addpath(fullfile(current_script_dir, 'src'));
addpath(fullfile(current_script_dir, 'MDprot', 'mdtoolbox'));
addpath(fullfile(current_script_dir, 'MDprot', 'Useful_tools'));
addpath(fullfile(current_script_dir, 'util'));
addpath(fullfile(current_script_dir, 'MDprot', 'GPCR_tools'));
addpath(fullfile(current_script_dir, 'alloPathCalc'));
addpath(fullfile(current_script_dir, 'MDprot', 'matdcd-1.0'));
addpath(fullfile(current_script_dir, 'MDprot'));
addpath(fullfile(current_script_dir, 'extra_scripts'));


% Path for the database used by md2path
% settings.databasePath = '../data/Databases';
settings.databasePath = fullfile(PROJECT_ROOT, 'data', 'Databases');
settings.systemName = 'DD2R'; % Points to the residue table in the database

% Name of the trajectory file in each run
settings.xtcName = 'traj.dcd';

% Number of frames to remove at the start of the simulation
settings.frames2skip = 500;

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
settings.numRuns = 5;

% How many frames to skip when saving xtcs to dcds
settings.stride = 1;

% Numbering based on PDB in mydir, NOT refdir
settings.helices = [ [1 31]
    [36 67]
    [72 107]
    [113 142]
    [155 195]
    [197 234]
    [239 264]
    ];

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

% Specific residues to be highlighted during KL1 visualziation and the
% corresponding text, numbering based on PDB in mydir, NOT refdir
settings.highlightRes =  [36 38 39 101 104 105 108 109 111 112 115 181 185 188 197 198 201 202 205 206 209 210 213 264 265];
settings.highlightText = 'GPI';

% Actiavtes GPCR specfic options
settings.isGPCR = true;
settings.includeLigPathwayCalc = false;