%% Prepare input files: common input files for all "meta" scripts
% Every "meta" script needs a different set of input, as shown below:

databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
% Father directory for all the runs
metadir = 'D:\Telework_library\dopamine_phase_3';
foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT'};

thisSysLabel = {'DA WT','BRC WT'};

name = 'DD2R_WT';

% (Optional)
% Reference PDB code that this simulation is based on
pdbCode = '6VMS';
pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
pdbCodeInactive = '6CM4';
pdbInactiveChains = 'A--';

refNdx = 1;
%% Parameters for metaPathCalcAnalysis.m
abc = char(65:65+25); % Can you sing the abc?

% Analysis parameters:
nPipes = 10; % number of pipelines to consider in analysis
topHubs = 25; % Number of top scoring global hubs to consider
isGPCR = true;

md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc_Culled_data';
%% Parameters for metaContactMap.m
% Generic GPCR numbering from https://gpcrdb.org/residue/residuetable
% You need to DL the table manually from the website for every receptor you
% want to study and save it as:
%
% gpcrdbRefName_residue_table.xlsx
%
% Replace gpcrdbRefName with the system name defined below, of course!
gpcrdbRefName = 'DD2R';

useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment

%% Parameters for metapca_4_multisys and metaDihedralAnalysisProtorype.m

numRunsSys = [5 5]; % Number of runs available for each system studied
frames2skip = 500; % remove first 50 ns of the simulation
nCentersLig = 10; % Number of points closest to the mean point of a system 
nCenters = 10;
% consider for analysis and pdb saving
pdbName = 'prot.pdb';
chains ='A--';
attempt2AlignTrajs = true; % Attempts to align trajectories from runs with 
% different number of residues, prone to error!!! ALWAYS DOUBLE CHECK,
% currently works if: 
% 1- residue numbers start from 1
% 2- residue numbers of same chain are sorted in increasing order (no funny 
% PDB files where residue 500 is before residue 1 in the same chain)
