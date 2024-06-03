%% Prepare input files: common input files for all "meta" scripts
% Every "meta" script needs a different set of input, as shown below:

databasePath = 'path/to/md2pathDatabase/'; 
% Father directory for all the runs
metadir = 'path/to/simulation/meta/dir';
foldersToStudy = {'1-D2_DA_WT','2-D2_BRC_WT'};

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
isGPCR = true;
saveData = true;
plotrmsd = false ;

%% Parameters for metaBondDist:

isGPCR = true;
frames2skip = 500; % 50 ns
useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment
attempt2AlignTrajs = false;


%% Parameters for metapca_4_multisys: Now in its own script!!!!!
