% Meta loads MI_delta_plot:

%% Settings:

% Load trajectories of WT and mutants, and then perform a common PCA to all
% of the common trajs

% Father directory for all the runs
metadir = 'D:\Telework_library\dopamine_phase_3';
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
pdbSaveDir = 'D:\Telework_library\dopamine_phase_3\md2pathMeta\dMI_norm_pdbs_240524';

pdbName = 'prot.pdb';

% Gi DA: 
foldersToStudy = {'2-d2_dop_T174M-C220L', ...
    '10-d2_dop_L214M','12-d2_dop_L92G','14-d2_dop_F217M','16-d2_dop_I125N', ...
    '19a-d2_dop_L92H'};
numRunsSys = [10 5 5 5 5 5 5]; % number of runs for each system
thisSysLabel = {'T5.54M-C6.47L','L6.41M','L3.41G','F6.44M','I4.46N', ...
    'L3.41H'};


% Mega BRC:
% foldersToStudy = {'3-d2_bromo_WT','4-d2_bromo_T174M-C220L','11-d2_brc_L214M', ...
%      '13-d2_brc_L92G','15-d2_brc_F217M','17-d2_brc_I125N', '18b-d2_brc_F217I', ...
%      '19b-d2_brc_L92H'};
% numRunsSys = [10 5 5 5 5 5 5]; % number of runs for each system
% thisSysLabel = {'WT','T5.54M-C6.47L','L6.41M','L3.41G','F6.44M','I4.46N', ...
%     'F6.44I','L3.41H'};

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
chains = 'A-B';

frames2skip = 500; % remove first 50 ns of the simulation
nCentersLig = 10; % Number of points closest to the mean point of a system 
% consider for analysis and pdb saving
% name = 'dd2_brc_June_29_2023';

% Analysis parameters:
nPipes = 10; % number of pipelines to consider in analysis
topHubs = 25; % Number of top scoring global hubs to consider
isGPCR = true;
% Generic GPCR numbering from https://gpcrdb.org/residue/residuetable
gpcrdbRefName = 'DD2R';

abc = char(65:65+25); % Can you sing the abc?

% Reference PDB code that this simulation is based on
pdbCode = '6VMS';
pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
pdbCodeInactive = '6CM4';
pdbInactiveChains = 'A--';

md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc_Culled_data';
%% D1 stuff:
% Load trajectories of WT and mutants, and then perform a common PCA to all
% of the common trajs

% Father directory for all the runs
metadir = 'D:\Telework_library\gpcrdb_extra_simulations\d1_dpa';
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
pdbSaveDir ='D:\Telework_library\dopamine_phase_3\md2pathMeta\dMI_norm_pdbs_240524\D1';

pdbName = 'prot.pdb';

% DA-Gs
foldersToStudy = {'c-d1r_dpa_pam_gp','d-dpa_I125N_gp','e-dpa-F228M_gp'};
numRunsSys = [5 5 5 5 5 5 5]; % number of runs for each system
thisSysLabel = {'DA-WT-PAM-Gs','DA-I4.46N-Gs','DA-F6.44M-Gs', ...
   };

% BRC-Gs
% foldersToStudy = {'g-d1r_brc_gp/b-d1_BRC_from_MD','h-d1r_brc_F228M_gp','i-d1r_brc_F228I_gp' };
% numRunsSys = [5 5 4 5 5 5 5]; % number of runs for each system
% thisSysLabel = { 'BRC-WT-Gs','BRC-F6.44M-Gs','BRC-F6.44I-Gs'};




% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
chains = 'ACB';

frames2skip = 500; % remove first 50 ns of the simulation
nCentersLig = 10; % Number of points closest to the mean point of a system 
% consider for analysis and pdb saving
% name = 'dd1_all_July_19_2023';

% Analysis parameters:
nPipes = 10; % number of pipelines to consider in analysis
topHubs = 25; % Number of top scoring global hubs to consider
isGPCR = true;
% Generic GPCR numbering from https://gpcrdb.org/residue/residuetable
gpcrdbRefName = 'DD1R';

abc = char(65:65+25); % Can you sing the abc?

% Reference PDB code that this simulation is based on
pdbCode = '7CKZ';
pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
pdbCodeInactive = '6CM4';
pdbInactiveChains = 'A--';

md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc';
useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment

%% 

for thisSys = 1:length(foldersToStudy)
    
    if ~isfield(settings,'frames2skipRef') 
        settings.frames2skipRef=settings.frames2skip;
    end
    1

    dirHere = fullfile(metadir,foldersToStudy{thisSys})

    run(fullfile(dirHere,"input_kldiv.m"));
    
    %% Initialize variables:

    
    if ~isfield(settings,'includeLigPathwayCalc') 
        settings.includeLigPathwayCalc = false;
    end
    
    if ~isfield(settings,'pdbCodeExtra') 
        settings.pdbCodeExtra = [];
    end
    
    klDivFileName = ['klDiv_' settings.refName  '.mat'];
    %% Load, fetch and align PDB files and trajectories
    database = Database(settings.databasePath);
    
    % Chains: [receptor, G protein, ligand]
    database.read(fullfile(settings.mydir, "prot.pdb"), settings.chains, settings.mainName);
    % Load the KLDiv reference structure:
    database.read(fullfile(settings.refdir, "prot.pdb"), settings.chains, settings.refName );
    
    % Read reference structures
    if ~isempty(settings.pdbCode)
        database.fetch(settings.pdbCode, settings.pdbChains);
    
        if ~isempty(settings.pdbCodeInactive)
            database.fetch(settings.pdbCodeInactive, settings.pdbInactiveChains);
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
    refEntry = database.entries{2}; % This entry spot usually reserved for active ref pdb
    
    mainChain = mainEntry.chains{Chains.receptor};
    refChain = refEntry.chains{Chains.receptor};
    
    if mainEntry.hasChain(Chains.ligand)
        ligandChain = mainEntry.chains{Chains.ligand};
    end
    
    % Transform xtcs to dcds if they're not there yet:
    dirHere = settings.mydir;
    for j=1:2
        areThereDCDs = dir(fullfile(dirHere, "run*", "traj.dcd"));
        
        if length(areThereDCDs) < settings.numRuns
            % run VMD from command line:  vmd -dispdev text -e
            pathToScript = fullfile(pwd(), "load_save.tcl");  % assumes script is in current directory
            cmdStr       = "vmd -dispdev text -e " + pathToScript + " -args " + dirHere + " " + num2str(settings.numRuns) + " " + num2str(settings.stride) + " " + database.entries{j}.path + " " + settings.xtcName;
            system(cmdStr);
    %         add2log(md2pathdir, "Transformed .xtc files to .dcd files with the following command: """ + cmdStr + """");
        end
        dirHere = settings.refdir;
    end
    % Load trajs into sim class
    mainSim = mainEntry.addSimulation(fullfile(settings.mydir, "run*"));
    refSim = refEntry.addSimulation(fullfile(settings.refdir, "run*"));
    
    md2pathdir = fullfile(settings.mydir, "md2pathdev");
    
    % Create md2pathdir if it does not exist
    if ~exist(md2pathdir, 'dir')
        mkdir(md2pathdir);
    end
    %% Add labels
    
    % Add labels with Ballesteros-Weinstein notation, aligned with the 2nd
    % database entry
    % Exported from GPCRdb (https://gpcrdb.org/residue/residuetable)
    if length(database.entries) > 1 && isfield(settings, 'systemName') && settings.isGPCR
        residueTablePath = fullfile(database.dir, settings.systemName + "_residue_table.xlsx");
        database.label(3, residueTablePath); % Read 3rd entry (supposedly reference)
    end
    %% Calculate dihedrals from trajs, assumes same receptor IDs between entries!
    
    % Test system
    if settings.includeLigPathwayCalc
        receptorResIds = mainChain.concatResIds(ligandChain);
        receptorResIdsRef = refChain.concatResIds(ligandChain);
    else
        receptorResIds = mainChain.resIds;
        receptorResIdsRef = refChain.resIds;
    end
    
    % if settings.alignTestRefDih % When Test and ref sims have different residue numbers/length
    %     receptorResIds = database.residues{1}.output1(find(database.residues{1}.output2));
    % end
    
    mainSim.computeDihedrals( ...
        'Path', fullfile(md2pathdir, "dihedrals.mat"), ...
        'ReSortPath', fullfile(md2pathdir, "reSort.txt"), ...
        'ResIds', receptorResIds, ...
        'StartFrame', settings.frames2skip + 1 ...
    );
    
    % Reference system
    refSim.computeDihedrals( 'Path', fullfile(settings.refdir, "md2pathdev","dihedrals.mat"),...
        'ReSortPath', fullfile(md2pathdir, "reSortRef.txt"), ...
        'ResIds', receptorResIdsRef, ...
        'StartFrame', settings.frames2skipRef  + 1 ...
    );

    %% Load the MI
    if ~exist('alloPathCalcName','var')
        alloPathCalcName = 'alloPathCalc_Culled_data';
    end
    temp = load(fullfile(settings.mydir,'md2pathdev',alloPathCalcName,'workspace.mat'),'-mat','MIres');
    MIres = temp.MIres;
    temp = load(fullfile(settings.refdir,'md2pathdev',alloPathCalcName,'workspace.mat'),'-mat','MIres');
    MIresWT = temp.MIres;
    helices = settings.helices;
    
    paperColors; % Loads paperColorPalette
    


    sMIres = sum(MIres);
    sMIresWT = sum(MIresWT);
    
    deltaData = sMIres-sMIresWT;
    deltaDataNorm = zeros(length(sMIresWT),1);
    for i = 1:length(sMIresWT)
        deltaDataNorm(i) =  (sMIres(i) - sMIresWT(i))/(0.5*(sMIresWT(i) +sMIres(i)));
    end
    %% Save deltaMI as B-factor:
    
    bfactor = zeros(1, mainEntry.atomCount);
    bfactordNorm = zeros(1, mainEntry.atomCount);
    
    if settings.includeLigPathwayCalc % Include ligand in b-factor array
        chainSelect = selectname(mainEntry.pdb.chainid, settings.chains(Chains.receptor)) | selectname(mainEntry.pdb.chainid, settings.chains(Chains.ligand));
        selection = selectname( mainEntry.pdb.name,'CA') & chainSelect;
    else
        selection = selectname( mainEntry.pdb.name,'CA') & selectname(mainEntry.pdb.chainid, settings.chains(Chains.receptor));
    end
    % resFromReSort = selectid(mainEntry.pdb.resseq,unique(reSortCommon(:,1))); % Only residues contained in reSortCommon
    % selection = selection & resFromReSort;
    
    bfactor(selection) = deltaData;
    writepdb(fullfile(pdbSaveDir, ['prot_dMI_' mainEntry.name '_' refEntry.name '.pdb']), mainEntry.pdb, [], 'default', bfactor);
    
    bfactordNorm(selection) = deltaDataNorm;
    writepdb(fullfile(pdbSaveDir, ['prot_dMINorm_' mainEntry.name '_' refEntry.name '.pdb']), mainEntry.pdb, [], 'default', bfactordNorm);

end