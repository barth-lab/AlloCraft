% metaKLDiv: meta-analysis of KLDiv from several directories
% 
% % Input section
% settings.databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
% % Father directory for all the runs
% settings.metadir = 'D:\Telework_library\beta1_project';
% % 
% % 
% % % For initial comparison:
% % settings.foldersToStudy = {'4-4amj_TS_cvd','2-beta1_iso_nb80_TS_A5.58Y-L7.53Y','1-beta1_iso_nb80_TS','6-4amj_TS_alp', '7-4amj_TS_cyn','10-6h7l_dbt'};
% % 
% % settings.thisSysLabel = {'CVD','ISO-YY-Nb','ISO-Nb','ALP','CYN','DBT-Nb'};
% % % 
% % 
% % % For a more complete comparison:
% % % settings.foldersToStudy = {'1-beta1_iso_nb80_TS','2-beta1_iso_nb80_TS_A5.58Y-L7.53Y', ...
% % %     '4-4amj_TS_cvd','5-4amj_TS_iso','6-4amj_TS_alp', '7-4amj_TS_cyn','8-4amj_TS_dbt','9-2y03_TS_iso','11-6h7l_YY_dbt','12-6h7o_cyn','13-6h7o_YY_cyn'};
% % 
% % % settings.thisSysLabel = {'ISO-nb80','ISO-nb80YY','CVD','ISO','ALP','CYN','DBT', 'ISO-2y03', ...
% % %     'DBT-nb80YY','CYN-nb80','CYN-nb80YY'};
% % 
% % 
% % For comparison with Grahl et al. (2020): 41467_2020_15864_MOESM3_ESM.xlsx
% 
% % TS only no Nb
% % settings.foldersToStudy = {'6-4amj_TS_alp','4-4amj_TS_cvd','18-4bvn_TS_cyn','16-2y01_TS_dbt','9-2y03_TS_iso'};
% % 
% % settings.thisSysLabel = {'ALP','CVD','CYN','DBT','ISO'};
% % settings.thisSysNumExp = [1 3 4 5 6 ]; % For experimental comparison
% 
% % % Using Nb bound structures (obsolete)
% % % settings.foldersToStudy = {'6-4amj_TS_alp','4-4amj_TS_cvd','12-6h7o_cyn','10-6h7l_dbt','1-beta1_iso_nb80_TS', ...
% % %     '13-6h7o_YY_cyn','11-6h7l_YY_dbt','2-beta1_iso_nb80_TS_A5.58Y-L7.53Y'};
% % % 
% % % settings.thisSysLabel = {'ALP','CVD','CYN-nb','DBT-nb','ISO-nb','CYN-nb-YY','DBT-nb-YY','ISO-nb-YY'};
% % % settings.thisSysNumExp = [1 3 4 5 6 10 11 12]; % For experimental comparison
% % 
% % % Using inactive structures
% % % settings.foldersToStudy = {'6-4amj_TS_alp','4-4amj_TS_cvd','18-4bvn_TS_cyn','16-2y01_TS_dbt','9-2y03_TS_iso', ...
% % %    '20-4amj_YY_alp','15-4amj_YY_cvd', '19-4bvn_YY_cyn' ,'17-2y01_YY_dbt','14-2y03_YY_iso'};
% % % 
% % % settings.thisSysLabel = {'ALP','CVD','CYN','DBT','ISO','ALP-YY','CVD-YY','CYN-YY','DBT-YY','ISO-YY'};
% % % settings.thisSysNumExp = [1 3 4 5 6 7 9 10 11 12]; % For experimental comparison
% %
% % % YY variants only:
% % % settings.foldersToStudy = { ...
% % %    '20-4amj_YY_alp','15-4amj_YY_cvd', '19-4bvn_YY_cyn' ,'17-2y01_YY_dbt','14-2y03_YY_iso'};
% % % 
% % % settings.thisSysLabel = {'ALP-YY','CVD-YY','CYN-YY','DBT-YY','ISO-YY'};
% % % settings.thisSysNumExp = [ 7 9 10 11 12]; % For experimental comparison
% %
% % % For comparison with more recent assignments: 15NV_b1ARYY_peaklists.xlsx
% % % settings.foldersToStudy = {'2-beta1_iso_nb80_TS_A5.58Y-L7.53Y','14-2y03_YY_iso','15-4amj_YY_cvd','20-4amj_YY_alp', '19-4bvn_YY_cyn','17-2y01_YY_dbt',};
% % % 
% % % settings.thisSysLabel = {'ISO-nb80-YY','ISO-YY','CVD-YY','ALP-YY','CYN-YY','DBT-YY'};
% % % settings.thisSysNumExp = [1 3 4 5 7 8]; % For experimental comparison
% % 
% % 
% % % (almost) full set for PCA:
% % settings.foldersToStudy = {'6-4amj_TS_alp','4-4amj_TS_cvd','18-4bvn_TS_cyn','16-2y01_TS_dbt','9-2y03_TS_iso', ...
% %    '20-4amj_YY_alp','15-4amj_YY_cvd', '19-4bvn_YY_cyn' ,'17-2y01_YY_dbt','14-2y03_YY_iso','13-6h7o_YY_cyn', ...
% %    '11-6h7l_YY_dbt','2-beta1_iso_nb80_TS_A5.58Y-L7.53Y'};
% % 
% % settings.thisSysLabel = {'ALP','CVD','CYN','DBT','ISO','ALP-YY','CVD-YY','CYN-YY','DBT-YY','ISO-YY', ...
% %     'CYN-YY-nb','DBT-YY-nb','ISO-YY-nb'};
% % % 
% % 
% % settings.name = 'b1AR_TS';
% 
% settings.gpcrdbRefName = 'b1AR_wild_turkey';
% 
% % Reference system for the KLDiv calculations
% settings.refdir =  'D:\Telework_library\beta1_project\3-4amj_TS_apo';
% settings.refName = 'Apo';
% 
% % Chains in this order: [receptor, G protein, ligand]
% % Missing chains can be set as '-'.
% settings.chains = 'A--';
% 
% % Reference PDB code that this simulation is based on
% settings.pdbCode = '6H7J';
% settings.pdbChains = 'AC-';
% 
% % Inactive state reference PDB, used for dihedral comparison
% settings.pdbCodeInactive = '4AMJ';
% settings.pdbInactiveChains = 'B--';
% 
% settings.isGPCR = true;
% settings.includeLigPathwayCalc = false;
% settings.md2pathName = 'md2path';

% Which types of dihedrals to be considered for KL1 calculation:
% phi = 1; psi = 2; chi1 = 3; chi2 = 4; ... chi5 = 7
% settings.dihedralTypes = [1 2 3 4 5 6 7]; 
% settings.dihedralTypes = 7; % For now called from input script
%% Make the meta-direcory
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end
md2pathdirMeta = [settings.metadir slash 'md2pathMeta' slash];
% md2pathdir = [metadir slash 'md2pathMeta' slash name slash];
md2pathdir = [settings.metadir slash 'md2pathMeta' slash];

% Create meta directory
if ~exist(md2pathdirMeta, 'dir')
    mkdir(md2pathdirMeta)
    add2log(md2pathdirMeta, ['Creating md2path directory in ' settings.metadir]);
end

% Create subdirectories for every bundle of systems studied
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir)
    add2log(md2pathdir, ['Creating md2path directory in ' md2pathdirMeta]);
end

%% Load, fetch and align PDB files and trajectories
database = Database(settings.databasePath);

% Read PDB of every system to be studied:
for thisSys = 1:length(settings.foldersToStudy)
    % Get local directory
    mydir = ([settings.metadir slash settings.foldersToStudy{thisSys}]); 
    database.read(fullfile(mydir, "prot.pdb"), settings.chains, settings.thisSysLabel{thisSys});

end

% Load the KLDiv reference structure:
database.read(fullfile(settings.refdir, "prot.pdb"), settings.chains, settings.refName );
refEntry = database.entries{(length(settings.foldersToStudy)+1)}; % This entry spot usually reserved for active ref pdb

% Store names for each chain
Chains.receptor = 1;
Chains.gprotein = 2;
Chains.ligand = 3;

refChain = refEntry.chains{Chains.receptor};
refEntryNdx = length(settings.foldersToStudy) + 1;

if refEntry.hasChain(Chains.ligand)
    ligandChain = refEntry.chains{Chains.ligand};
end

% Get list of residues
if settings.includeLigPathwayCalc
    receptorResIds = refChain.concatResIds(ligandChain);
else
    receptorResIds = refChain.resIds;
end

% Read reference PDB structures
if ~isempty(settings.pdbCode)
    database.fetch(settings.pdbCode, settings.pdbChains);

    if ~isempty(settings.pdbCodeInactive)
        database.fetch(settings.pdbCodeInactive, settings.pdbInactiveChains);
    end
end

% Align sequences and produce the database.residues table
database.align();

% Align structures to the first database entry
database.alignStructures();

% Add labels with Ballesteros-Weinstein notation, aligned with the ref.
% database entry from the PDB
% Exported from GPCRdb (https://gpcrdb.org/residue/residuetable)
if length(database.entries) > (length(settings.foldersToStudy)+1) && isfield(settings, 'gpcrdbRefName') && settings.isGPCR
    residueTablePath = fullfile(database.dir, settings.gpcrdbRefName + "_residue_table.xlsx");
    database.label(length(settings.foldersToStudy)+2, residueTablePath); 
end

%% Read KLDiv from different systems:
klDivFileName = ['klDiv_' settings.refName  '.mat'];

kl1 = cell(length(settings.foldersToStudy),1);
kl1Res = cell(length(settings.foldersToStudy),1); % Per residue 
reSortCommon = cell(length(settings.foldersToStudy),1);

for thisSys = 1:length(settings.foldersToStudy)
    mydir = ([settings.metadir slash settings.foldersToStudy{thisSys} slash settings.md2pathName]); 
%     mydir = ([metadir slash foldersToStudy{thisSys} slash]); 
    temp =  importdata([mydir slash klDivFileName]);
%     reSortCommon{thisSys} = importdata([mydir slash 'reSortKLDiv.mat']);
     tempStruct = load([mydir slash 'reSortKLDiv.mat'],'reSortCommon');
    reSortCommon{thisSys} = tempStruct.reSortCommon;
    if isstruct(temp)
        kl1{thisSys} = temp.kl1;
    else
        kl1{thisSys} = temp;
    end
    
    % Get custom residue wise KL1: (according to specified dihedral types)
    kl1Res{thisSys} = zeros(length(receptorResIds),1);

    for i = 1:length(receptorResIds)
        resHere = receptorResIds(i); % Residue in ref. system
        dbRow = find(database.residues{1}{:,refEntryNdx} == resHere);
        resTest = database.residues{1}{dbRow,thisSys};


        ndxHere = find(reSortCommon{thisSys}(:,1) == resTest); % Index for dihedrals of this residue
%         reSortHere = (reSortCommon{thisSys}(ndxHere,2));
        
%         finalNdx = [];
%         % Dih Types:
%         for dihType = settings.dihedralTypes
%             finalNdx = [finalNdx ];
% 
%         end

        kl1Res{thisSys}(i) = sum(kl1{thisSys}(ndxHere(1:min(length(ndxHere) ,settings.dihedralTypes ))));
        
    end
end

