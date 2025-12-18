% metaLoadAlign: loads PDBs from all studied systems and reference and then
% align and add ref generic numbering where applicable:

%% Load, fetch and align PDB files
database = Database(databasePath);

% Read PDB of every system to be studied:
for thisSys = 1:length(foldersToStudy)
    % Get local directory
%     mydir = ([metadir slash foldersToStudy{thisSys}]); 
    mydir = fullfile(metadir,foldersToStudy{thisSys});
    database.read(fullfile(mydir, "prot.pdb"), chains, thisSysLabel{thisSys});

end
if ~exist('refNdx','var')
    refNdx = 1; % This entry spot usually reserved for active ref pdb
end
refEntry = database.entries{refNdx}; 

% Store names for each chain
Chains.receptor = 1;
Chains.gprotein = 2;
Chains.ligand = 3;

refChain = refEntry.chains{Chains.receptor};

if refEntry.hasChain(Chains.ligand)
    ligandChain = refEntry.chains{Chains.ligand};
end

% Get list of residues
receptorResIds = refChain.resIds;

% Read reference PDB structures
if ~isempty(pdbCode)
    database.fetch(pdbCode, pdbChains);

    if ~isempty(pdbCodeInactive)
        database.fetch(pdbCodeInactive, pdbInactiveChains);
    end
end

% Align sequences and produce the database.residues table
database.align();

% Align structures to the first database entry
database.alignStructures();

% Add labels with Ballesteros-Weinstein notation, aligned with the ref.
% database entry from the PDB
% Exported from GPCRdb (https://gpcrdb.org/residue/residuetable)
if length(database.entries) > (length(foldersToStudy)+1) && isGPCR
    residueTablePath = fullfile(database.dir, gpcrdbRefName + "_residue_table.xlsx");
    database.label(length(foldersToStudy)+1, residueTablePath); 
end