%% Prepare input files:

% % Father directory for all the runs
% databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
% metadir = 'D:\Telework_library\dopamine_phase_3';
% 
% % Assumes 1st directory is reference in alignment and naming!
% foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT','44-D2_WT_DPA_Gi_7JVR', ...
%     '42-D2_WT_BRC_Gi_7JVR', '34-d2_dpa_barr2_WT_for_Real','35-d2_brc_barr2_WT_forREAL', ...
%     '21-d2_ris_6cm4'}; 
% % 
% % foldersToStudy = {'1-d2_dop_WT','2-d2_dop_T174M-C220L','10-d2_dop_L214M','12-d2_dop_L92G','14-d2_dop_F217M' ...
% %     '3-d2_bromo_WT','4-d2_bromo_T174M-C220L', ...
% %     '11-d2_brc_L214M','13-d2_brc_L92G','15-d2_brc_F217M'};
% 
% % foldersToStudy = {'3-d2_bromo_WT','4-d2_bromo_T174M-C220L','9-d2_bromo_V130F', ...
% %     '11-d2_brc_L214M','13-d2_brc_L92G','15-d2_brc_F217M'};
% 
% thisSysLabel = {'DA 6VMS','BRC 6VMS','DA 7JVR', ...
%     'BRC 7JVR', 'DA barr2 AF2','BRC barr2 AF2','RIS 6CM4'};
% % thisSysLabel = {'DA WT','DA T5.54M-C6.47L','DA L6.41M','DA L3.41G','DA F6.44M', ...
% %     'BRC WT','BRC T5.54M-C6.47L','BRC L6.41M','BRC L3.41G','BRC F6.44M' };
% name = 'DD2R_WT_extended';
% 
% Generic GPCR numbering from https://gpcrdb.org/residue/residuetable
gpcrdbRefName = 'DD2R';

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
chains = 'A-B';

% Reference PDB code that this simulation is based on
pdbCode = '6VMS';
pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
pdbCodeInactive = '6CM4';
pdbInactiveChains = 'A--';

isGPCR = true;
useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment
saveData = true;
plotrmsd = false ;
%% Make the meta-direcory
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end
md2pathdirMeta = [metadir slash 'md2pathMeta' slash];
md2pathdir = [metadir slash 'md2pathMeta' slash name slash];

% Create meta directory
if ~exist(md2pathdirMeta, 'dir')
    mkdir(md2pathdirMeta)
    add2log(md2pathdirMeta, ['Creating md2path directory in ' metadir]);
end

% Create subdirectories for every bundle of systems studied
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir)
    add2log(md2pathdir, ['Creating md2path directory in ' md2pathdirMeta]);
end

%% Load, fetch and align PDB files

metaLoadAlign

%% Load contact maps from separate folders

% I take the maximum contact from every contact map, which is not STRICTLY
% CORRECT, since it may miss some overlapping contacts between different
% atoms from the ligand to the same protein residue 

contactCutMeta = 0.3; % Cutoff for signficant contacts, ALL rows should 
% show contact less than cutoff or else it won't zero them out

meta_contacts = [];

simResList =  database.residues{1}{:,1:length(foldersToStudy)}; % list from simulation set
% Find rows that are perfectly aligned (no zeros!)
allAligned = sum(simResList==0,2)==0;

for thisSys = 1:length(foldersToStudy)
    L = load([metadir slash foldersToStudy{thisSys} slash 'md2pathdev' slash  'workspace.mat']);
    contactVector = max(mean(L.contactsPerRun,3));

    if useDatabaseAlignment
        meta_contacts = [meta_contacts ; contactVector(simResList(allAligned,thisSys))]; 
    else
        meta_contacts = [meta_contacts ; contactVector]; 
    end
end


% meta_contacts(meta_contacts<contactCutMeta)=0; % Zero non-significant contacts
toCut = sum(meta_contacts<contactCutMeta); % If ALL variants have less than 
% cutoff contact frequency, then I'll zero the column

for i = 1:length(toCut)
   if toCut(i) == length(foldersToStudy)
       meta_contacts(:,i) = zeros(length(foldersToStudy),1);
   end
end

[row, col] = find(meta_contacts); % Find nonzero contacts over all runs (averaged)
contact_map = meta_contacts(unique(row),unique(col));

if useDatabaseAlignment
    receptorResIds = simResList(allAligned,refNdx);
end

liProtRes = receptorResIds(unique(col)); % Contacting receptor residues

%% Plot the figure:

% Grab BW numbering of ligand binding residues from database/chain
resname_num = refChain.formatResidues(liProtRes,'SecondaryEntry',database.entries{length(foldersToStudy)+1},'BWonly', true);
X = round(contact_map,2,"significant");
figure
h = heatmap(resname_num,thisSysLabel,X,'Colormap',parula)

% xlabel('Protein residue', 'FontSize', fontsz)
h.Title = 'Ligand contact frequency from simulations';
h.XLabel = 'Protein residue';
h.YLabel = 'Ligand';
set(gca,'FontSize',12)

%% Save!!
if saveData
savefig([md2pathdir 'contact_map_' name]);
print2pdf([md2pathdir 'contact_map_' name]);

% Save table too:
tabTemp = array2table(contact_map);
tabTemp.Properties.VariableNames = resname_num;
catTable = [table(thisSysLabel') tabTemp]; % Meow
writetable(catTable, [md2pathdir 'contact_map_' name '.xlsx'])
end
%% Why not do RMSD too:
if plotrmsd
legend_count = 1;
errStrideRMSD = 50;
figure
for thisSys = 1:length(foldersToStudy)
    L = load([metadir slash foldersToStudy{thisSys} slash 'md2pathdev' slash  'workspace.mat']);
    
    RMSD_mean_L= mean(L.RMSD_L,2);
    RMSD_std_L= std(L.RMSD_L,0,2);
    plot(RMSD_mean_L,'LineWidth',1);
    hold on
    errorbar(1:errStrideRMSD:length(RMSD_mean_L),RMSD_mean_L(1:errStrideRMSD:length(RMSD_mean_L)),...
        RMSD_std_L(1:errStrideRMSD:length(RMSD_mean_L))/sqrt( size(L.RMSD,2)),'o','MarkerSize',5)
    legend_entries{legend_count} = thisSysLabel{thisSys};
    legend_count = legend_count + 1;
    legend_entries{legend_count} = '';
    legend_count = legend_count + 1;
end

% Add annotations
xlabel('Frame'); ylabel('RMSD [Angstrom]')
legend(legend_entries,'location','best');
legend boxoff
title('RMSD of non-H ligand atoms')
set(gca,'FontSize',16)


%% Save!!

savefig([md2pathdir 'ligand_RMSD_' name]);
print2pdf([md2pathdir 'ligand_RMSD_' name]);
end