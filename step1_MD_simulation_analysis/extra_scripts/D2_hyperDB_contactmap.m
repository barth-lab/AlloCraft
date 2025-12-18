% all data load:
%% Experimental data:
tabExp = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data\experimental_data_labeled.xlsx','Sheet','D2');

tabExpDeltas = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data\delta_experimental_data.xlsx','Sheet','D2');
% Load colors
paperColors; % Loads paperColorPalette

% KLDiv: 
%% For metaKLDiv:

% Reference PDB code that this simulation is based on
settings.pdbCode = '6VMS';
settings.pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
settings.pdbCodeInactive = '6CM4';
settings.pdbInactiveChains = 'A--';

settings.isGPCR = true;
settings.includeLigPathwayCalc = false;

% Which types of dihedrals to be considered for KL1 calculation:
% phi = 1; psi = 2; chi1 = 3; chi2 = 4; ... chi5 = 7
% settings.dihedralTypes = [1 2 3 4 5 6 7]; 
settings.dihedralTypes = 7;
% metaKLDiv: meta-analysis of KLDiv from several directories

% Input section
settings.databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
% Father directory for all the runs
settings.metadir = 'D:\Telework_library\dopamine_phase_3';
settings.md2pathName = 'md2pathdev';

%% (almost) full set for PCA: RIS inactive reference
settings.foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT', '2-d2_dop_T174M-C220L', '4-d2_bromo_T174M-C220L', ...
     '10-d2_dop_L214M','11-d2_brc_L214M','12-d2_dop_L92G', '13-d2_brc_L92G', ...
      '14-d2_dop_F217M', '15-d2_brc_F217M',  '16-d2_dop_I125N', '17-d2_brc_I125N','18a-d2_dop_F217I', '18b-d2_brc_F217I', ...
      '19a-d2_dop_L92H','19b-d2_brc_L92H', '34-d2_dpa_barr2_WT_for_Real', '35-d2_brc_barr2_WT_forREAL',...
       '31-d2_dpa_barr2', '32-d2_brc_barr2', '36b-d2_dpa_barr2_F217M', '36a-d2_brc_barr2_F217M', '37b-d2_dpa_barr2_I125N', ...
      '37a-d2_brc_barr2_I125N', '38b-d2_dpa_barr2_L214M','38a-d2_brc_barr2_L214M','39b-d2_dpa_barr2_L92H', ...
      '51b-d2_dpa_barr2_T174M-C220L'}; %, ...
      %'42-D2_WT_BRC_Gi_7JVR','44-D2_WT_DPA_Gi_7JVR'};

settings.thisSysLabel = {'DA-WT-Gi','BRC-WT-Gi','DA-T5.54M-C6.47L-Gi','BRC-T5.54M-C6.47L-Gi', ...
    'DA-L6.41M-Gi','BRC-L6.41M-Gi','DA-L3.41G-Gi','BRC-L3.41G-Gi', ...
    'DA-F6.44M-Gi','BRC-F6.44M-Gi','DA-I4.46N-Gi','BRC-I4.46N-Gi','DA-F6.44I-Gi','BRC-F6.44I-Gi', ...
    'DA-L3.41H-Gi','BRC-L3.41H-Gi','DA-WT-Barr2', 'BRC-WT-Barr2','DA-L3.41G-Barr2', ...
    'BRC-L3.41G-Barr2','DA-F6.44M-Barr2','BRC-F6.44M-Barr2','DA-I4.46N-Barr2', ...
    'BRC-I4.46N-Barr2','DA-L6.41M-Barr2','BRC-L6.41M-Barr2','DA-L3.41H-Barr2', ...
    'DA-T5.54M-C6.47L-Barr2'}; %, ...
%     '7JVR-DA-WT-Gi','7JVR-BRC-WT-Gi'};


settings.name = 'DD2R-June_29_2023';

settings.gpcrdbRefName = 'DD2R';

% Reference system for the KLDiv calculations
settings.refdir =  'D:\Telework_library\dopamine_phase_3\21-d2_ris_6cm4';
settings.refName = 'RIS-inactive'; % MAKE SURE THIS IS THE SAME NAME AS WHEN YOU CALCULATED THE KLDIVs!

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'A-B';

% Run:
metaKLDiv;
kl1RisInactive = kl1Res;
kl1RisInactiveSysLabel = settings.thisSysLabel;

%% WT reference (for every system, its equivalent WT) 
% so for example: DA-F6.44M-Gi will have DA-WT-Gi as its reference WT

% (almost) full set for PCA:
settings.foldersToStudy = { '2-d2_dop_T174M-C220L', '4-d2_bromo_T174M-C220L', ...
     '10-d2_dop_L214M','11-d2_brc_L214M','12-d2_dop_L92G', '13-d2_brc_L92G', ...
      '14-d2_dop_F217M', '15-d2_brc_F217M',  '16-d2_dop_I125N', '17-d2_brc_I125N','18a-d2_dop_F217I', '18b-d2_brc_F217I', ...
      '19a-d2_dop_L92H','19b-d2_brc_L92H', '31-d2_dpa_barr2', '32-d2_brc_barr2', ...
       '36a-d2_brc_barr2_F217M', '36b-d2_dpa_barr2_F217M', '37a-d2_brc_barr2_I125N', ...
      '38a-d2_brc_barr2_L214M','37b-d2_dpa_barr2_I125N', '38b-d2_dpa_barr2_L214M','39b-d2_dpa_barr2_L92H', ...
      '51b-d2_dpa_barr2_T174M-C220L'};

settings.thisSysLabel = {'DA-T5.54M-C6.47L-Gi','BRC-T5.54M-C6.47L-Gi', ...
    'DA-L6.41M-Gi','BRC-L6.41M-Gi','DA-L3.41G-Gi','BRC-L3.41G-Gi', ...
    'DA-F6.44M-Gi','BRC-F6.44M-Gi','DA-I4.46N-Gi','BRC-I4.46N-Gi','DA-F6.44I-Gi','BRC-F6.44I-Gi', ...
    'DA-L3.41H-Gi','BRC-L3.41H-Gi','DA-L3.41G-Barr2','BRC-L3.41G-Barr2', ...
    'BRC-F6.44M-Barr2','DA-F6.44M-Barr2','BRC-I4.46N-Barr2', ...
    'BRC-L6.41M-Barr2','DA-I4.46N-Barr2','DA-L6.41M-Barr2','DA-L3.41H-Barr2', ...
    'DA-T5.54M-C6.47L-Barr2'};


settings.name = 'DD2R-June_29_2023';

settings.gpcrdbRefName = 'DD2R';

% Reference system for the KLDiv calculations
settings.refdir =  'D:\Telework_library\dopamine_phase_3\21-d2_ris_6cm4';
settings.refName = 'WT'; % MAKE SURE THIS IS THE SAME NAME AS WHEN YOU CALCULATED THE KLDIVs!

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
settings.chains = 'A-B';

% Run:
metaKLDiv;
kl1WT = kl1Res;
kl1WTSysLabel = settings.thisSysLabel;

% MI
%% MI per residue
metadir = 'D:\Telework_library\dopamine_phase_3';
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
pdbName = 'prot.pdb';

foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT', '2-d2_dop_T174M-C220L', '4-d2_bromo_T174M-C220L', ...
     '10-d2_dop_L214M','11-d2_brc_L214M','12-d2_dop_L92G', '13-d2_brc_L92G', ...
      '14-d2_dop_F217M', '15-d2_brc_F217M',  '16-d2_dop_I125N', '17-d2_brc_I125N','18a-d2_dop_F217I', '18b-d2_brc_F217I', ...
      '19a-d2_dop_L92H','19b-d2_brc_L92H', '34-d2_dpa_barr2_WT_for_Real', '35-d2_brc_barr2_WT_forREAL',...
       '31-d2_dpa_barr2', '32-d2_brc_barr2', '36b-d2_dpa_barr2_F217M', '36a-d2_brc_barr2_F217M', '37b-d2_dpa_barr2_I125N', ...
      '37a-d2_brc_barr2_I125N', '38b-d2_dpa_barr2_L214M','38a-d2_brc_barr2_L214M','39b-d2_dpa_barr2_L92H', ...
      '51b-d2_dpa_barr2_T174M-C220L','21-d2_ris_6cm4'};

thisSysLabel = {'DA-WT-Gi','BRC-WT-Gi','DA-T5.54M-C6.47L-Gi','BRC-T5.54M-C6.47L-Gi', ...
    'DA-L6.41M-Gi','BRC-L6.41M-Gi','DA-L3.41G-Gi','BRC-L3.41G-Gi', ...
    'DA-F6.44M-Gi','BRC-F6.44M-Gi','DA-I4.46N-Gi','BRC-I4.46N-Gi','DA-F6.44I-Gi','BRC-F6.44I-Gi', ...
    'DA-L3.41H-Gi','BRC-L3.41H-Gi','DA-WT-Barr2', 'BRC-WT-Barr2','DA-L3.41G-Barr2', ...
    'BRC-L3.41G-Barr2','DA-F6.44M-Barr2','BRC-F6.44M-Barr2','DA-I4.46N-Barr2', ...
    'BRC-I4.46N-Barr2','DA-L6.41M-Barr2','BRC-L6.41M-Barr2','DA-L3.41H-Barr2', ...
    'DA-T5.54M-C6.47L-Barr2','RIS-inactive'};

md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc_Culled_data';

MIPerRes = [];

for i=1:length(foldersToStudy)
    dirHere = fullfile(metadir,foldersToStudy{i},md2pathName,alloPathCalcName,'MIperresidue.txt');
    if exist(dirHere,'file')
        temp = importdata(dirHere);
        MIPerRes = [MIPerRes temp];
    end
end
MIPerResSysLabel = thisSysLabel;
%% MIDomains


% Hubscore related stuff:
%% All hubscores:
tabHubscores = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Compiled_AlloDy_data\hubscores_compiled_latest.xlsx','Sheet','D2');
hubscoreLabel = {'DA-WT-Gi','BRC-WT-Gi','DA-T5.54M-C6.47L-Gi','BRC-T5.54M-C6.47L-Gi', ...
    'DA-L6.41M-Gi','BRC-L6.41M-Gi','DA-L3.41G-Gi','BRC-L3.41G-Gi', 'DA-F6.44M-Gi','BRC-F6.44M-Gi', ...
    'DA-V4.51F-Gi','BRC-V4.51F-Gi','DA-I4.46N-Gi','BRC-I4.46N-Gi','DA-F6.44I-Gi','BRC-F6.44I-Gi', ...
    'DA-L3.41H-Gi','BRC-L3.41H-Gi','DA-WT-Barr2', 'BRC-WT-Barr2','DA-L3.41G-Barr2', ...
    'BRC-L3.41G-Barr2','DA-F6.44M-Barr2','BRC-F6.44M-Barr2','DA-I4.46N-Barr2', ...
    'BRC-I4.46N-Barr2','DA-L6.41M-Barr2','BRC-L6.41M-Barr2','DA-L3.41H-Barr2', ...
    'DA-T5.54M-C6.47L-Barr2','RIS-inactive'};

%% Hubscores of conserved residues:
tabConserved = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Compiled_AlloDy_data\conserved_motifs_hubscores_latest.xlsx','Sheet','D2');
% same labels as tabHubscores


% Contact map:
%% MetaContactMap:
useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment
name = 'DD2R_WT_extended';

% Gi only
foldersToStudy = {'1-d2_dop_WT', '2-d2_dop_T174M-C220L','10-d2_dop_L214M', ...
    '12-d2_dop_L92G','14-d2_dop_F217M', '16-d2_dop_I125N','18a-d2_dop_F217I', '19a-d2_dop_L92H',  ...
      '3-d2_bromo_WT', '4-d2_bromo_T174M-C220L', '11-d2_brc_L214M', '13-d2_brc_L92G', ...
       '15-d2_brc_F217M',   '17-d2_brc_I125N', '18b-d2_brc_F217I', ...
      '19b-d2_brc_L92H','21-d2_ris_6cm4'};

thisSysLabel = {'DA-WT-Gi','DA-T5.54M-C6.47L-Gi','DA-L6.41M-Gi','DA-L3.41G-Gi', 'DA-F6.44M-Gi', ...
    'DA-I4.46N-Gi','DA-F6.44I-Gi', 'DA-L3.41H-Gi', 'BRC-WT-Gi','BRC-T5.54M-C6.47L-Gi',  ...
    'BRC-L6.41M-Gi','BRC-L3.41G-Gi','BRC-F6.44M-Gi','BRC-I4.46N-Gi','BRC-F6.44I-Gi', ...
    'BRC-L3.41H-Gi','RIS-inactive'};
metaContactMap; % Outputs contact_map and liProtRes

contactLabel = thisSysLabel;


%% More contact stuff: sum of atom-wise contact with residues:

useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment
name = 'DD2R_WT_extended';

% Gi only
% foldersToStudy = {'1-d2_dop_WT', '2-d2_dop_T174M-C220L','10-d2_dop_L214M', ...
%     '12-d2_dop_L92G','14-d2_dop_F217M', '16-d2_dop_I125N','18a-d2_dop_F217I', '19a-d2_dop_L92H',  ...
%     };
foldersToStudy = {      '3-d2_bromo_WT', '4-d2_bromo_T174M-C220L', '11-d2_brc_L214M', '13-d2_brc_L92G', ...
       '15-d2_brc_F217M',   '17-d2_brc_I125N', '18b-d2_brc_F217I', ...
      '19b-d2_brc_L92H'};

% thisSysLabel = {'DA-WT-Gi','DA-T5.54M-C6.47L-Gi','DA-L6.41M-Gi','DA-L3.41G-Gi', 'DA-F6.44M-Gi', ...
%     'DA-I4.46N-Gi','DA-F6.44I-Gi', 'DA-L3.41H-Gi'};
thisSysLabel = {'BRC-WT-Gi','BRC-T5.54M-C6.47L-Gi',  ...
    'BRC-L6.41M-Gi','BRC-L3.41G-Gi','BRC-F6.44M-Gi','BRC-I4.46N-Gi','BRC-F6.44I-Gi', ...
    'BRC-L3.41H-Gi'};

meta_contacts = [];
meta_contacts_err = [];
simResList =  database.residues{1}{:,1:length(foldersToStudy)}; % list from simulation set
% Find rows that are perfectly aligned (no zeros!)
allAligned = sum(simResList==0,2)==0;

for thisSys = 1:length(foldersToStudy)
    L = load([metadir slash foldersToStudy{thisSys} slash 'md2pathdev' slash  'workspace.mat']);
    contactVector = (mean(L.contactsPerRun,3));
    contactVectorErr = (std(L.contactsPerRun,[],3))/sqrt(size(L.contactsPerRun,3));

    meta_contacts = [meta_contacts sum(contactVector,2)];
    meta_contacts_err= [meta_contacts_err arrayfun(@(x) sqrt(sumsqr(contactVectorErr(x,:))),1:size(contactVector,1))'];
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
% contact_map = meta_contacts(unique(row),unique(col));

if useDatabaseAlignment
    receptorResIds = simResList(allAligned,refNdx);
end

liProtRes = receptorResIds(unique(col)); % Contacting receptor residues

% Get ligand atom names:
ligandAtomNames = database.entries{1}.pdb.name(ligandChain.getLigandAtoms,:);
ligandAtomContactLabel = thisSysLabel;
%% Combining data from different entries: plot heatmap with hubscores/MI for ligand binding residues:
hdb = HyperDataBase;
hdb.read(contact_map', contactLabel,'Contact Map','VarNames',resname_num)
hdb.read(MIPerRes,MIPerResSysLabel,'MI per residue');
hdb.read(tabHubscores{:,3:end-2},hubscoreLabel,'Hubscores','VarNames',tabHubscores.(2));
% hdb.read(interTMDomainMISumCell', settings.thisSysLabel,'DomainMI','VarNames',[tmDomain.name])
hdb.read(cell2mat(kl1RisInactive'),kl1RisInactiveSysLabel,'KLDiv wrt RIS');
hdb.read(cell2mat(kl1WT'),kl1WTSysLabel,'KLDiv wrt WT');
% hdb.read([tabExp.Efficacy tabExp.LogEC50]',tabExp.Label,'Experimental data');
% hdb.read([tabExpDeltas.dEfficacy tabExpDeltas.dLogEC50]',tabExpDeltas.Label, 'Experiment d(Mut-WT)','VarNames',["dEfficacy" "dlog(EC50)"]);
hdb.read([tabExp.Efficacy tabExp.LogEC50 tabExp.Efficacy_SEM tabExp.LogEC50_SEM]',tabExp.Label,'Experimental data', ...
    'VarNames',["Efficacy" "log(EC50)" "Efficacy_SEM" "log(EC50)_SEM"]);
dexpData = [tabExpDeltas.dEfficacy tabExpDeltas.dLogEC50 tabExpDeltas.dEfficacy_normalized ...
    tabExpDeltas.dEfficacy_SEM tabExpDeltas.dLogEC50_SEM tabExpDeltas.dEfficacy_normalized_SEM]';
hdb.read(dexpData,tabExpDeltas.Label, 'Experiment d(Mut-WT)', ...
    'VarNames',["dEfficacy" "dlog(EC50)" "dEfficacy_norm" "dEfficacy_SEM" "dlog(EC50)_SEM" "dEfficacy_norm_SEM"]);

hdb.read(meta_contacts, ligandAtomContactLabel,'Contact Sum All Res','VarNames',ligandAtomNames)
hdb.read(meta_contacts_err, ligandAtomContactLabel,'Contact Sum All ERR','VarNames',ligandAtomNames)

hdb.align;
hdb.alignment

 % Run metaContactMap first
hDBNdx = 3;
contactCutMeta;
meta_contacts;
liProtRes;
alignHere = hdb.alignment.(string(hdb.entryName(hDBNdx)));

dataHere = hdb.hyperEntries{hDBNdx}.rawData(liProtRes,alignHere)';
dataHere(contact_map<contactCutMeta) = nan; % Remove data with low contact

X = round(dataHere,3,"significant");
figure
h = heatmap(hdb.hyperEntries{1}.VarNames,hdb.hyperEntries{1}.labels,X,'Colormap',turbo);

% xlabel('Protein residue', 'FontSize', fontsz)
h.Title =  hdb.hyperEntries{hDBNdx}.name;
h.XLabel = 'Protein residue';
h.YLabel = 'Ligand';
h.MissingDataColor=[0.8 0.8 0.8];
set(gca,'FontSize',12)

%% Plot delta-contacts from respective WTs
ligNames = [ "DA" "BRC"];

ddataContact = [];
labelsAll = [];
for i =1:length(ligNames)
    dataContact = hdb.hyperEntries{1}.fetchData(ligNames(i));
    ddataContactHere = [dataContact(:,1) (dataContact(:,2:end) - dataContact(:,1))];
    labelsHere = hdb.hyperEntries{1}.labels(hdb.hyperEntries{1}.fetchLabels(ligNames(i)));
    ddataContact = [ddataContact ddataContactHere];
    labelsAll = [labelsAll; labelsHere];
end

ddataContact(abs(ddataContact)<0.01) = 0; % zero out very small differences 
X = round(ddataContact',2,"significant");
figure
h = heatmap(hdb.hyperEntries{1}.VarNames,labelsAll,X,'Colormap',myColorLavenderGrey);

caxis([-max(abs(caxis)) max(abs(caxis))]); % Make colorbar symmetrical
% xlabel('Protein residue', 'FontSize', fontsz)
h.Title =  '\Delta Contact from WT';
h.XLabel = 'Protein residue';
h.YLabel = 'Ligand';
h.MissingDataColor=[0.8 0.8 0.8];
set(gca,'FontSize',12)


%% Correlate ligand contact sum with potency
ligandName = "BRC"; % "DA" or "BRC"
propertyName = "Potency"; % "Efficacy" or "Potency"
if strcmp(propertyName,"Efficacy")
    dataNdx = 1;
elseif strcmp(propertyName,"Potency")
    dataNdx = 2;
end
hdb.align("refEntry",8);
hdb.alignment
dataTable = hdb.dbFetchData([ligandName "Gi"],"hyperEntryNdx",[6 8 9],'refEntry',8);
labelNamesShort = extractBetween(dataTable.Label,ligandName + "-","-Gi");


x = sum(dataTable.("Contact Sum All Res"),2);
y = dataTable.("Experimental data")(:,dataNdx);
figure
s = scatter(x,y,50,'MarkerFaceColor',paperColorPalette.Hex(2),'MarkerEdgeColor',paperColorPalette.Hex(1));

% Add data labels:
row = dataTipTextRow('Label',labelNamesShort);
s.DataTipTemplate.DataTipRows(end+1) = row;
xRange = range(x);
yRange = range(y);
text(x+xRange/100,y+yRange/50, labelNamesShort,'FontSize',10);  

hold on

% Add error bars:
yneg = dataTable.("Experimental data")(:,dataNdx+2);
ypos = yneg;
xneg = mean(dataTable.("Contact Sum All ERR"),2);
xpos = xneg;
errorbar(x,y,yneg,ypos,xneg,xpos,"LineStyle","none","Marker","none","Color",paperColorPalette.Hex(1))
hold on

% Fitting: 
if strcmp(propertyName,"Efficacy")
    y(isnan(y)) = 0;
elseif strcmp(propertyName,"Potency")
    y(isnan(y)) = mean(y,"omitnan");
end
X = [ones(length(x),1) x];
r = X\y;
yCalc2 = X*r;
plot(x,yCalc2,'-','Color',paperColorPalette.Hex(3),'LineWidth',1.5)
Rsq = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
[~,p] = corrcoef(x,y);
p = p(1,2);
legend('Data', '' ,['Fit, R^2=' num2str(Rsq,2) ', p=' num2str(p,2)],'Location','best');

xlabel('Summed contact')
ylabel(sprintf('Variant %s',propertyName))
title(sprintf('%s %s',ligandName, propertyName))
formatplot2