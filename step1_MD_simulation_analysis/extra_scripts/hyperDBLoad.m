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

%% (almost) full set for PCA: RIS inactive reference
settings.foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT', '2-d2_dop_T174M-C220L', '4-d2_bromo_T174M-C220L', ...
     '10-d2_dop_L214M','11-d2_brc_L214M','12-d2_dop_L92G', '13-d2_brc_L92G', ...
      '14-d2_dop_F217M', '15-d2_brc_F217M',  '16-d2_dop_I125N', '17-d2_brc_I125N', '18a-d2_dop_F217I', '18b-d2_brc_F217I', ...
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


settings.name = 'DD2R-Feb-1-2024';
settings.md2pathName = 'md2pathdev';
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
      '14-d2_dop_F217M', '15-d2_brc_F217M',  '16-d2_dop_I125N', '17-d2_brc_I125N', '18a-d2_dop_F217I', '18b-d2_brc_F217I', ...
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


% settings.name = 'DD2R-June_29_2023';

% settings.gpcrdbRefName = 'DD2R';

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
      '14-d2_dop_F217M', '15-d2_brc_F217M','8-d2_dop_V130F','9-d2_bromo_V130F', ...
      '16-d2_dop_I125N', '17-d2_brc_I125N','18a-d2_dop_F217I', '18b-d2_brc_F217I', ...
      '19a-d2_dop_L92H','19b-d2_brc_L92H', '34-d2_dpa_barr2_WT_for_Real', '35-d2_brc_barr2_WT_forREAL',...
       '31-d2_dpa_barr2', '32-d2_brc_barr2', '36b-d2_dpa_barr2_F217M', '36a-d2_brc_barr2_F217M', '37b-d2_dpa_barr2_I125N', ...
      '37a-d2_brc_barr2_I125N', '38b-d2_dpa_barr2_L214M','38a-d2_brc_barr2_L214M','39b-d2_dpa_barr2_L92H', ...
      '51b-d2_dpa_barr2_T174M-C220L','21-d2_ris_6cm4'};

thisSysLabel = {'DA-WT-Gi','BRC-WT-Gi','DA-T5.54M-C6.47L-Gi','BRC-T5.54M-C6.47L-Gi', ...
    'DA-L6.41M-Gi','BRC-L6.41M-Gi','DA-L3.41G-Gi','BRC-L3.41G-Gi', 'DA-F6.44M-Gi','BRC-F6.44M-Gi', ...
    'DA-V4.51F-Gi','BRC-V4.51F-Gi','DA-I4.46N-Gi','BRC-I4.46N-Gi','DA-F6.44I-Gi','BRC-F6.44I-Gi', ...
    'DA-L3.41H-Gi','BRC-L3.41H-Gi','DA-WT-Barr2', 'BRC-WT-Barr2','DA-L3.41G-Barr2', ...
    'BRC-L3.41G-Barr2','DA-F6.44M-Barr2','BRC-F6.44M-Barr2','DA-I4.46N-Barr2', ...
    'BRC-I4.46N-Barr2','DA-L6.41M-Barr2','BRC-L6.41M-Barr2','DA-L3.41H-Barr2', ...
    'DA-T5.54M-C6.47L-Barr2','RIS-inactive'};

md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc_Culled_data';

MIPerRes = [];
MIresAll = cell(length(foldersToStudy),1);
for i=1:length(foldersToStudy)
    % Averaged MI per residue
%     dirHere = fullfile(metadir,foldersToStudy{i},md2pathName,alloPathCalcName,'MIperresidue.txt'); 
%     if exist(dirHere,'file')
%         temp = importdata(dirHere);
%         MIPerRes = [MIPerRes temp];
%     end

    % Summed MI per residue:
    dirHere = fullfile(metadir,foldersToStudy{i},md2pathName,alloPathCalcName,'workspace.mat');
    temp = load(dirHere,'-mat','MIres');
    MIres = temp.MIres;
    
    if exist(dirHere,'file')
        temp = load(dirHere,'-mat','MIres');
        MIres = temp.MIres;
        MIPerResHere = sum(MIres,2);
        MIPerRes = [MIPerRes MIPerResHere];
        MIresAll{i} = MIres;
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
name = 'temp';
metaContactMap; % Outputs contact_map and liProtRes

contactLabel = thisSysLabel;
%% Read the data into the hyperDB
hdb = HyperDataBase;
% hdb.read([tabExp.Efficacy tabExp.LogEC50]',tabExp.Label,'Experimental data');
% hdb.read([tabExpDeltas.dEfficacy tabExpDeltas.dLogEC50 tabExpDeltas.dEfficacy_normalized]',tabExpDeltas.Label, 'Experiment d(Mut-WT)','VarNames',["dEfficacy" "dlog(EC50)" "dEfficacy_norm"]);
hdb.read([tabExp.Efficacy tabExp.LogEC50 tabExp.Efficacy_SEM tabExp.LogEC50_SEM]',tabExp.Label,'Experimental data', ...
    'VarNames',["Efficacy" "log(EC50)" "Efficacy_SEM" "log(EC50)_SEM"]);
dexpData = [tabExpDeltas.dEfficacy tabExpDeltas.dLogEC50 tabExpDeltas.dEfficacy_normalized ...
    tabExpDeltas.dEfficacy_SEM tabExpDeltas.dLogEC50_SEM tabExpDeltas.dEfficacy_normalized_SEM]';
hdb.read(dexpData,tabExpDeltas.Label, 'Experiment d(Mut-WT)', ...
    'VarNames',["dEfficacy" "dlog(EC50)" "dEfficacy_norm" "dEfficacy_SEM" "dlog(EC50)_SEM" "dEfficacy_norm_SEM"]);


% hdb.read(contact_map', contactLabel,'Contact Map','VarNames',resname_num)
hdb.read(MIPerRes,MIPerResSysLabel,'MI per residue','VarNames',tabHubscores.(2));
hdb.read(MIresAll',MIPerResSysLabel,'MI','VarNames',tabHubscores.(2));
hdb.read(tabHubscores{:,3:end-2},hubscoreLabel,'Hubscores','VarNames',tabHubscores.(2));
% hdb.read(tabConserved{:,2:end},hubscoreLabel,'Conserved motifs HS','VarNames',tabConserved.Motif);

% hdb.read(interTMDomainMISumCell', settings.thisSysLabel,'DomainMI','VarNames',[tmDomain.name])
hdb.read(cell2mat(kl1RisInactive'),kl1RisInactiveSysLabel,'KLDiv wrt RIS');
hdb.read(cell2mat(kl1WT'),kl1WTSysLabel,'KLDiv wrt WT');

hdb.align;
hdb.alignment

%% Example of fetching data

dataTable = hdb.dbFetchData(["DA" "Gi"],"hyperEntryNdx",[1 2 3]);
%% Using labels and alignments: calculating delta delta efficacy (BRC - DA)(mutant - WT)

mutList  = unique(tabExp.Mutation,'stable');
mutList = mutList(1:end-1); % Remove C6.47L
mutList(3) = []; % remove V4.51F

deltadelta = zeros(length(mutList),4);
var2plot = [2 3 5 6]; %[2 1 5 4] to plot raw efficacies and [2 3 5 6] for normalized efficacies
for i = 1:length(mutList)
    dataTable = hdb.dbFetchData([string(mutList(i)) "Gi"],"hyperEntryNdx",[1 2]);
    if size(dataTable,1) > 1
    deltadelta(i,:) = dataTable.("Experiment d(Mut-WT)")(2,var2plot) - dataTable.("Experiment d(Mut-WT)")(1,var2plot);
    end
end

deltadelta(isnan(deltadelta))=0;
table(mutList,deltadelta)

%
pos = [112.2000  109.0000  772.8000  640.0000];
figure('Position',pos)
xRange = range(deltadelta(:,1));
yRange = range(deltadelta(:,2));
s =scatter(deltadelta(:,1),deltadelta(:,2),100,'LineWidth',1.5,'MarkerFaceColor',paperColorPalette.Hex(8), ...
   'MarkerEdgeColor',paperColorPalette.Hex(7));
hold on
xlabel('ddLog(EC50)')
ylabel('ddEfficacy')
row = dataTipTextRow('Label',mutList);
s.DataTipTemplate.DataTipRows(end+1) = row;

ytextPos = deltadelta(:,2)+yRange/25;
ytextPos(4:5) = deltadelta(4:5,2)-yRange/25; % Fix text for 3.41H and 3.41G

text(deltadelta(:,1)+xRange/40,ytextPos, mutList,'FontSize',14);     
% Add error bars:
yneg = deltadelta(:,4);
ypos = yneg;
xneg = deltadelta(:,3);
xpos = xneg;
errorbar(deltadelta(:,1),deltadelta(:,2),yneg,ypos,xneg,xpos,"LineStyle","none","Marker","none","Color",paperColorPalette.Hex(7))

title('\Delta\Delta (BRC - DA) (Mut - WT)')
formatplot2
%% Combining data from different entries: plot heatmap with hubscores/MI for ligand binding residues:
hdb = HyperDataBase;
hdb.read(contact_map', contactLabel,'Contact Map','VarNames',resname_num)
hdb.read(MIPerRes,MIPerResSysLabel,'MI per residue');
hdb.read(tabHubscores{:,3:end-2},hubscoreLabel,'Hubscores','VarNames',tabHubscores.(2));
% hdb.read(interTMDomainMISumCell', settings.thisSysLabel,'DomainMI','VarNames',[tmDomain.name])
hdb.read(cell2mat(kl1RisInactive'),kl1RisInactiveSysLabel,'KLDiv wrt RIS');
hdb.read(cell2mat(kl1WT'),kl1WTSysLabel,'KLDiv wrt WT');
hdb.read([tabExp.Efficacy tabExp.LogEC50]',tabExp.Label,'Experimental data');
hdb.read([tabExpDeltas.dEfficacy tabExpDeltas.dLogEC50]',tabExpDeltas.Label, 'Experiment d(Mut-WT)','VarNames',["dEfficacy" "dlog(EC50)"]);

hdb.align;
hdb.alignment

 % Run metaContactMap first
hDBNdx = 2;
contactCutMeta;
meta_contacts;
liProtRes;
alignHere = hdb.alignment.(string(hdb.entryName(hDBNdx)));

dataHere = hdb.hyperEntries{hDBNdx}.rawData(liProtRes,alignHere)';
dataHere(contact_map<contactCutMeta) = nan; % Remove data with low contact

X = round(dataHere,2,"significant");
figure
h = heatmap(hdb.hyperEntries{1}.VarNames,hdb.hyperEntries{1}.labels,X,'Colormap',turbo);

% xlabel('Protein residue', 'FontSize', fontsz)
h.Title =  hdb.hyperEntries{hDBNdx}.name;
h.XLabel = 'Protein residue';
h.YLabel = 'Ligand';
h.MissingDataColor=[0.8 0.8 0.8];
set(gca,'FontSize',12)

%% Combining data from different entries 2: Correlate experimental data with summed MI 
ligandName = "BRC"; % "DA" or "BRC"
calcMode = "noH8"; % "all" "noH8" "TM"
propertyName = "Efficacy"; % "Efficacy" or "Potency"
if strcmp(propertyName,"Efficacy")
    dataNdx = 1;
elseif strcmp(propertyName,"Potency")
    dataNdx = 3;
end

helices = settings.helices;
helicalResidues = [];
for i = 1:length(helices)
    helicalResidues = [helicalResidues helices(i,1):helices(i,2)];
end

dataTable = hdb.dbFetchData([ligandName "Gi"],"hyperEntryNdx",[1 2 4]);
labelNamesShort = extractBetween(dataTable.Label,ligandName + "-","-Gi");

x=zeros(length(dataTable.MI),1);
if strcmp(calcMode,'all')
    x = arrayfun(@(a) sum(a{1},'all'),dataTable.MI);
elseif strcmp(calcMode,'noH8')
    x = arrayfun(@(a) sum(a{1}(1:settings.helices(end)+1,1:settings.helices(end)+1),'all'),dataTable.MI);
elseif  strcmp(calcMode,'TM')
    x = arrayfun(@(a) sum(a{1}(helicalResidues,helicalResidues),'all'),dataTable.MI);
else
    sprintf('Invalid mode!')
end

%legacy
% dataTable = hdb.dbFetchData([ligandName "Gi"],"hyperEntryNdx",[1 2 3]);
% x = sum(dataTable.("MI per residue"),2); 
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
yneg = dataTable.("Experimental data")(:,dataNdx+1);
ypos = yneg;
xneg = zeros(length(y),1);
xpos = xneg;
errorbar(x,y,yneg,ypos,xneg,xpos,"LineStyle","none","Marker","none","Color",paperColorPalette.Hex(1))
hold on

% Fitting: 
y(isnan(y)) = 0;
X = [ones(length(x),1) x];
r = X\y;
yCalc2 = X*r;
plot(x,yCalc2,'-','Color',paperColorPalette.Hex(3),'LineWidth',1.5)
Rsq = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);

legend('Data', '' ,['Fit, R^2=' num2str(Rsq)],'Location','best');

xlabel(sprintf('Summed MI %s',calcMode))
ylabel(sprintf('Variant %s',propertyName))
title(sprintf('%s %s',ligandName, propertyName))
formatplot2
%% Combining data from hyperdatabase and external input: calculate hubscores/MI of different helices:

% Helices from D2-6VMS model
helices = [     2    30
    37    67
    73   107
   118   143
   156   189
   198   234
   240   264];

dataTable = hdb.dbFetchData("","hyperEntryNdx",[1 3 4]);

effector = "Gi";
ligandNames =["DA" "BRC"];
scoreHelixSum = zeros(sum(contains(dataTable.Label,effector)),size(helices,1));
scoreHelixLabels = cell(sum(contains(dataTable.Label,effector)),1);
count = 1;
for i = 1:length(ligandNames)
    ndxHere =  contains(dataTable.Label,effector) & contains(dataTable.Label,ligandNames(i)); % rows belonging to this category
    ndxHere = find(ndxHere);
    for j = 1:length(ndxHere)
        for k = 1:size(helices,1)
            tempData = dataTable.(3);
            scoreHelixSum(count,k) = sum(tempData(ndxHere(j),helices(k,1):helices(k,2)));
            scoreHelixLabels{count} = dataTable.Label(ndxHere(j));
        end
        count = count + 1;
    end
end

resultTable = table(string(scoreHelixLabels),scoreHelixSum);
% writetable(resultTable,"MIsum_helices_Gi.xlsx")
