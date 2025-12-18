% MI
%% MI per residue
metadir = 'D:\Telework_library\dopamine_phase_3';
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
pdbName = 'prot.pdb';

foldersToStudy = {'1-d2_dop_WT','70b-D2_DA_V444I_VC_rescale','71b-D2_DA_I181F_BW561_VC_rescale'};

thisSysLabel = {'DA-WT-Gi','DA-V4.44I-Gi','DA-I5.61F-Gi'};

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
tabHubscores = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Compiled_AlloDy_data\hubscores_DA_SNPS.xlsx','Sheet','D2');
hubscoreLabel = {'DA-WT-Gi','DA-V4.44I-Gi','DA-I5.61F-Gi'};


% Contact map:
%% MetaContactMap:
useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment
foldersToStudy = {'1-d2_dop_WT','70b-D2_DA_V444I_VC_rescale','71b-D2_DA_I181F_BW561_VC_rescale'};

thisSysLabel ={'DA-WT-Gi','DA-V4.44I-Gi','DA-I5.61F-Gi'};
name = 'temp';
metaContactMap; % Outputs contact_map and liProtRes

contactLabel = thisSysLabel;

%% Read the data into the hyperDB
hdb = HyperDataBase;
% hdb.read([tabExp.Efficacy tabExp.LogEC50]',tabExp.Label,'Experimental data');
% hdb.read([tabExpDeltas.dEfficacy tabExpDeltas.dLogEC50 tabExpDeltas.dEfficacy_normalized]',tabExpDeltas.Label, 'Experiment d(Mut-WT)','VarNames',["dEfficacy" "dlog(EC50)" "dEfficacy_norm"]);


hdb.read(contact_map', contactLabel,'Contact Map','VarNames',resname_num)
hdb.read(MIPerRes,MIPerResSysLabel,'MI per residue','VarNames',tabHubscores.(2));
hdb.read(MIresAll',MIPerResSysLabel,'MI','VarNames',tabHubscores.(2));
hdb.read(tabHubscores{:,3:end-2},hubscoreLabel,'Hubscores','VarNames',tabHubscores.(2));
% hdb.read(tabConserved{:,2:end},hubscoreLabel,'Conserved motifs HS','VarNames',tabConserved.Motif);

hdb.align;
hdb.alignment

%% Combining data from different entries: plot heatmap with hubscores/MI for ligand binding residues:


 % Run metaContactMap first
hDBNdx = 4;
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

%% Plot delta-contacts from respective WTs
ligNames = [ "DA"];

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