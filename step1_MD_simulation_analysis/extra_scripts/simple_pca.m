% Input data?

% Hubscore tables?
tabHere = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Compiled_AlloDy_data\hubscores_compiled_latest.xlsx','Sheet','D2');
% tabHere = readtable('D:\Telework_library\gpcrdb_extra_simulations\d1_dpa\md2pathMeta\hub_scores_June_7_2023_forPCA.xlsx');


% MI per residue?
metadir = 'D:\Telework_library\dopamine_phase_3';
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
pdbName = 'prot.pdb';

% All
foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT', '2-d2_dop_T174M-C220L', '4-d2_bromo_T174M-C220L', ...
     '10-d2_dop_L214M','11-d2_brc_L214M','12-d2_dop_L92G', '13-d2_brc_L92G', ...
      '14-d2_dop_F217M', '15-d2_brc_F217M',  '16-d2_dop_I125N', '17-d2_brc_I125N', '18b-d2_brc_F217I', ...
      '19a-d2_dop_L92H','19b-d2_brc_L92H', '34-d2_dpa_barr2_WT_for_Real', '35-d2_brc_barr2_WT_forREAL',...
       '31-d2_dpa_barr2', '32-d2_brc_barr2', '36b-d2_dpa_barr2_F217M', '36a-d2_brc_barr2_F217M', '37b-d2_dpa_barr2_I125N', ...
      '37a-d2_brc_barr2_I125N', '38b-d2_dpa_barr2_L214M','38a-d2_brc_barr2_L214M','39b-d2_dpa_barr2_L92H', ...
      '51b-d2_dpa_barr2_T174M-C220L','21-d2_ris_6cm4'};

thisSysLabel = {'DA-WT-Gi','BRC-WT-Gi','DA-T5.54M-C6.47L-Gi','BRC-T5.54M-C6.47L-Gi', ...
    'DA-L6.41M-Gi','BRC-L6.41M-Gi','DA-L3.41G-Gi','BRC-L3.41G-Gi', ...
    'DA-F6.44M-Gi','BRC-F6.44M-Gi','DA-I4.46N-Gi','BRC-4.46N-Gi','BRC-F6.44I-Gi', ...
    'DA-L3.41H-Gi','BRC-L3.41H-Gi','DA-WT-Barr2', 'BRC-WT-Barr2','DA-L3.41G-Barr2', ...
    'BRC-L3.41G-Barr2','DA-F6.44M-Barr2','BRC-F6.44M-Barr2','DA-I4.46N-Barr2', ...
    'BRC-I4.46N-Barr2','DA-L6.41M-Barr2','BRC-L6.41M-Barr2','DA-L3.41H-Barr2', ...
    'DA-T5.54M-C6.47L-Barr2','RIS-inactive'};

% Gi only (with V4.51F)
foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT', '2-d2_dop_T174M-C220L', '4-d2_bromo_T174M-C220L', ...
     '8-d2_dop_V130F','9-d2_bromo_V130F','10-d2_dop_L214M','11-d2_brc_L214M','12-d2_dop_L92G', '13-d2_brc_L92G', ...
      '14-d2_dop_F217M', '15-d2_brc_F217M',  '16-d2_dop_I125N', '17-d2_brc_I125N', '18b-d2_brc_F217I', ...
      '19a-d2_dop_L92H','19b-d2_brc_L92H','21-d2_ris_6cm4'};

thisSysLabel = {'DA-WT-Gi','BRC-WT-Gi','DA-T5.54M-C6.47L-Gi','BRC-T5.54M-C6.47L-Gi', ...
    'DA-V4.51F-Gi','BRC-V4.51F-Gi','DA-L6.41M-Gi','BRC-L6.41M-Gi','DA-L3.41G-Gi','BRC-L3.41G-Gi', ...
    'DA-F6.44M-Gi','BRC-F6.44M-Gi','DA-I4.46N-Gi','BRC-4.46N-Gi','BRC-F6.44I-Gi', ...
    'DA-L3.41H-Gi','BRC-L3.41H-Gi','RIS-inactive'};
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

% Maybe KLDiv will show us something?
% Run metaKLDiv first!
%
% Notes:
% KLDiv wrt inactive state only shows two clusters separating Gi-bound and
% Barr2-bound simulations, which is not very informative I would say.
% Combining KLDiv with MIperres gives very similar result to KLDiv only
%
% Try KLDiv wrt WT same ligand same IC binding partner
%% Do the actual PCA + plots

% Normal selections
% data = sqrt(tabHere{:,3:end});
% systemLabels = tabHere.Properties.VariableNames(3:end);
% data = MIPerRes;
% systemLabels = thisSysLabel;
% data = [kl1Res{:}];
% systemLabels = settings.thisSysLabel;
% dataLabels = tabHere.Var2;

%% Use HDB for selections
dbNdx = 7;
hdbName = hdb;
% data =  sqrt(hdbName.hyperEntries{dbNdx}.rawData);
% data =  (hdbName.hyperEntries{dbNdx}.rawData);
% systemLabels =  hdbName.hyperEntries{dbNdx}.labels;
% dataLabels = hdbName.hyperEntries{dbNdx}.VarNames ;
plotName = hdbName.hyperEntries{dbNdx}.name;

tabhere = hdbName.dbFetchData("Gi","hyperEntryNdx",[1 dbNdx]);
data = tabhere.(plotName)';
systemLabels = tabhere.Label;

%% Try deltas between BRC and DA: (run HDB selection first)
mutationNames = [ "T5.54M-C6.47L" "L6.41M" "L3.41G" "L3.41H" "F6.44M" "I4.46N"];
deltaData = zeros(length(mutationNames),size(tabhere.(plotName),2));
tempData = tabhere.(plotName);
for i = 1:length(mutationNames)
    ndxDA = contains(tabhere.Label,["DA"]) & contains(tabhere.Label,mutationNames(i));
    ndxBRC = contains(tabhere.Label,["BRC"]) & contains(tabhere.Label,mutationNames(i));
    
    deltaData(i,:) = tempData(ndxBRC,:) - tempData(ndxDA,:);
end
data = deltaData';
systemLabels = mutationNames;

%% Combine deltas?
dbNdx = [  3 7];
hdbName = hdb;
mutationNames = [ "T5.54M-C6.47L" "L6.41M" "L3.41G" "L3.41H" "F6.44M" "I4.46N"];

data =[];
dataLabels = [];
plotName = cell(length(dbNdx),1);
tabhere = hdbName.dbFetchData("Gi","hyperEntryNdx",[1 dbNdx]);

for i = 1:length(dbNdx)
    plotName{i} = hdbName.hyperEntries{dbNdx(i)}.name;
    deltaData = zeros(length(mutationNames),size(tabhere.(plotName{i}),2));
    tempData = tabhere.(plotName{i});

    for j = 1:length(mutationNames)
        ndxDA = contains(tabhere.Label,["DA"]) & contains(tabhere.Label,mutationNames(j));
        ndxBRC = contains(tabhere.Label,["BRC"]) & contains(tabhere.Label,mutationNames(j));
        
        deltaData(j,:) = tempData(ndxBRC,:) - tempData(ndxDA,:);
    end
     data = [data deltaData/max(abs(deltaData),[],'All')];
     dataLabels = [dataLabels ; (string(tabHere.Var2 + "-" + plotName{i}))];
end
% data = data';
systemLabels = mutationNames;

%%%%%%%%%%%%%%%%%%%% For switching PCA data columns and rows
systemLabels = dataLabels;
dataLabels = mutationNames;
%% Try combining different data:
% Add more data as extra "systems"
dbNdx = [ 3 4 6 7];
hdbName = hdb;
% data =  sqrt(hdbName.hyperEntries{dbNdx}.rawData);
% data =  (hdbName.hyperEntries{dbNdx}.rawData);
% systemLabels =  hdbName.hyperEntries{dbNdx}.labels;
% dataLabels = hdbName.hyperEntries{dbNdx}.VarNames ;
data =[];
systemLabels = [];
tabhere = hdbName.dbFetchData("Gi","hyperEntryNdx",[1 dbNdx]);
plotName = cell(length(dbNdx),1);
for i = 1:length(dbNdx)
plotName{i} = hdbName.hyperEntries{dbNdx(i)}.name;
dataTemp = tabhere.(plotName{i})';
data = [data dataTemp/max(dataTemp,[],'All')]; % make data between 0 and 1
systemLabels = [systemLabels ;(tabhere.Label + "-" + plotName{i})];
end
plotName = string(plotName);
% systemLabels = tabhere.Label;

%% Add more data as longer vectors:
dbNdx = [  3 7];
hdbName = hdb;
% data =  sqrt(hdbName.hyperEntries{dbNdx}.rawData);
% data =  (hdbName.hyperEntries{dbNdx}.rawData);
% systemLabels =  hdbName.hyperEntries{dbNdx}.labels;
% dataLabels = hdbName.hyperEntries{dbNdx}.VarNames ;
data =[];
dataLabels = [];
tabhere = hdbName.dbFetchData("Gi","hyperEntryNdx",[1 dbNdx]);
plotName = cell(length(dbNdx),1);
for i = 1:length(dbNdx)
plotName{i} = hdbName.hyperEntries{dbNdx(i)}.name;
dataTemp = tabhere.(plotName{i})';
data = [data ;dataTemp/max(dataTemp,[],'All')]; % make data between 0 and 1
dataLabels = [dataLabels ; (string(tabHere.Var2 + "-" + plotName{i}))];
end

plotName = string(plotName);
systemLabels = tabhere.Label;

% Transpose data?
data = data';
systemLabels = dataLabels;
dataLabels = tabhere.Label;

% Helical residues only?
% data = data(helicalResidues,:)';
% systemLabels = dataLabels(helicalResidues);
% dataLabels = tabhere.Label;
%%
[coeff,p,latent,tsquared,explained,mu] = pca(data);

figure
tiledlayout('flow')
nexttile
pcaRange = range(p,'all');
s = scatter(p(:, 1), p(:, 2), 15,1:size(p,1), 'filled');
% text(p(:, 1)+pcaRange/50,p(:, 2)+pcaRange/50,settings.thisSysLabel,'FontSize',14)
% num2str(resAll)
% text(p(:, 1)+pcaRange/100,p(:, 2)+pcaRange/100,num2str(resAll),'FontSize',10)
row = dataTipTextRow('Residue',dataLabels);
s.DataTipTemplate.DataTipRows(end+1) = row;

xlabel(['PC 1, ' num2str(explained(1),'%.2f') ' explained'], 'fontsize', 20);
ylabel(['PC 2, ' num2str(explained(2),'%.2f') ' explained'], 'fontsize', 20);
  
text(p(:, 1)+pcaRange/50,p(:, 2)+pcaRange/50,dataLabels,'FontSize',11)

nexttile
biplot(coeff(:,1:2),'scores',p(:,1:2),'varlabels',systemLabels);
set(gca,'FontSize',20)

if length(plotName) == 1
 sgtitle(['Data = ' plotName],'FontSize',20)
else
    sgtitle(plotName,'FontSize',20)
end
% sgtitle('Data = MIPerRes','FontSize',20)
% sgtitle(['Data = KLDiv|' settings.refName],'FontSize',20)
