%% D1 stuff:
%% Experimental data:
tabExpD1 = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data\experimental_data_labeled.xlsx','Sheet','D1');

tabExpDeltasD1 = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data\delta_experimental_data.xlsx','Sheet','D1');
paperColors;
%% Hubscores
tabHubscores = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Compiled_AlloDy_data\hubscores_compiled_latest.xlsx','Sheet','D1');

hubscoreLabel = {'DA-WT-Gs','DA-WT-PAM-Gs','DA-I4.46N-Gs','DA-F6.44M-Gs', ...
    'BRC-WT-Gs','BRC-F6.44M-Gs','BRC-F6.44I-Gs','BRC-I4.46N-Gs'};

%% MIperres

metadir = 'D:\Telework_library\gpcrdb_extra_simulations\d1_dpa';
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 

pdbName = 'prot.pdb';
foldersToStudy = {'a-d1r_dpa_gp','c-d1r_dpa_pam_gp','d-dpa_I125N_gp','e-dpa-F228M_gp', ...
   'g-d1r_brc_gp/b-d1_BRC_from_MD','h-d1r_brc_F228M_gp','i-d1r_brc_F228I_gp','j-d1r_brc_I125N_gp'  };
thisSysLabel = {'DA-WT-Gs','DA-WT-PAM-Gs','DA-I4.46N-Gs','DA-F6.44M-Gs', ...
    'BRC-WT-Gs','BRC-F6.44M-Gs','BRC-F6.44I-Gs','BRC-I4.46N-Gs'};


md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc';
MIPerRes = [];

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
    end
end
MIPerResSysLabel = thisSysLabel;

%% MetaContactMap:
useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment

foldersToStudy = {'a-d1r_dpa_gp','c-d1r_dpa_pam_gp','d-dpa_I125N_gp','e-dpa-F228M_gp', ...
   'g-d1r_brc_gp/b-d1_BRC_from_MD','h-d1r_brc_F228M_gp','i-d1r_brc_F228I_gp','j-d1r_brc_I125N_gp' };
numRunsSys = [5 5 5 5 5 5 5 5]; % number of runs for each system

thisSysLabel = {'DA-WT-Gs','DA-WT-PAM-Gs','DA-I4.46N-Gs','DA-F6.44M-Gs', ...
    'BRC-WT-Gs','BRC-F6.44M-Gs','BRC-F6.44I-Gs','BRC-I4.46N-Gs'};
metaContactMap; % Outputs contact_map and liProtRes

contactLabel = thisSysLabel;
%% Load the D1 hDB:

D1hdb = HyperDataBase;
D1hdb.read([tabExpD1.Efficacy tabExpD1.LogEC50 tabExpD1.Efficacy_SEM tabExpD1.LogEC50_SEM]',tabExpD1.Label,'Experimental data', ...
    'VarNames',["Efficacy" "log(EC50)" "Efficacy_SEM" "log(EC50)_SEM"]);
dexpData = [tabExpDeltasD1.dEfficacy tabExpDeltasD1.dLogEC50 tabExpDeltasD1.dEfficacy_normalized ...
    tabExpDeltasD1.dEfficacy_SEM tabExpDeltasD1.dLogEC50_SEM tabExpDeltasD1.dEfficacy_normalized_SEM]';
D1hdb.read(dexpData,tabExpDeltasD1.Label, 'Experiment d(Mut-WT)', ...
    'VarNames',["dEfficacy" "dlog(EC50)" "dEfficacy_norm" "dEfficacy_SEM" "dlog(EC50)_SEM" "dEfficacy_norm_SEM"]);
D1hdb.read(MIPerRes,MIPerResSysLabel,'MI per residue','VarNames',tabHubscores.Var2);
D1hdb.read(tabHubscores{:,3:end},hubscoreLabel,'Hubscores','VarNames',tabHubscores.Var2);
D1hdb.read(contact_map', contactLabel,'Contact Map','VarNames',resname_num)

D1hdb.align('refEntry',1);
D1hdb.alignment

%% Plot heatmap with hubscores/MI for ligand binding residues:
 % Run metaContactMap first
hDBNdx = 3;
contactCutMeta;
meta_contacts;
liProtRes;

D1hdbMIref = HyperDataBase;
D1hdbMIref.read([tabExpD1.Efficacy tabExpD1.LogEC50]',tabExpD1.Label,'Experimental data');
D1hdbMIref.read([tabExpDeltasD1.dEfficacy tabExpDeltasD1.dLogEC50 tabExpDeltasD1.dEfficacy_normalized]',tabExpDeltasD1.Label, 'Experiment d(Mut-WT)','VarNames',["dEfficacy" "dlog(EC50)" "dEfficacy_norm"]);
D1hdbMIref.read(MIPerRes,MIPerResSysLabel,'MI per residue','VarNames',tabHubscores.Var2);
D1hdbMIref.read(tabHubscores{:,3:end},hubscoreLabel,'Hubscores','VarNames',tabHubscores.Var2);
D1hdbMIref.read(contact_map', contactLabel,'Contact Map','VarNames',resname_num)

D1hdbMIref.align('refEntry',3);
D1hdbMIref.alignment

dataHere = D1hdbMIref.hyperEntries{hDBNdx}.rawData(liProtRes,:)';
dataHere(contact_map<contactCutMeta) = nan; % Remove data with low contact

X = round(dataHere,3,"significant");
figure
h = heatmap(D1hdbMIref.hyperEntries{hDBNdx}.VarNames(liProtRes),D1hdbMIref.hyperEntries{hDBNdx}.labels,X,'Colormap',turbo);

% xlabel('Protein residue', 'FontSize', fontsz)
h.Title =  D1hdbMIref.hyperEntries{hDBNdx}.name;
h.XLabel = 'Protein residue';
h.YLabel = 'Ligand';
h.MissingDataColor=[0.8 0.8 0.8];

set(gca,'FontSize',12)

%% Plot delta-contacts from respective WTs
ligNames = [ "DA" "BRC"];
entryNdx = 5; % For contacts

ddataContact = [];
labelsAll = [];
for i =1:length(ligNames)
    dataContact = D1hdbMIref.hyperEntries{entryNdx}.fetchData(ligNames(i));
    ddataContactHere = [dataContact(:,1) (dataContact(:,2:end) - dataContact(:,1))];
    labelsHere = D1hdbMIref.hyperEntries{entryNdx}.labels(D1hdbMIref.hyperEntries{entryNdx}.fetchLabels(ligNames(i)));
    ddataContact = [ddataContact ddataContactHere];
    labelsAll = [labelsAll; labelsHere];
end

ddataContact(abs(ddataContact)<0.01) = 0; % zero out very small differences 
X = round(ddataContact',2,"significant");
figure
h = heatmap(D1hdbMIref.hyperEntries{entryNdx}.VarNames,labelsAll,X,'Colormap',myColorLavenderGrey);

caxis([-max(abs(caxis)) max(abs(caxis))]); % Make colorbar symmetrical
% xlabel('Protein residue', 'FontSize', fontsz)
h.Title =  '\Delta Contact from WT';
h.XLabel = 'Protein residue';
h.YLabel = 'Ligand';
h.MissingDataColor=[0.8 0.8 0.8];
set(gca,'FontSize',12)

%% Regress property with experimental data: (This should be a class function)

ligNames = ["DA" "BRC"];
effectorNames = ["Gs"];
markers ={'o','d','s','^'};
% propHere =  "MI per residue";
propHere =  "Hubscores";

figure
tiledlayout('flow')
for ligHere = 1:length(ligNames)
%         ligNdx = contains(thisSysLabel(ndxCommon(:,2)),ligNames(ligHere));
        for effectorHere = 1:length(effectorNames)
            % Plot the stuff
            nexttile
            tabHere = D1hdb.dbFetchData([ligNames(ligHere) effectorNames(effectorHere)]);
%             effectorNdx = contains(thisSysLabel(ndxCommon(:,2)),effectorNames(effectorHere));

            Xhere = sum(tabHere.(propHere),2);
            Yhere = tabHere.("Experimental data")(:,1) ;

            s = scatter(Xhere, Yhere,50, 'MarkerFaceColor', paperColorPalette.Hex{2*ligHere},...
             'MarkerEdgeColor',paperColorPalette.Hex{2*ligHere-1} ,'Marker',markers{ligHere},'LineWidth',1.5, ...
               'DisplayName', ligNames(ligHere)+"-"+effectorNames(effectorHere));
            hold on
            row = dataTipTextRow('Label',tabHere.Label);
            s.DataTipTemplate.DataTipRows(end+1) = row;
            text(Xhere(:,1)+xRange/50,Yhere(:,1)+yRange/50, tabHere.Label);   
            
            % Regression:
            [b,bint,r,rint,stats] = regress(Yhere,[ones(size(Xhere)) Xhere]);
            yFit = b(1) + b(2)*linspace(min(Xhere),max(Xhere),100);
            
            plot(linspace(min(Xhere),max(Xhere),100),yFit,'--r')
            legend('Data',['Linear fit, R^2=' num2str(stats(1))],'Location','best')
%             legend boxoff
            % Labels and formatting
            xlabel('\Sigma(MI)')
            ylabel('Efficacy (normalized to WT)')
%             ylabel('logEC50')
            title( ligNames(ligHere)+"-"+effectorNames(effectorHere))
            set(gca,'FontSize',16)
        end
end
