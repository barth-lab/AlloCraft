% Figure 1

%% Panel a: DB scan clustering
% Data extracted from Daniel's plot
dbscanData = readtable('C:\Users\mahdi\Desktop\PhD EPFL\Manuscripts\D2 paper 2022\Data and figures\Daniel\dbscan_data_Extracted_from_prism.xlsx');

clusterNdx = zeros(size(dbscanData,1),1);
pc2Data =  zeros(size(dbscanData,1),1);
for i = 1:8 % Number of clusters
    datHere = dbscanData.("Group"+num2str(i));
    clusterNdx(~isnan(datHere)) = i;
    pc2Data(~isnan(datHere)) = datHere(~isnan(datHere));
end

figure;
s = scatter(dbscanData.PCA1,pc2Data,[],clusterNdx,'filled' );
row = [dataTipTextRow('Ligand',dbscanData.Var1) dataTipTextRow('Cluster',clusterNdx)];

s.DataTipTemplate.DataTipRows(end+1:end+2) = row;

%% Plot in another way:

colorOrder = [3 1 9 5 7];
colorsExtra = ["#4DBEEE" "#F58216" "#A2142F" ];
figure;
for i = 1:8 % Number of clusters
    datHere = dbscanData.("Group"+num2str(i));
    ndxHere = ~isnan(datHere);

%     clusterNdx(~isnan(datHere)) = i;
%     pc2Data(~isnan(datHere)) = datHere(~isnan(datHere));
    if i <= length(colorOrder)
        colorHere = paperColorPalette.Hex(colorOrder(i));
    else
        colorHere = colorsExtra(i-5);
    end

    s = scatter(dbscanData.PCA1(ndxHere),datHere(~isnan(datHere)),200,'filled', ...
        'MarkerFaceColor',colorHere,'MarkerEdgeColor','k');
    hold on
    row = [dataTipTextRow('Ligand',dbscanData.Var1(ndxHere)) dataTipTextRow('Cluster',clusterNdx(ndxHere))];
    
    s.DataTipTemplate.DataTipRows(end+1:end+2) = row;
end
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('PC 1','FontSize',40)
ylabel('PC 2','FontSize',40)
% title('Dopamine receptor agonist clustering','FontSize',20)
% set(gca,'FontSize',20)


%% OR do the DBscan yourself:
ligandData = readtable(fullfile('C:\Users\mahdi\Desktop\PhD EPFL\Manuscripts\D2 paper 2022\Data and figures\Daniel',"d2_project2_ligand_properties_joelib_new.csv"));
epsDB = 26; % 4.08 used by daniel
minpts = 1;
datHere = ligandData{:,5:end-1};

idx = dbscan(datHere,epsDB,minpts);

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(datHere);

figure
 s = scatter(SCORE(:,1),SCORE(:,2),200,idx,'filled', ...
       'MarkerEdgeColor','k');
    hold on
    row = [dataTipTextRow('Ligand',ligandData.name) dataTipTextRow('Cluster',idx)];
    
    s.DataTipTemplate.DataTipRows(end+1:end+2) = row;
xlabel('PC 1','FontSize',40)
ylabel('PC 2','FontSize',40)
title([num2str(max(idx)) ', \epsilon=' num2str(epsDB) ],'FontSize',40)

%% Plot the DBscan yourself:

colorOrder = [3 1 9 5 7];
colorsExtra = ["#4DBEEE" "#F58216" "#A2142F" ];
figure;
xDat = SCORE(:,1);
yDat = SCORE(:,2);
for i = 1:8 % Number of clusters

    ndxHere = idx==i;

%     clusterNdx(~isnan(datHere)) = i;
%     pc2Data(~isnan(datHere)) = datHere(~isnan(datHere));
    if i <= length(colorOrder)
        colorHere = paperColorPalette.Hex(colorOrder(i));
    else
        colorHere = colorsExtra(i-5);
    end

    s = scatter(xDat(ndxHere),yDat(ndxHere),200,'filled', ...
        'MarkerFaceColor',colorHere,'MarkerEdgeColor','k');
    hold on
    row = [dataTipTextRow('Ligand',ligandData.name(ndxHere)) dataTipTextRow('Cluster',i*ones(sum(ndxHere),1))];
    
    s.DataTipTemplate.DataTipRows(end+1:end+2) = row;
end
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('PC 1','FontSize',40)
ylabel('PC 2','FontSize',40)

%% Next
hyperDBLoad;
paperColors; % Loads paperColorPalette
figDir = "D:\Telework_library\dopamine_phase_3\a-analysis\";
% Difference in hubscores between BRC-WT-Gi and DA-WT-Gi
hubDiff = hdb.hyperEntries{4}.fetchData(["BRC" "WT" "Gi"]) - hdb.hyperEntries{4}.fetchData(["DA" "WT" "Gi"]);


%% Load mainEntry and mainChain from DA WT Gi
receptorAtoms = mainChain.getAtoms();
bfactor = zeros(1, mainEntry.atomCount);

bfactor(1, receptorAtoms) = hubDiff;
writepdb(fullfile(figDir, "prot_hubDiffWTGi.pdb"), mainEntry.pdb, [], 'default', bfactor);

%% Try to plot the structure in Matlab:
myColorMap = colormapCustom('div','inputColor1', paperColorPalette.Hex(1) ,'inputColor2', paperColorPalette.Hex(3),'isHex',true);
CAcoord = mainEntry.pdb.xyz(receptorAtoms,:);
Nres = length(receptorAtoms);
figure
p3 = plot3(CAcoord(:,2),CAcoord(:,1),CAcoord(:,3),...
'color',[0.8 0.5 0.],'Linewidth',3);
hold on

% Scatter the CAs colored by sequence
s3 = scatter3(CAcoord(:,2),CAcoord(:,1),CAcoord(:,3),...
    100,hubDiff,'o','filled','Linewidth',2);

resText = mainEntry.chains{1}.formatResidues([1:Nres],'BWonly',true);
row = [dataTipTextRow('Residue',resText) dataTipTextRow('\Delta hubScore',hubDiff)];

p3.DataTipTemplate.DataTipRows(end+1:end+2) = row;
s3.DataTipTemplate.DataTipRows(end+1:end+2) = row;

colormap(myColorMap)
colorbar
axis xy
set(gca,'ydir','reverse')
axis('equal')
axis off

%% Or use the nice new function that I made:
helices = settings.helices;
helicalResidues = [];
for i = 1:length(helices)
    helicalResidues = [helicalResidues helices(i,1):helices(i,2)];
end
plotGPCRsnake(mainChain,'colorData',hubDiff(helicalResidues),'colorMap',myColorMap,'colorDataName','\Delta hubscore','resText',false)