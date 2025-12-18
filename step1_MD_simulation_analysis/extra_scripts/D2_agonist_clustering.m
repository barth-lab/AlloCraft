% Ligand clustering:

% load ligand data:
ligData = readtable('C:\Users\mahdi\Desktop\PhD EPFL\Manuscripts\D2 paper 2022\Data and figures\Daniel\d2_project2_ligand_properties_joelib_new.csv');

% Do PCA
firstCol = 4;
lastCol = size(ligData,2) - 1;
data = ligData{:,firstCol:lastCol};
[coeff,p,latent,tsquared,explained,mu] = pca(data);

% Cluster using DBscan
epsilon = 4.08;
minpts = 1;
labels = dbscan(normalize(data),epsilon,minpts);
% labels = dbscan(p(:, 1:2),epsilon,minpts);

% Plot
figure
tiledlayout('flow')
nexttile
pcaRange = range(p,'all');
s = scatter(p(:, 1), p(:, 2), 15,labels, 'filled');
% text(p(:, 1)+pcaRange/50,p(:, 2)+pcaRange/50,settings.thisSysLabel,'FontSize',14)
% num2str(resAll)
% text(p(:, 1)+pcaRange/100,p(:, 2)+pcaRange/100,num2str(resAll),'FontSize',10)
row = dataTipTextRow('Ligand Name',ligData.name);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('Cluster',labels);
s.DataTipTemplate.DataTipRows(end+1) = row;

xlabel(['PC 1, ' num2str(explained(1),'%.2f') ' explained'], 'fontsize', 20);
ylabel(['PC 2, ' num2str(explained(2),'%.2f') ' explained'], 'fontsize', 20);
  
nexttile
biplot(coeff(:,1:2),'scores',p(:,1:2),'varlabels',ligData.Properties.VariableNames(firstCol:lastCol));
set(gca,'FontSize',20)

% sgtitle(['Data = ' plotName],'FontSize',20)

