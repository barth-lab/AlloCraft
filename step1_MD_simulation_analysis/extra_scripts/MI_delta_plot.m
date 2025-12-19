% Load the MI
if ~exist('alloPathCalcName','var')
    alloPathCalcName = 'alloPathCalc_Culled_data';
end
md2pathName = 'md2pathdev';

temp = load(fullfile(settings.mydir,'md2pathdev',alloPathCalcName,'workspace.mat'),'-mat','MIres');
MIres = temp.MIres;
temp = load(fullfile(settings.refdir,'md2pathdev',alloPathCalcName,'workspace.mat'),'-mat','MIres');
MIresWT = temp.MIres;
helices = settings.helices;

paperColors; % Loads paperColorPalette
%% MI microscope
resHere = [125]; % in system 1
refNdx = 2;

alignmentArray = table2array(database.residues{1}(:,[1 refNdx])); % Assuming target system is always indexed 1
resHereRef = alignmentArray(alignmentArray(:,1)==resHere,refNdx);
assert(resHereRef~=0, "No alignment found for this residue!")

pos = [289.8000  234.6000  948.8000  478.4000];
figure('Position',pos); 
% f = figure;
p1 = plot(MIres(resHere,:));
hold on
p2 = plot(MIresWT(resHereRef,:));
row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds));
rowRef = dataTipTextRow('Residue',refChain.formatResidues(refChain.resIds));
p1.DataTipTemplate.DataTipRows(end+1) = row;
p2.DataTipTemplate.DataTipRows(end+1) = rowRef;
xlabel('Residue')
ylabel('MI')
title( settings.systemName + " " + mainEntry.name + " Gi, " + mainChain.formatResidues(resHere,'BWonly',1))

MIresHere = MIres(resHere,:);
legend_entries{1} = [mainEntry.name ', MI=' num2str(sum(MIres(resHere,:)))];
legend_entries{2} = [refEntry.name ', MI=' num2str(sum(MIresWT(resHereRef,:)))];
legend_count = 3;

% Add Annotations to the figure
if exist(fullfile(md2pathdir,'BS_residues.txt'),'file')
    bsRes.main = importdata(fullfile(md2pathdir,'BS_residues.txt'));
    scatter(bsRes.main,MIresHere(receptorResIdsNdx(bsRes.main)), 'LineWidth',1.5)
    legend_entries{legend_count} = 'Binding residues';
    legend_count = legend_count + 1;
end
if exist(fullfile(settings.refdir,md2pathName,'BS_residues.txt'),'file')
    bsRes.ref = importdata(fullfile(settings.refdir,md2pathName,'BS_residues.txt'));

    % Translate bsRef to main numbering before plotting:
    bsResRefRenum = zeros(size(bsRes.ref));
    for resi = 1:length(bsRes.ref)
        bsResRefRenum(resi) = alignmentArray(alignmentArray(:,2)==bsRes.ref(resi),1);
    end
    scatter(bsResRefRenum,MIresHere(bsResRefRenum), 25, 'filled')
    legend_entries{legend_count} = 'Ref binding residues';
    legend_count = legend_count + 1;
end

if ~isempty(mutations)
    s3 = scatter(mutRes, MIresHere(receptorResIdsNdx(mutRes)),80, 'LineWidth',1.5);
    text(mutRes,MIresHere(receptorResIdsNdx(mutRes)) +0.05*max(MIresHere(resHere)),mutations)
    row = dataTipTextRow('Mut',mutations);
    s3.DataTipTemplate.DataTipRows(end+1) = row;

    legend_entries{legend_count} = 'Mutation sites';
    legend_count = legend_count + 1;
end

if ~isempty(settings.highlightRes)
    scatter(settings.highlightRes, MIresHere(receptorResIdsNdx(settings.highlightRes)),60,'x', 'LineWidth',1.5);
    legend_entries{legend_count} = settings.highlightText;
    legend_count = legend_count + 1;
end

legend(legend_entries)
drawTMhelices(MIres(resHere,:),helices, mainChain.resIds)
set(gca,'FontSize',16)
% formatplot2
% saveas(gcf,fullfile("D:\Telework_library\dopamine_phase_3\a-analysis\mutations\I446N","MI D1 " + mainEntry.name + " Gs " + num2str(resHere) +".svg"),'svg')
% saveas(gcf,fullfile("D:\Telework_library\dopamine_phase_3\a-analysis\D1\PAM","MI D2 " + mainEntry.name + " Gi " + num2str(resHere) +".svg"),'svg')


%% Delta MI over the whole receptor

if ~exist('receptorResIds','var')
    receptorResIds = mainChain.resIds;
end

figure('Position',1.0e3*[0.3994    0.0762    1.1056    0.7000]);
tiledlayout('flow')
nexttile
ligHere = 1;
 s1 =scatter(mainChain.resIds, sum(MIres),25, ...
               'LineWidth',1.5, ...
               'MarkerEdgeColor',paperColorPalette.Hex(1), ...
               'MarkerFaceColor',paperColorPalette.Hex(2));
hold on
 s2 =scatter(mainChain.resIds, sum(MIresWT),25, ...
               'LineWidth',1.5, ...
               'MarkerEdgeColor',paperColorPalette.Hex(3), ...
               'MarkerFaceColor',paperColorPalette.Hex(4));
row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds));
s1.DataTipTemplate.DataTipRows(end+1) = row;
s2.DataTipTemplate.DataTipRows(end+1) = row;

% row = dataTipTextRow('Label',tabExp.Label(ligNdx&effectorNdx));
% s.DataTipTemplate.DataTipRows(end+1) = row
drawTMhelices(sum(MIres),helices, mainChain.resIds)
xlabel('Residue')
ylabel('MI')
legend([mainEntry.name ', \Sigma MI(i,j) = ' num2str(round(sum(sum(MIres))))], ...
    [refEntry.name ', \Sigma MI(i,j) = ' num2str(round(sum(sum(MIresWT))))],'Location','bestoutside')
title(['Summed MI per residue, ' mainEntry.name ', ' refEntry.name ])
set(gca,'FontSize',16)

% Difference of MI plot:
% figure; 
nexttile
sMIres = sum(MIres);
sMIresWT = sum(MIresWT);

deltaData = sMIres-sMIresWT;
deltaDataNorm = zeros(length(sMIresWT),1);
for i = 1:length(sMIresWT)
    deltaDataNorm(i) =  (sMIres(i) - sMIresWT(i))/(0.5*(sMIresWT(i) +sMIres(i)));
end

p = plot(deltaData,'-o');
hold on
drawTMhelices(deltaData,helices, mainChain.resIds)
row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds));
p.DataTipTemplate.DataTipRows(end+1) = row;
grid on

legend_entries{1} = '\Delta MI';
hold on
legend_count = 2;

receptorResIdsNdx = zeros(max(receptorResIds),1); % Helpful for indexing
receptorResIdsNdx(receptorResIds)=1:length(receptorResIds); % Hoping this would do the trick

% Add mutations (from reference structure)
mutPos = find(database.residues{1}.Name(:,1) ~= database.residues{1}.Name(:,2));
% Remove any difference in length where alignment is not present, those are
% probably not mutations
mutPos(isspace(database.residues{1}.Name(mutPos,2)) | isspace(database.residues{1}.Name(mutPos,1))) = [];
mutRes = table2array(database.residues{1}(mutPos,1)); % Mutated residues
mutations = cellstr([database.residues{1}.Name(mutPos,2) num2str(mutRes) database.residues{1}.Name(mutPos,1)  ]);

if settings.includeLigPathwayCalc % do same for peptide ligand
    mutPosLig = find(database.residues{3}.Name(:,1) ~= database.residues{3}.Name(:,2));
    mutations = [mutations ; cellstr( [database.residues{3}.Name(mutPosLig,2) num2str(table2array(database.residues{3}(mutPosLig,1))) database.residues{3}.Name(mutPosLig,1)])  ];
    mutRes = [mutRes ; table2array(database.residues{3}(mutPosLig,1))];
end

% Add Annotations to the figure
if exist(fullfile(md2pathdir,'BS_residues.txt'),'file')
    bsRes.main = importdata(fullfile(md2pathdir,'BS_residues.txt'));
    scatter(bsRes.main,deltaData(receptorResIdsNdx(bsRes.main)), 'LineWidth',1.5)
    legend_entries{legend_count} = 'Binding residues';
    legend_count = legend_count + 1;
end
if exist(fullfile(settings.refdir,md2pathName,'BS_residues.txt'),'file')
    bsRes.ref = importdata(fullfile(settings.refdir,md2pathName,'BS_residues.txt'));

    % Translate bsRef to main numbering before plotting:
    ndxMainRef = [database.residues{1}.output1 database.residues{1}.output2];
    bsResRefRenum = zeros(size(bsRes.ref));
    for resi = 1:length(bsRes.ref)
        bsResRefRenum(resi) = ndxMainRef(ndxMainRef(:,2)==bsRes.ref(resi),1);
    end
    scatter(bsResRefRenum,deltaData(bsResRefRenum), 25, 'filled')
    legend_entries{legend_count} = 'Ref binding residues';
    legend_count = legend_count + 1;
end

if ~isempty(mutations)
    s3 = scatter(mutRes, deltaData(receptorResIdsNdx(mutRes)),80, 'LineWidth',1.5);
    text(mutRes,deltaData(receptorResIdsNdx(mutRes)) +0.05*max(deltaData),mutations)
    row = dataTipTextRow('Mut',mutations);
    s3.DataTipTemplate.DataTipRows(end+1) = row;

    legend_entries{legend_count} = 'Mutation sites';
    legend_count = legend_count + 1;
end

if ~isempty(settings.highlightRes)
    scatter(settings.highlightRes, deltaData(receptorResIdsNdx(settings.highlightRes)),60,'x', 'LineWidth',1.5);
    legend_entries{legend_count} = settings.highlightText;
    legend_count = legend_count + 1;
end

legend(legend_entries,'Location','bestoutside')
xlabel('Residue')
ylabel('\Delta MI (Mut - WT)')
title(['Summed MI per residue, ' mainEntry.name ' - ' refEntry.name ])
set(gca,'FontSize',16)
% 
% figPath = fullfile('C:\Users\mahdi\Desktop\PhD EPFL\Manuscripts\D2 paper 2022\New submission 2023\deliverables\MI\D1',['deltaMI_' mainEntry.name '_' refEntry.name]);
figPath = fullfile(md2pathdir,['deltaMI_' mainEntry.name '_' refEntry.name]);
savefig(figPath + ".fig");
print2pdf(figPath+ ".pdf");
saveas(gcf,figPath + ".svg",'svg')





%% Save deltaMI as B-factor:

bfactor = zeros(1, mainEntry.atomCount);
bfactordNorm = zeros(1, mainEntry.atomCount);

if settings.includeLigPathwayCalc % Include ligand in b-factor array
    chainSelect = selectname(mainEntry.pdb.chainid, settings.chains(Chains.receptor)) | selectname(mainEntry.pdb.chainid, settings.chains(Chains.ligand));
    selection = selectname( mainEntry.pdb.name,'CA') & chainSelect;
else
    selection = selectname( mainEntry.pdb.name,'CA') & selectname(mainEntry.pdb.chainid, settings.chains(Chains.receptor));
end
% resFromReSort = selectid(mainEntry.pdb.resseq,unique(reSortCommon(:,1))); % Only residues contained in reSortCommon
% selection = selection & resFromReSort;

bfactor(selection) = deltaData;
writepdb(fullfile(md2pathdir, ['prot_dMI_' mainEntry.name '_' refEntry.name '.pdb']), mainEntry.pdb, [], 'default', bfactor);

bfactordNorm(selection) = deltaDataNorm;
writepdb(fullfile(md2pathdir, ['prot_dMINorm_' mainEntry.name '_' refEntry.name '.pdb']), mainEntry.pdb, [], 'default', bfactordNorm);

% 

%% Delta MI, heatmap:
figure('Position',[127.4 128.2 739.2 574.4]);
imagesc(MIres - MIresWT)

ylabel('Residue')
xlabel('Residue')
title(['\Delta MI per residue, ' mainEntry.name ' - ' refEntry.name ])
set(gca,'FontSize',16)
colorbar

drawTMhelices( mainChain.resIds,settings.helices, mainChain.resIds)
drawTMhelices( mainChain.resIds,settings.helices, mainChain.resIds,'Y')

% colormapCustom
colormap(myColorLavenderGrey)
% colormap(myColorOlivePurple)
% caxis([-max(abs(caxis)) max(abs(caxis))])
caxis([-1 1])

%% Try to make it more informative?


% t = tiledlayout(3,2);
% ax1 = nexttile([1 2]);
%  s = scatter(mainChain.resIds, deltaData,25,normalize(deltaData, 'range'),'filled');
% colormap(myColorLavenderGrey)
% axis([-26.7000  278.5000 -1.2*abs(min(deltaData)) 1.2*max(deltaData)])
% ylabel('\Delta \Sigma MI (Mut - WT)')
% set(gca,'FontSize',16)


figure('Position',[488.0000  137.8000  735.4000  624.2000]);

% t = tiledlayout(2,3);


t = tiledlayout(3,3);

ax1 = nexttile([1 2]);
% barplot showing MI per helix:
MIresSummed = sum(MIres);
MIresWTSummed= sum(MIresWT);
MIhelix = zeros(size(settings.helices,1),1);
MIhelixWT = zeros(size(settings.helices,1),1);
barCenters = zeros(size(settings.helices,1),1);
barCenterLabel = categorical(["TM1" "TM2" "TM3" "TM4" "TM5" "TM6" "TM7"]);

 sum(MIresSummed(settings.helices(1,1):settings.helices(1,2)))

 for i = 1: size(settings.helices,1)
    MIhelix(i) =  sum(MIresSummed(settings.helices(i,1):settings.helices(i,2)));
    MIhelixWT(i) =  sum(MIresWTSummed(settings.helices(i,1):settings.helices(i,2)));
    barCenters(i) = round((settings.helices(i,1) +settings.helices(i,2))/2);
 end


 b = bar(barCenters,[MIhelixWT MIhelix],'FaceColor','flat','FaceAlpha',0.8 );
b(1).CData = 1;
b(2).CData = 2;
colormap(myColorLavenderGrey)
set(gca,'FontSize',16)
ylabel('TM summed MI')
legend(refEntry.name,  mainEntry.name,'Location','best','fontsize',8);

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData));
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(round(b(2).YData));
text(xtips1,1.05*ytips1,labels1,'HorizontalAlignment','left',...
    'VerticalAlignment','middle')
text(xtips2,1.05*ytips2,labels2,'HorizontalAlignment','left',...
    'VerticalAlignment','middle')
H=findobj(gca,'Type','text');
set(H,'Rotation',90); % tilt
yticklabels(ax1,{})
xticklabels(ax1,{})
temp = axis;
axis([  -26.7000  length(MIres)+.5000 temp(3) 1.1*temp(4)])
box off

nexttile([2 2])

imagesc(MIres - MIresWT)

ylabel('Residue')
xlabel('Residue')
% title(['\Delta MI per residue, ' mainEntry.name ' - ' refEntry.name ])
set(gca,'FontSize',16)
% colorbar

drawTMhelices( mainChain.resIds,settings.helices, mainChain.resIds)
drawTMhelices( mainChain.resIds,settings.helices, mainChain.resIds,'Y')
axis([-26.7000  length(MIres)+.5000  -26.7000  length(MIres)+.5000]);
% colormapCustom
colormap(myColorLavenderGrey)
% colormap(myColorOlivePurple)
% caxis([-max(abs(caxis)) max(abs(caxis))])
caxis([-1 1])


ax2 = nexttile(6,[2 1]);
barh(deltaData,'FaceColor',[0.8 0.8 0.8])
box off

hold on
 s = scatter(deltaData, mainChain.resIds ,25,normalize(deltaData, 'range'),'filled');
 row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds));
 s.DataTipTemplate.DataTipRows(end+1) = row;
colormap(myColorLavenderGrey)
hold on
if ~isempty(mutations)
    s3 = scatter(deltaData(receptorResIdsNdx(mutRes)),mutRes,50, 'LineWidth',1.5);
%     text(deltaData(receptorResIdsNdx(mutRes)) +0.2*max(deltaData),mutRes,mutations)
    row = dataTipTextRow('Mut',mutations);
    s3.DataTipTemplate.DataTipRows(end+1) = row;
end

axis([ -1.2*abs(min(deltaData)) 1.2*max(deltaData) -26.7000  length(MIres)+.5000])
xlabel('\Delta \Sigma MI (Mut - WT)')
set(gca,'FontSize',16)
set(gca, 'YDir','reverse')
cb = colorbar;
cb.TickLabels = {};
% Move plots closer together
yticklabels(ax2,{})
t.TileSpacing = 'tight';

% for reference:
sgtitle(['\Delta MI, ' mainEntry.name ' - ' refEntry.name ])

figPath = fullfile(md2pathdir,['deltaMI_' mainEntry.name '_' refEntry.name '_map_bar']);
savefig(figPath + ".fig");
print2pdf(figPath+ ".pdf");
saveas(gcf,figPath + ".svg",'svg')

%% barplot showing MI per helix:

MIresSummed = sum(MIres);
MIresWTSummed= sum(MIresWT);
MIhelix = zeros(size(settings.helices,1),1);
MIhelixWT = zeros(size(settings.helices,1),1);
barCenters = zeros(size(settings.helices,1),1);
barCenterLabel = categorical(["TM1" "TM2" "TM3" "TM4" "TM5" "TM6" "TM7"]);

 sum(MIresSummed(settings.helices(1,1):settings.helices(1,2)))

 for i = 1: size(settings.helices,1)
    MIhelix(i) =  sum(MIresSummed(settings.helices(i,1):settings.helices(i,2)));
    MIhelixWT(i) =  sum(MIresWTSummed(settings.helices(i,1):settings.helices(i,2)));
    barCenters(i) = round((settings.helices(i,1) +settings.helices(i,2))/2);
 end

 figure
 b = bar(barCenterLabel,[MIhelixWT MIhelix],'FaceColor','flat','FaceAlpha',0.8 );
b(1).CData = 1;
b(2).CData = 2;
colormap(myColorLavenderGrey)
legend(refEntry.name,  mainEntry.name,'Location','best');
%% Show delta on snake plot:
helices = settings.helices;
helicalResidues = [];
for i = 1:length(helices)
    helicalResidues = [helicalResidues helices(i,1):helices(i,2)];
end
mainChain.helices = helices;

snakePlotGPCR(mainChain,'colorData',deltaData(helicalResidues),'colorDataName','\Delta MI')
set(gcf,'position',[485.0000  268.2000  832.8000  420]);
colorbar
colormap(myColorLavenderGrey)
title(['\Delta MI per residue, ' mainEntry.name ' - ' refEntry.name ])

figPath = fullfile(md2pathdir,['snake_deltaMI_' mainEntry.name '_' refEntry.name]);
savefig(figPath + ".fig");
print2pdf(figPath+ ".pdf");
saveas(gcf,figPath + ".svg",'svg')
%% Gp region sum MI?
gpiResidues = settings.highlightRes;
otherResidues = mainChain.resIds;

for i = 1:length(gpiResidues)
    otherResidues(otherResidues == gpiResidues(i)) = [];
end
% intersect(otherResidues,gpiResidues) % check whether it works

sum(MIres(gpiResidues,otherResidues),'all')

sum(MIresWT(gpiResidues,otherResidues),'all')
