%% setting up part:

options.SavePath = 'C:\Users\mahdi\Desktop\PhD EPFL\Manuscripts\D2 paper 2022\New submission 2023\deliverables\hubscores_and_pipelines\hubscore_pairwise_comparison';
options.propNdx = 4;
NdxRef = 15;
NdxTest = 17;

obj = hdb.hyperEntries{options.propNdx};
%% Plotting part
figure('Position',1.0e3*[0.3994    0.0762    1.1056    0.7000]);
tiledlayout('flow')
nexttile
ligHere = 1;
 s1 =scatter(mainChain.resIds, obj.rawData(:,NdxTest),25, ...
               'LineWidth',1.5, ...
               'MarkerEdgeColor',paperColorPalette.Hex(1), ...
               'MarkerFaceColor',paperColorPalette.Hex(2));
hold on
 s2 =scatter(mainChain.resIds, obj.rawData(:,NdxRef),25, ...
               'LineWidth',1.5, ...
               'MarkerEdgeColor',paperColorPalette.Hex(3), ...
               'MarkerFaceColor',paperColorPalette.Hex(4));
row = dataTipTextRow('VarNames',obj.VarNames);
s1.DataTipTemplate.DataTipRows(end+1) = row;
s2.DataTipTemplate.DataTipRows(end+1) = row;

% row = dataTipTextRow('Label',tabExp.Label(ligNdx&effectorNdx));
% s.DataTipTemplate.DataTipRows(end+1) = row
drawTMhelices(obj.rawData(:,NdxTest),helices, mainChain.resIds)
xlabel('Residue')
ylabel(obj.name)
legend([mainEntry.name ', \Sigma MI(i,j) = ' num2str(round(sum(obj.rawData(:,NdxTest))))], ...
      [refEntry.name ', \Sigma MI(i,j) = ' num2str(round(sum(obj.rawData(:,NdxRef))))],'Location','bestoutside')

title(sprintf("%s, %s, %s",obj.name,obj.labels{NdxTest},obj.labels{NdxRef}))
set(gca,'FontSize',16)

% Difference of MI plot:
% figure; 
nexttile
deltaData = obj.rawData(:,NdxTest)-obj.rawData(:,NdxRef);
p = plot(deltaData,'-o');
hold on
drawTMhelices(deltaData,helices, mainChain.resIds)
% row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds));
p.DataTipTemplate.DataTipRows(end+1) = row;
grid on

legend_entries{1} = ['\Delta ' obj.name];
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
if exist(fullfile(settings.refdir,'md2path','BS_residues.txt'),'file')
    bsRes.ref = importdata(fullfile(settings.refdir,'md2path','BS_residues.txt'));

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
ylabel(sprintf("delta %s, (%s - %s)",obj.name,obj.labels{NdxTest},obj.labels{NdxRef}))
title([obj.name ', ' obj.labels{NdxTest} ' - ' obj.labels{NdxRef} ])
set(gca,'FontSize',16)

if ~isempty(options.SavePath)
figPath = fullfile(options.SavePath,sprintf("delta%s_%s_%s",obj.name,obj.labels{NdxTest},obj.labels{NdxRef}));
savefig(figPath + ".fig");
print2pdf(figPath+ ".pdf");
saveas(gcf,figPath + ".svg",'svg')
end