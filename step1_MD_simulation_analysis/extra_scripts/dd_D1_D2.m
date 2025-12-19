% mutantMI = hdb.hyperEntries{4}.fetchData(["I4.46N" "Gi"]);
% WTMI = hdb.hyperEntries{4}.fetchData(["WT" "Gi"]);

%% Compare D1 and D2 systems:
% Load D1 and D2 hyperDBs
D1HyperDBLoad; % Outputs D1hdb
hyperDBLoad; % Outputs hdb
% Load colors
paperColors; % Loads paperColorPalette
%%


%% Delta D2 - D1:
ndxD2 = strcmp(hdb.entryName,'Experiment d(Mut-WT)');
ndxD1 = strcmp(D1hdb.entryName,'Experiment d(Mut-WT)');


mutList  = unique(tabExpD1.Mutation,'stable');
mutList{1} = 'WT';
deltadeltaDA = zeros(length(mutList),4);
deltadeltaBRC =  zeros(length(mutList),4);
deltadeltaD1 = zeros(length(mutList),4);
deltadeltaD2 = zeros(length(mutList),4);

var2plot = [2 3 5 6]; %[2 1 5 4] to plot raw efficacies and [2 3 5 6] for normalized efficacies
for i = 1:length(mutList)
    dataTableD2 = hdb.dbFetchData([string(mutList(i)) "Gi"],"hyperEntryNdx",[1 2]);
    if strcmp(string(mutList(i)),"C6.47L")% Care for C6.47L!
        dataTableD2([1,2],:) = [];
    end
    dataTableD1 = D1hdb.dbFetchData([string(mutList(i)) "Gs"],"hyperEntryNdx",[1 2]);

    deltadeltaDA(i,:) = dataTableD1.("Experiment d(Mut-WT)")(1,var2plot) - dataTableD2.("Experiment d(Mut-WT)")(1,var2plot);
    deltadeltaBRC(i,:) = dataTableD1.("Experiment d(Mut-WT)")(2,var2plot) - dataTableD2.("Experiment d(Mut-WT)")(2,var2plot);
    deltadeltaD1(i,:) = dataTableD1.("Experiment d(Mut-WT)")(2,var2plot) - dataTableD1.("Experiment d(Mut-WT)")(1,var2plot);
    deltadeltaD2(i,:) = dataTableD2.("Experiment d(Mut-WT)")(2,var2plot) - dataTableD2.("Experiment d(Mut-WT)")(1,var2plot);

end

deltadeltaBRC(isnan(deltadeltaBRC))=0;
deltadeltaD2(isnan(deltadeltaD2))=0;
table(mutList,deltadeltaDA)

%
figure;
tiledlayout flow
nexttile
xRange = range([deltadeltaD1(:,1); deltadeltaD2(:,1)]);
yRange = range([deltadeltaD1(:,2); deltadeltaD2(:,2)]);

for i =1:2
    if i ==1
    deltadeltaHere = deltadeltaD1;
    else
    deltadeltaHere = deltadeltaD2;
    end   
    s =scatter(deltadeltaHere(:,1),deltadeltaHere(:,2),150,'LineWidth',1.5,'MarkerFaceColor',paperColorPalette.Hex(2*i+4), ...
   'MarkerEdgeColor',paperColorPalette.Hex(2*i-1+4));
    hold on
    row = dataTipTextRow('Label',mutList);
    s.DataTipTemplate.DataTipRows(end+1) = row;
    
    ytextPos = deltadeltaHere(:,2)+yRange/25;
    ytextPos(end-1) = deltadeltaHere(end-1,2)-yRange/25; % Fix text for 3.41H
    text(deltadeltaHere(:,1)+xRange/40,ytextPos, mutList,'FontSize',14);   
     % Add error bars:
    yneg = deltadeltaHere(:,4);
    ypos = yneg;
    xneg = deltadeltaHere(:,3);
    xpos = xneg;
    errorbar(deltadeltaHere(:,1),deltadeltaHere(:,2),yneg,ypos,xneg,xpos,"LineStyle","none","Color",paperColorPalette.Hex(2*i-1+4))
end
legend('D1','','D2','Location','best')
xlabel('ddLog(EC50)')
ylabel('ddEfficacy')
title('\Delta\Delta (BRC - DA) (Mut - WT)')
formatplot2

figure
xRange = range([deltadeltaDA(:,1); deltadeltaBRC(:,1)]);
yRange = range([deltadeltaDA(:,2); deltadeltaBRC(:,2)]);
for i =1:2
    if i ==1
    deltadeltaHere = deltadeltaDA;
    else
    deltadeltaHere = deltadeltaBRC;
    end    
    s =scatter(deltadeltaHere(:,1),deltadeltaHere(:,2),150,'LineWidth',1.5,'MarkerFaceColor',paperColorPalette.Hex(2*i), ...
   'MarkerEdgeColor',paperColorPalette.Hex(2*i-1));
    hold on
    row = dataTipTextRow('Label',mutList);
    s.DataTipTemplate.DataTipRows(end+1) = row;
    
    text(deltadeltaHere(:,1)+xRange/50,deltadeltaHere(:,2)+yRange/50, mutList,'FontSize',14);   

    % Add error bars:
    yneg = deltadeltaHere(:,4);
    ypos = yneg;
    xneg = deltadeltaHere(:,3);
    xpos = xneg;
    errorbar(deltadeltaHere(:,1),deltadeltaHere(:,2),yneg,ypos,xneg,xpos,"LineStyle","none","Color",paperColorPalette.Hex(2*i-1))
end
legend('DA','','BRC','Location','best')
xlabel('ddLog(EC50)')
ylabel('ddEfficacy')
title('\Delta\Delta (D1 - D2) (Mut - WT)')
formatplot2
