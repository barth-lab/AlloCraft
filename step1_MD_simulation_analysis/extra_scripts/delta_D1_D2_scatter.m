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
deltaD1DA = zeros(length(mutList),4);
deltaD1BRC =  zeros(length(mutList),4);
deltaD2DA = zeros(length(mutList),4);
deltaD2BRC = zeros(length(mutList),4);

var2plot = [2 1 5 4]; %[2 1 5 4] to plot raw efficacies and [2 3 5 6] for normalized efficacies
for i = 1:length(mutList)
    dataTableD2 = hdb.dbFetchData([string(mutList(i)) "Gi"],"hyperEntryNdx",[1 2]);
    if strcmp(string(mutList(i)),"C6.47L")% Care for C6.47L!
        dataTableD2([1,2],:) = [];
    end
    dataTableD1 = D1hdb.dbFetchData([string(mutList(i)) "Gs"],"hyperEntryNdx",[1 2]);

    deltaD1DA(i,:) = dataTableD1.("Experiment d(Mut-WT)")(1,var2plot);
    deltaD1BRC(i,:) = dataTableD1.("Experiment d(Mut-WT)")(2,var2plot);
    deltaD2DA(i,:) = dataTableD2.("Experiment d(Mut-WT)")(1,var2plot);
    deltaD2BRC(i,:) = dataTableD2.("Experiment d(Mut-WT)")(2,var2plot);
end

deltaD1BRC(isnan(deltaD1BRC))=0;
deltaD2BRC(isnan(deltaD2BRC))=0;
table(mutList,deltaD1DA)

%
figure;
tiledlayout flow
nexttile
xRange = range([deltaD2DA(:,1); deltaD2BRC(:,1)]);
yRange = range([deltaD2DA(:,2); deltaD2BRC(:,2)]);
markersList = 'odod';

for i =1:4
    if i ==1
    deltadeltaHere = deltaD2DA;
    colorNdx = 2;
    elseif i==2
    deltadeltaHere = deltaD2BRC;
    colorNdx = 2;
    elseif i==3
    deltadeltaHere = deltaD1DA;
    colorNdx = 1;
    elseif i==4
    deltadeltaHere = deltaD1BRC;
    colorNdx = 1;
    end   
    s =scatter(deltadeltaHere(:,1),deltadeltaHere(:,2),150,'LineWidth',1.5,'MarkerFaceColor',paperColorPalette.Hex(2*colorNdx+4), ...
   'MarkerEdgeColor',paperColorPalette.Hex(2*colorNdx-1+4),'Marker',markersList(i));
    hold on
    row = dataTipTextRow('Label',mutList);
    s.DataTipTemplate.DataTipRows(end+1) = row;
    
    xtextPos = deltadeltaHere(:,1)+xRange/40;
    ytextPos = deltadeltaHere(:,2)+yRange/15;
    if i == 1
    ytextPos(end-1) = deltadeltaHere(end-1,2)-yRange/10; % Fix text for 3.41H
    xtextPos(end-1) = deltadeltaHere(end-1,1)-xRange/15; % Fix text for 3.41H
    elseif i == 2
    xtextPos(end-1) = deltadeltaHere(end-1,1)-xRange/15; % Fix text for 3.41H 
    ytextPos(end-1) = deltadeltaHere(end-1,2)+yRange/8; % Fix text for 3.41H
    end
    text(xtextPos,ytextPos, mutList,'FontSize',14);   
     % Add error bars:
    yneg = deltadeltaHere(:,4);
    ypos = yneg;
    xneg = deltadeltaHere(:,3);
    xpos = xneg;
    errorbar(deltadeltaHere(:,1),deltadeltaHere(:,2),yneg,ypos,xneg,xpos,"LineStyle","none","Color",paperColorPalette.Hex(2*colorNdx-1+4))
end
legend('D2 DA','','D2 BRC','','D1 DA','','D1 BRC','Location','best')
xlabel('\DeltaLog(EC50)')
ylabel('\DeltaEfficacy')
title('\Delta (Mut - WT)')
formatplot2
set(gcf,'Position',[302.6000   97.0000  780.8000  607.2000]);
% 
% figure
% xRange = range([deltaD1DA(:,1); deltaD1BRC(:,1)]);
% yRange = range([deltaD1DA(:,2); deltaD1BRC(:,2)]);
% for i =1:2
%     if i ==1
%     deltadeltaHere = deltaD1DA;
%     else
%     deltadeltaHere = deltaD1BRC;
%     end    
%     s =scatter(deltadeltaHere(:,1),deltadeltaHere(:,2),150,'LineWidth',1.5,'MarkerFaceColor',paperColorPalette.Hex(2*i), ...
%    'MarkerEdgeColor',paperColorPalette.Hex(2*i-1));
%     hold on
%     row = dataTipTextRow('Label',mutList);
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     
%     text(deltadeltaHere(:,1)+xRange/50,deltadeltaHere(:,2)+yRange/50, mutList,'FontSize',14);   
% 
%     % Add error bars:
%     yneg = deltadeltaHere(:,4);
%     ypos = yneg;
%     xneg = deltadeltaHere(:,3);
%     xpos = xneg;
%     errorbar(deltadeltaHere(:,1),deltadeltaHere(:,2),yneg,ypos,xneg,xpos,"LineStyle","none","Color",paperColorPalette.Hex(2*i-1))
% end
% legend('DA','','BRC','Location','best')
% xlabel('ddLog(EC50)')
% ylabel('ddEfficacy')
% title('\Delta\Delta (D1 - D2) (Mut - WT)')
% formatplot2
