tabExp = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data\experimental_data_labeled.xlsx');

% Load colors
paperColors; % Loads paperColorPalette
%%
figure;


s =scatter(tabExp.LogEC50,tabExp.Efficacy,25,'filled');
xlabel('Log(EC50)')
ylabel('Efficacy')
row = dataTipTextRow('Label',tabExp.Label);
s.DataTipTemplate.DataTipRows(end+1) = row;
formatplot2


%%
dxDataAll = [];
dyDataAll = [];
labelsdData = [];
ligNames = unique(string(tabExp.Ligand));
effectorNames = unique(string(tabExp.ICBindingPartner));
markers ={'o','d','s','^'};
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
xRange = range(tabExp.LogEC50);
yRange = range(tabExp.Efficacy);

pos = [112.2000  109.0000  772.8000  640.0000];
% figure('Position',pos)
% tiledlayout('flow')
for plotHere = 1:2
%     nexttile
figure('Position',pos)
    for ligHere = 1:length(ligNames)
        ligNdx = strcmp(string(tabExp.Ligand),ligNames(ligHere));
        for effectorHere = 1%:length(effectorNames)
            effectorNdx = strcmp(string(tabExp.ICBindingPartner),effectorNames(effectorHere));            
            xData = tabExp.LogEC50(ligNdx&effectorNdx);
            if strcmp(ligNames(ligHere),"BRC") 
                xData(isnan(xData)) = -6.33;
            end

            yData = tabExp.Efficacy(ligNdx&effectorNdx);

            WTNdx =  find(strcmp(string(tabExp.Mutation(ligNdx&effectorNdx)),'WT'));
            if plotHere == 1 % Plot raw values
            s =scatter(xData,yData,50, ...
               colors(ligHere,:),'Marker',markers{ligHere},'LineWidth',1.5, ...
               'DisplayName', ligNames(ligHere)+"-"+effectorNames(effectorHere),'MarkerEdgeColor',paperColorPalette.Hex(2*ligHere-1), ...
               'MarkerFaceColor',paperColorPalette.Hex(2*ligHere));
            % Add mutations as labels to each data point
            text(xData+xRange/100,yData+yRange/50, tabExp.Mutation(ligNdx&effectorNdx));   
            hold on
            else % Plot deltas
            dxData = xData - xData(WTNdx);
            dyData = yData - yData(WTNdx);
            
            dxDataAll = [dxDataAll ;dxData];
            dyDataAll = [dyDataAll ;dyData];
            labelsdData = [labelsdData ;tabExp.Label(ligNdx&effectorNdx)];

            s =scatter(dxData,dyData,100, ...
               'Marker',markers{ligHere},'LineWidth',1.5, ...
               'DisplayName', ligNames(ligHere)+"-"+effectorNames(effectorHere), ...
               'MarkerEdgeColor',paperColorPalette.Hex(2*ligHere-1), ...
               'MarkerFaceColor',paperColorPalette.Hex(2*ligHere));
            hold on
            for i = 1:length(dxData)
                drawArrow([0 dxData(i)],[0 dyData(i)],'color',paperColorPalette.Hex(2*ligHere-1),'LineWidth',1, ...
                   'ShowArrowHead','off','HandleVisibility','off')
            end
            % Add mutations as labels to each data point
            text(xData - xData(WTNdx)+xRange/50,yData - yData(WTNdx)+yRange/50, tabExp.Mutation(ligNdx&effectorNdx));   
            end
            row = dataTipTextRow('Label',tabExp.Label(ligNdx&effectorNdx));
            s.DataTipTemplate.DataTipRows(end+1) = row;
  
        end
    end
    legend('Location','best')
    legend boxoff
    if plotHere == 1
        xlabel('Log(EC50)','FontSize',20)
        ylabel('Efficacy','FontSize',20)
        title('Experimental efficacies and EC50s')
    else
        grid on
        xlabel('\Delta Log(EC50)','FontSize',20)
        ylabel('\Delta Efficacy','FontSize',20)
        title('\Delta(value - WT)')
    end
end

formatplot2
%% Additional functions:

drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );
