% Helix 5 shit:
a = [694	0.402552204
2224	1.201512696
776	0.681898067
192	0.20622986
802	1
589	0.512619669
531	0.820710974
682	0.868789809
719	0.971621622
696	1
368	0.447144593
907	1.098062954
176	0.126618705
561	0.899038462
496	1];

aLabels = ["TM5 major sum" "Ratio (strongest)"];

b = [0	0	0
-13.26	-1.661	-13.26
6.69	-0.394	6.69
4.49	0.285	4.49
9.59	-1.281	9.59
-8.15	-0.852	-8.15
-18.82	0.532	-18.82
0	0	0
-10.45	0.022	-18.29161561
25.95	-0.193	45.42272011
22.59	0.121	39.54139681
5.43	-0.057	9.504638544
-57.13	nan	-100
4	0.162	7.001575354
-53.842	-0.551	-94.24470506];

bLabels = ["dEfficacy" "dlogEC50" "dEfficacy (Normalized)"];

paperColors;
sheetNb = 3; % 3 for D2, 5 for D1

tabFeatures = readtable(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\Compiled_AlloDy_data','features_exp_comparison.xlsx'),'Sheet',sheetNb);
a = [ tabFeatures.Sum tabFeatures.Ratio_strongest_ ];
b = [tabFeatures.dEfficacy tabFeatures.dLogEC50 tabFeatures.dEfficacy_normalized ...
    tabFeatures.dEfficacy_SEM tabFeatures.dLogEC50_SEM tabFeatures.dEfficacy_normalized_SEM];
%% Plotting and regression: combined DA and BRC
xDatNdx = 1;
figure
count = 1;
for yDatNdx = 2:3
    subplot(1,2,count)
    count = count + 1;
    x = a(:,xDatNdx);
    y = b(:,yDatNdx);
    %replace nan with 0:
    y(isnan(y)) = 0;
    X = [ones(length(x),1) x];
    r = X\y;
    
    
    scatter(x(contains(tabFeatures.Variant,"BRC")),y(contains(tabFeatures.Variant,"BRC")),50,'MarkerFaceColor',paperColorPalette.Hex(2),'MarkerEdgeColor',paperColorPalette.Hex(1))
    hold on
    scatter(x(contains(tabFeatures.Variant,"DA")),y(contains(tabFeatures.Variant,"DA")),50,'MarkerFaceColor',paperColorPalette.Hex(4),'MarkerEdgeColor',paperColorPalette.Hex(3))
    
    yCalc2 = X*r;
    plot(x,yCalc2,'-','Color',paperColorPalette.Hex(3),'LineWidth',1.5)
    Rsq = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
    
    xlabel(aLabels(xDatNdx));
    ylabel(bLabels(yDatNdx));
     legend('Data BRC','Data DA',['Fit, R^2=' num2str(Rsq)],'Location','best');
     set(gca,'FontSize',16)

end

%% Separate DA and BRC:
nameHere = 'DA';
rowNdx = (contains(tabFeatures.Variant,nameHere));
xDatNdx = 1;
figure
count = 1;
for yDatNdx = 2:3
    subplot(1,2,count)
    count = count + 1;
    x = a(rowNdx,xDatNdx);
    y = b(rowNdx,yDatNdx);
    yErr = b(rowNdx,yDatNdx+3);
    xRange = range(x);
    yRange = range(y);
    %replace nan with 0:
    y(isnan(y)) = 0;
    X = [ones(length(x),1) x];
    r = X\y;
    
    if strcmp(nameHere,'BRC')
    scatter(x,y,50,'MarkerFaceColor',paperColorPalette.Hex(2),'MarkerEdgeColor',paperColorPalette.Hex(1))
    hold on
    % Add error bars:
    yneg = yErr;
    ypos = yneg;
    xneg = zeros(length(y),1);
    xpos = xneg;
    errorbar(x,y,yneg,ypos,xneg,xpos,"LineStyle","none","Marker","none","Color",paperColorPalette.Hex(1))

    elseif strcmp(nameHere,'DA')
    scatter(x,y,50,'MarkerFaceColor',paperColorPalette.Hex(4),'MarkerEdgeColor',paperColorPalette.Hex(3))
    hold on
        % Add error bars:
    yneg = yErr;
    ypos = yneg;
    xneg = zeros(length(y),1);
    xpos = xneg;
    errorbar(x,y,yneg,ypos,xneg,xpos,"LineStyle","none","Marker","none","Color",paperColorPalette.Hex(3))
    end
    text(x+xRange/100,y+yRange/50, tabFeatures.Variant(rowNdx),'FontSize',10);  
     hold on
%     scatter(x(contains(tabFeatures.Variant,"DA")),y(contains(tabFeatures.Variant,"DA")),50,'MarkerFaceColor',paperColorPalette.Hex(4),'MarkerEdgeColor',paperColorPalette.Hex(3))
%     
    yCalc2 = X*r;
    plot(x,yCalc2,'-','Color',paperColorPalette.Hex(3),'LineWidth',1.5)
    Rsq = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
    
    xlabel(aLabels(xDatNdx));
    ylabel(bLabels(yDatNdx));
    legend(['Data ' nameHere], '' ,['Fit, R^2=' num2str(Rsq)],'Location','best');
    set(gca,'FontSize',20)
%     view([90 -90])
end


%% Separate DA and BRC: Flip x and y axes
nameHere = 'DA';
rowNdx = (contains(tabFeatures.Variant,nameHere));
xDatNdx = 1;
labelNames = tabFeatures.Variant(rowNdx);
labelNamesShort = extractBetween(labelNames,nameHere + "-","-Gi");
    
    
figure
count = 1;
for yDatNdx = 2:3
    subplot(1,2,count)
    count = count + 1;
    x = a(rowNdx,xDatNdx);
    y = b(rowNdx,yDatNdx);
    yErr = b(rowNdx,yDatNdx+3);
    xRange = range(x);
    yRange = range(y);
    %replace nan with 0:
    y(isnan(y)) = 0;
    X = [ones(length(x),1) x];
    r = X\y;
    
    if strcmp(nameHere,'BRC')
    scatter(y,x,100,'MarkerFaceColor',paperColorPalette.Hex(2),'MarkerEdgeColor',paperColorPalette.Hex(1))
    hold on
    % Add error bars:
    yneg = yErr;
    ypos = yneg;
    xneg = zeros(length(y),1);
    xpos = xneg;
    errorbar(y,x,xneg,xpos,yneg,ypos,"LineStyle","none","Marker","none","Color",paperColorPalette.Hex(1))

    elseif strcmp(nameHere,'DA')
    scatter(y,x,100,'MarkerFaceColor',paperColorPalette.Hex(4),'MarkerEdgeColor',paperColorPalette.Hex(3))
    hold on
        % Add error bars:
    yneg = yErr;
    ypos = yneg;
    xneg = zeros(length(y),1);
    xpos = xneg;
    errorbar(y,x,xneg,xpos,yneg,ypos,"LineStyle","none","Marker","none","Color",paperColorPalette.Hex(3))
    end
%     text(y+yRange/50,x+xRange/30, labelNamesShort,'FontSize',10);  
     hold on
%     scatter(x(contains(tabFeatures.Variant,"DA")),y(contains(tabFeatures.Variant,"DA")),50,'MarkerFaceColor',paperColorPalette.Hex(4),'MarkerEdgeColor',paperColorPalette.Hex(3))
%     
    yCalc2 = X*r;
    plot(yCalc2,x,'-','Color',paperColorPalette.Hex(3),'LineWidth',1.5)
    Rsq = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2);
    
    ylabel(aLabels(xDatNdx));
    xlabel(bLabels(yDatNdx));
    legend(['Data ' nameHere], '' ,['Fit, R^2=' num2str(Rsq)],'Location','best');
    set(gca,'FontSize',20)
%     view([90 -90])

end