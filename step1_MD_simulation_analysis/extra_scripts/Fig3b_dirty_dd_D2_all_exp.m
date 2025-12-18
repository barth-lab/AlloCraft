%% Using labels and alignments: calculating delta delta efficacy (BRC - DA)(mutant - WT)
tabDA1 = readtable(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data','delta_exp_all_differently_formatted'), ...
    'Range','A2:C24','Sheet','DA_dlogEC50');
tabDA2 = readtable(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data','delta_exp_all_differently_formatted'), ...
    'Range','A2:C24','Sheet','DA_dEFF');
tabBRC1 = readtable(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data','delta_exp_all_differently_formatted'), ...
    'Range','A2:C24','Sheet','BRC_dlogEC50');
tabBRC2 = readtable(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data','delta_exp_all_differently_formatted'), ...
    'Range','A2:C24','Sheet','BRC_dEFF');

mutList  = tabDA1.Label;


deltadeltadLogEC50 = tabBRC1.dMean - tabDA1.dMean;
deltadeltadEFF = tabBRC2.dMean - tabDA2.dMean;



deltadeltadLogEC50(isnan(deltadeltadLogEC50))=0;
% table(mutList,deltadelta)

%
pos = [112.2000  109.0000  772.8000  640.0000];
figure('Position',pos)
xRange = range(deltadeltadLogEC50);
yRange = range(deltadeltadEFF);
s =scatter(deltadeltadLogEC50,deltadeltadEFF,100,'LineWidth',1.5,'MarkerFaceColor',paperColorPalette.Hex(8), ...
   'MarkerEdgeColor',paperColorPalette.Hex(7));
hold on
% Mark the WT:
scatter(0,0,200,'x','linewidth',4,'MarkerEdgeColor',paperColorPalette.Hex(9))
xlabel('ddLog(EC50)')
ylabel('ddEfficacy')
row = dataTipTextRow('Label',mutList);
s.DataTipTemplate.DataTipRows(end+1) = row;

ytextPos = deltadeltadEFF+yRange/25;
% ytextPos(4:5) = deltadelta(4:5,2)-yRange/25; % Fix text for 3.41H and 3.41G

textfit(deltadeltadLogEC50+xRange/40,ytextPos, mutList,'FontSize',14);     
% Add error bars:
% yneg = deltadelta(:,4);
% ypos = yneg;
% xneg = deltadelta(:,3);
% xpos = xneg;
% errorbar(deltadelta(:,1),deltadelta(:,2),yneg,ypos,xneg,xpos,"LineStyle","none","Marker","none","Color",paperColorPalette.Hex(7))

title('\Delta\Delta (BRC - DA) (Mut - WT)')
formatplot2