%   res1 = "2.41";
%   res2 = "4.46";
%   res3 = "3.49";

  res1 = "2.41";
  res2 = "4.46";
  res3 = "3.50";
  res4 = "7.53";
  res1Ndx = database.findResidue(res1);
  res2Ndx = database.findResidue(res2); % 6.39 or 6.42 for class B2?
  res3Ndx = database.findResidue(res3); %7.56 for class B2?
  res4Ndx = database.findResidue(res4); 

  dis1 = database.calcDistance(res1Ndx, res2Ndx);
  dis2 = database.calcDistance(res3Ndx, res4Ndx);
  paperColors;
  %%
figure('Position',[488.0000   84.2000  804.2000  677.8000])
errStrideRMSD = 50;
% factor to transform from frames to nanoseconds:
timeFactor = 10;
for j = 1:2 % TM36 and TM37 distances
subplot(2,1,j)

    for i =1:2
        if j == 1
            RMSD_mean= mean([dis1{i}{:}],2);
            RMSD_std= std([dis1{i}{:}],0,2);
            runCount = length(dis1{i});
            ylabelText = sprintf("%s-%s distance",res1,res2);
            lineRef1 = dis1{3};
            lineRef2 = dis1{4};
        elseif j==2
            RMSD_mean= mean([dis2{i}{:}],2);
            RMSD_std= std([dis2{i}{:}],0,2);
            runCount = length(dis2{i});
            ylabelText = sprintf("%s-%s distance",res3,res4);
            lineRef1 = dis2{3};
            lineRef2 = dis2{4};
        end
        plot((1:length(RMSD_mean))/timeFactor,RMSD_mean,'LineWidth',1,'Color',paperColorPalette.Hex(2*i));
        hold on
        errorbar((1:errStrideRMSD:length(RMSD_mean))/timeFactor,RMSD_mean(1:errStrideRMSD:length(RMSD_mean)), ...
            RMSD_std(1:errStrideRMSD:length(RMSD_mean))/sqrt(runCount),'o','MarkerSize',5,'Color',paperColorPalette.Hex(2*i-1))
    end
    yline(lineRef1, '--','Active reference','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
    yline(lineRef2, '--','Inactive reference','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
    legend(settings.mainName,'',settings.refName,'',settings.pdbCode, settings.pdbCodeInactive)
    legend boxoff
    xlabel('Time (ns)'); ylabel(ylabelText)
%     title(paramsWT.(fieldsHere{j}).label)
    set(gca,'FontSize',16)
end
sgtitle(sprintf("%s %s,%s-%s,%s-%s distance",settings.systemName,settings.mainName,res1,res2,res3,res4),'fontsize',20)

fileName = sprintf("%s_%s,%s-%s,%s-%s distance",settings.systemName,settings.mainName,res1,res2,res3,res4);
savefig(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\figures\I446N',fileName + ".fig"))
print2pdf(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\figures\I446N',fileName + ".pdf"))

%% Do some dihedral time series too:
resHere = 239;
    % Get dihedrals belonging to residue:
    ndxHere = find(reSortCommon(:,1) == resHere);
    ndxRef = find(reSortCommonRef(:,1) == resHere);

    dihedralsTest(:,ndxHere(j))

    dihedralsRef(:,ndxHere(j))