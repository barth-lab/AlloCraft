paramsF217M = load('D:\Telework_library\dopamine_phase_3\23-d2_F217M_Act2Ina\md2pathdev\workspace.mat','params');
paramsF217M = paramsF217M.params;

paramsWT = load('D:\Telework_library\dopamine_phase_3\25-d2_WT_ICL_pre6VMS\md2pathdev\workspace.mat','params');
paramsWT = paramsWT.params;
paperColors
%% Plots plots plots

figure; 
errStrideRMSD = 50;
fieldsHere = fields(paramsWT);
% factor to transform from frames to nanoseconds:
timeFactor = 10;
for j = 2:length(fieldsHere) % TM36 and TM37 distances
subplot(2,1,j-1)

    for i =1:2
        if i == 1
            RMSD_mean= mean([paramsWT.(fieldsHere{j}).data{1}{:}],2);
            RMSD_std= std([paramsWT.(fieldsHere{j}).data{1}{:}],0,2);
            runCount = length(paramsWT.(fieldsHere{j}).data{1});
        elseif i==2
            RMSD_mean= mean([paramsF217M.(fieldsHere{j}).data{1}{:}],2);
            RMSD_std= std([paramsF217M.(fieldsHere{j}).data{1}{:}],0,2);
            runCount = length(paramsF217M.(fieldsHere{j}).data{1});
        end
        plot((1:length(RMSD_mean))/timeFactor,RMSD_mean,'LineWidth',1,'Color',paperColorPalette.Hex(2*i));
        hold on
        errorbar((1:errStrideRMSD:length(RMSD_mean))/timeFactor,RMSD_mean(1:errStrideRMSD:length(RMSD_mean)), ...
            RMSD_std(1:errStrideRMSD:length(RMSD_mean))/sqrt(runCount),'o','MarkerSize',5,'Color',paperColorPalette.Hex(2*i-1))
    end
    yline(paramsWT.(fieldsHere{j}).data{2}, '--','Active reference','LineWidth',1.5,'Color',[0 0.4470 0.7410]);
    yline(paramsWT.(fieldsHere{j}).data{3}, '--','Inactive reference','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
    legend('D2 WT apo','','D2 F6.44M apo','','6VMS','6CM4')
    legend boxoff
    xlabel('Time (ns)'); ylabel(paramsWT.(fieldsHere{j}).label)
    title(paramsWT.(fieldsHere{j}).label)
    set(gca,'FontSize',16)
end