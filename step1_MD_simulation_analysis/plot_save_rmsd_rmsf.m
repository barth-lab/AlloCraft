% Variables needed in this script:
% RMSD, RMSF
% RMSD_L, RMSF_L
% RMSD_G, RMSF_G

errStrideRMSD = ceil(mainSim.trajLength/60);

% Plot RMSD of receptor, ligand, and effector protein
figure

RMSD_mean= mean(RMSD,2);
RMSD_std= std(RMSD,0,2);


plot(RMSD_mean,'LineWidth',1);
hold on
errorbar(1:errStrideRMSD:length(RMSD_mean),RMSD_mean(1:errStrideRMSD:length(RMSD_mean)), ...
    RMSD_std(1:errStrideRMSD:length(RMSD_mean))/sqrt(mainSim.runCount),'o','MarkerSize',5)

legend_entries{1} = name;
legend_entries{2} = '';

% Ligand (if it exists)
legend_count = 3;
if ~isempty(mainEntry.chains{Chains.ligand})
    RMSD_mean_L= mean(RMSD_L,2);
    RMSD_std_L= std(RMSD_L,0,2);
    plot(RMSD_mean_L,'LineWidth',1);
    hold on
    errorbar(1:errStrideRMSD:length(RMSD_mean_L),RMSD_mean_L(1:errStrideRMSD:length(RMSD_mean_L)),...
        RMSD_std_L(1:errStrideRMSD:length(RMSD_mean_L))/sqrt(mainSim.runCount),'o','MarkerSize',5)
    legend_entries{legend_count} = [name ' ligand'];
    legend_count = legend_count + 1;
    legend_entries{legend_count} = '';
    legend_count = legend_count + 1;
end

% Effector protein (if it exists)
if ~isempty(mainEntry.chains{Chains.gprotein})
    RMSD_mean_G= mean(RMSD_G,2);
    RMSD_std_G= std(RMSD_G,0,2);
    yyaxis right
    plot(RMSD_mean_G,'LineWidth',1);
    hold on
    errorbar(1:errStrideRMSD:length(RMSD_mean_G),RMSD_mean_G(1:errStrideRMSD:length(RMSD_mean_G)),...
        RMSD_std_L(1:errStrideRMSD:length(RMSD_mean_G))/sqrt(mainSim.runCount),'o','MarkerSize',5)
    legend_entries{legend_count} = [name ' effector'];
    legend_count = legend_count + 1;
    legend_entries{legend_count} = '';
    legend_count = legend_count + 1;
end

% Add line to show where we took the equilibration cutoff
axisLimits = axis;
line([settings.frames2skip settings.frames2skip],axisLimits(3:4),'Color','red','LineStyle','--')
legend_entries{legend_count} = 'Equilibration cutoff';
% Add annotations
xlabel('Frame'); ylabel('RMSD [Angstrom]')
legend(legend_entries,'location','best');
legend boxoff
title('RMSD of C\alpha atoms')
set(gca,'FontSize',16)

% save figure as pdf and as .fig
figPath = fullfile(md2pathdir, "RMSD_receptor_ligand");
savefig(figPath);
print2pdf(figPath);

add2log(md2pathdir, 'RMSD plotted and saved successfully');


%% Plot RMSF too!

% Receptor
figure

RMSF_mean_CA = mean(RMSF,2);
RMSF_std_CA = std(RMSF,0,2);

plot(RMSF_mean_CA,'LineWidth',1)
hold on
errorbar(1:errStrideRMSF:length(RMSF_mean_CA),RMSF_mean_CA(1:errStrideRMSF:length(RMSF_mean_CA)),RMSF_std_CA(1:errStrideRMSF:length(RMSF_mean_CA))/sqrt(mainSim.runCount) ...
,'o','MarkerSize',5)
legend_entries{1} =name;
legend_entries{2} = '';

% Draw helices if helices are defined
if ~isempty(settings.helices)
    drawTMhelices(RMSF_mean_CA, settings.helices, mainChain.resIds)
end
xlabel('Residue'); ylabel('RMSF [Angstrom]')
% Add gaps between the legend entries

legend(legend_entries,'location','best');
legend boxoff
title('RMSF of C\alpha atoms - Receptor')
set(gca,'FontSize',16)

% save figure as pdf and as .fig
figPath = fullfile(md2pathdir, "RMSF_receptor");
savefig(figPath);
print2pdf(figPath);

add2log(md2pathdir, 'Receptor RMSF plotted and saved successfully');

% Ligand
if ~isempty(mainEntry.chains{Chains.ligand})
    figure

    RMSF_mean_LHA = mean(RMSF_L,2);
    RMSF_std_LHA = std(RMSF_L,0,2);

    plot(RMSF_mean_LHA,'LineWidth',1)
    hold on


    if ligandChain.isSmall
        xlabel('Ligand heavy atom');
        title('RMSF of heavy atoms - Ligand')
        errStrideRMSF = 1;
    else
        xlabel('Ligand C\alpha');
        title('RMSF of C\alpha - Ligand')
        if length(RMSF_mean_LHA) < 30 % ligand has "small" number of residues
            errStrideRMSF = 1;
        end
    end

    errorbar(1:errStrideRMSF:length(RMSF_mean_LHA),RMSF_mean_LHA(1:errStrideRMSF:length(RMSF_mean_LHA)),RMSF_std_LHA(1:errStrideRMSF:length(RMSF_mean_LHA))/sqrt(mainSim.runCount) ...
    ,'o','MarkerSize',5)
    legend_entries{1} =name;
    legend_entries{2} = '';

    ylabel('RMSF [Angstrom]')
    % Add gaps between the legend entries

    legend(legend_entries,'location','best');
    legend boxoff
    set(gca,'FontSize',16)

    % save figure as pdf and as .fig
    figPath = fullfile(md2pathdir, "RMSF_ligand");
    savefig(figPath);
    print2pdf(figPath);
    add2log(md2pathdir, {'Do not forget ligand RMSF! It has also been plotted and saved',''});
end

% Effector
if mainEntry.hasChain(Chains.gprotein)
    figure
    
    RMSF_mean_G = mean(RMSF_G,2);
    RMSF_std_G = std(RMSF_G,0,2);
    
    plot(RMSF_mean_G,'LineWidth',1)
    hold on
    errorbar(1:errStrideRMSF:length(RMSF_mean_G),RMSF_mean_G(1:errStrideRMSF:length(RMSF_mean_G)),RMSF_std_G(1:errStrideRMSF:length(RMSF_mean_G))/sqrt(mainSim.runCount) ...
    ,'o','MarkerSize',5)
    legend_entries{1} =name;
    legend_entries{2} = '';
    
    xlabel('Residue'); ylabel('RMSF [Angstrom]')
    % Add gaps between the legend entries
    
    legend(legend_entries,'location','best');
    legend boxoff
    title('RMSF of C\alpha atoms - Effector')
    set(gca,'FontSize',16)
    
    % save figure as pdf and as .fig
    figPath = fullfile(md2pathdir, "RMSF_effector");
    savefig(figPath);
    print2pdf(figPath);
    
    add2log(md2pathdir, 'Effector RMSF plotted and saved successfully');
end