%% Visualization of contact map
% Variables needed in this script:
% contact_lig: contact matrix between protein and ligand
% protRes: list of residues of protein according to PDB file
% noH index

% Contact map containing all non-zero elements


%%

function importantResIds = plot_save_contact_map(contactsPerRun, receptorChain, ligandChain, ligandAtomIndices, options)
    arguments
        contactsPerRun
        receptorChain
        ligandChain
        ligandAtomIndices
        options.Directory
        options.ImportantCutoff
        options.RCut
        options.SaveDir
        options.SaveName = "%s"
        options.SaveXlsName
        options.SaveXlsSheet = "Sheet 1"
    end

    contacts = mean(contactsPerRun, 3);
    entry = receptorChain.entry;

    [contactsPruned, row, col] = pruneArray(contacts);
    receptorResIds = receptorChain.resIds(col);

    importantResIds = receptorResIds(max(contactsPruned) > options.ImportantCutoff);

    FigH = figure('Position', get(0, 'Screensize'));
    imagesc(contactsPruned);

    ax = gca;
    fontsz = 12;
    Ang = char(197);

    hcb = colorbar;
    hcb.Label.FontSize = fontsz;
    hcb.Label.String = "Contact frequency";

    ligandName = ligandChain.formatName();

    % TODO: Update to use sprintf()
    title("Residues within " + num2str(options.RCut) + Ang + " of " + ligandName + ...
        ", average over " + num2str(size(contactsPerRun, 3)) + " runs", 'FontSize', fontsz);
    

    if ligandChain.isSmall
        ligandNames = ligandChain.formatAtoms(ligandAtomIndices(row));
        ylabel("Ligand atom", 'FontSize', fontsz);
    else
        ligandNames = ligandChain.formatResidues(ligandChain.resIds(row));
        ylabel("Target residue", 'FontSize', fontsz);
    end
    set(ax, 'FontSize', 16);
    
    receptorNames = receptorChain.formatResidues(receptorResIds);

    xlabel("Protein residue", 'FontSize', fontsz);
    xticks(1:length(receptorNames));
    xticklabels(receptorNames);
    xtickangle(45);

    yticks(1:length(ligandNames));
    yticklabels(ligandNames);


    if length(receptorNames) > 45
        ax.XAxis.FontSize = 12;
    else
        ax.XAxis.FontSize = 16;
    end

    if length(ligandNames) > 30
        ax.YAxis.FontSize = 12;
    end


    if isfield(options, 'SaveDir')
        % Write out ligand contacting residues to file
        writematrix(importantResIds', fullfile(options.SaveDir, sprintf(options.SaveName, "BS_residues")), 'Delimiter', 'space');

        % Save figure as pdf and as .fig
        figPath = fullfile(options.SaveDir, sprintf(options.SaveName, "contact_map_" + ligandName));
        savefig(figPath);
        print2pdf(figPath, 1);

        %% Save ligand contacting data to excel file
        if isfield(options, 'SaveXlsName')
            xlsName = options.SaveXlsName;
        else
            xlsName = options.SaveName;
        end

        xlsPath = fullfile(options.SaveDir, sprintf(xlsName, "ligandContactData") + ".xls");

        writematrix(ligandNames, xlsPath, 'Range', 'A2', 'Sheet', options.SaveXlsSheet);
        writematrix(contactsPruned, xlsPath, 'Range', 'B2', 'Sheet', options.SaveXlsSheet);
        writematrix(receptorNames', xlsPath, 'Range', 'B1', 'Sheet', options.SaveXlsSheet);
    end
end
