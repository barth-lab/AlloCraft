 %% This script will calculate the g-protein contacting residues from the
 % downloaded pdb: pdbCode and then translate that into input PDB residue
 % numbering
 %
 % Variables needed in this script:
 % Input PDB obviously
 % protRes: list of residues of protein according to PDB file
 %


%%


function [receptorGpResIdsInput, receptorGpResIdsRef] = calc_plot_gp_contact_map(inputEntry, refEntry, Chains, options)
    arguments
        inputEntry
        refEntry
        Chains
        options.SaveDir = []
        options.SaveName = "%s"
        options.StartFrame
    end

    saveName = sprintf(options.SaveName, "contact_map_" + inputEntry.name + "_G-protein");

    database = inputEntry.database;
    receptorResidues = database.residues{Chains.receptor};
    gproteinResidues = database.residues{Chains.gprotein};


    % (1) Simulation

    if inputEntry.hasChain(Chains.gprotein)
        receptorChain = inputEntry.chains{Chains.receptor};
        gproteinChain = inputEntry.chains{Chains.gprotein};

        contacts = mean(inputEntry.simulation.computeContacts(receptorChain, gproteinChain, 'StartFrame', options.StartFrame), 3);
        plotContacts(contacts, "Contact map of simulation", receptorChain, gproteinChain, 'SaveDir', options.SaveDir, 'SaveName', saveName + "_sim");

        if ~isempty(refEntry)
            receptorRes = receptorResidues(receptorResidues{:, inputEntry.databaseIndex} > 0, :);
            gproteinRes = gproteinResidues(gproteinResidues{:, inputEntry.databaseIndex} > 0, :);

            inp.subset = contacts(gproteinRes{:, refEntry.databaseIndex} > 0, receptorRes{:, refEntry.databaseIndex} > 0);
        end

        [~, ~, col] = pruneArray(contacts);
        receptorGpResIdsInput = receptorChain.resIds(col);
    else
        receptorGpResIdsInput = [];
    end


    % (2) Reference

    if ~isempty(refEntry)
        receptorChain = refEntry.chains{Chains.receptor};
        gproteinChain = refEntry.chains{Chains.gprotein};

        contacts = refEntry.computeContacts(receptorChain, gproteinChain);
        %plotContacts(contacts, "Contact map of reference", receptorChain, gproteinChain, 'SaveDir', options.SaveDir, 'SaveName', saveName + "_ref");
        contacts = contacts(:, 1:length(receptorChain.resIds));

        receptorRes = receptorResidues(receptorResidues{:, refEntry.databaseIndex} > 0, :);
        gproteinRes = gproteinResidues(gproteinResidues{:, refEntry.databaseIndex} > 0, :);

        if inputEntry.hasChain(Chains.gprotein)
            ref.subset = contacts(gproteinRes{:, inputEntry.databaseIndex} > 0, receptorRes{:, inputEntry.databaseIndex} > 0);
        end

        [~, ~, col] = pruneArray(contacts);
        col(col>size(receptorRes,1))=[]; % Clean any contacts that are with
        % residues not found in the input pdb
        receptorGpResIdsRef = receptorRes{col, inputEntry.databaseIndex};
    else
        receptorGpResIdsRef = [];
    end


    % (3) Difference

    if inputEntry.hasChain(Chains.gprotein) && ~isempty(refEntry)
        receptorResIds = receptorResidues{(receptorResidues{:, inputEntry.databaseIndex} > 0) & (receptorResidues{:, refEntry.databaseIndex} > 0), inputEntry.databaseIndex};
        gproteinResIds = gproteinResidues{(gproteinResidues{:, inputEntry.databaseIndex} > 0) & (gproteinResidues{:, refEntry.databaseIndex} > 0), inputEntry.databaseIndex};

        diff = inp.subset - ref.subset;

        receptorChain = inputEntry.chains{Chains.receptor};
        gproteinChain = inputEntry.chains{Chains.gprotein};

        plotContacts(diff, "Difference between contact maps of simulation and reference", receptorChain, gproteinChain, receptorResIds, gproteinResIds, 'SaveDir', options.SaveDir, 'SaveName', saveName + "_diff");
    end
end


function plotContacts(contacts, figTitle, xChain, yChain, xResIds, yResIds, options)
    arguments
        contacts
        figTitle
        xChain
        yChain
        xResIds = xChain.resIds
        yResIds = yChain.resIds
        options.Name
        options.SaveDir
        options.SaveName
    end

    database = xChain.entry.database;
    [contactsPruned, row, col] = pruneArray(contacts);

    figure('Position', get(0, 'Screensize'));
    imagesc(contactsPruned);
    title(figTitle);
    colorbar;

    % FigH = figure('Position', get(0, 'Screensize'));

    xLabels = xChain.formatResidues(xResIds(col));
    yLabels = yChain.formatResidues(yResIds(row));

    xticks(1:length(xLabels));
    xticklabels(xLabels);
    xtickangle(45);

    yticks(1:length(yLabels));
    yticklabels(yLabels);

    if (length(database.entries) > 2) && (xChain.entry ~= database.entries{2})
        xlabel("Receptor residues: input residue (ref. residue, BW numbering)");
        ylabel("G protein residues: input residue (ref. residue)");
    else
        xlabel("Receptor residues: input residue");
        ylabel("G protein residues");
    end

    set(gca, 'FontSize', 16);

    ax = gca;
    ax.XAxis.FontSize = 12;

    if length(yLabels) > 30
        ax.YAxis.FontSize = 12;
    end

    if ~isempty(options.SaveDir)
        figPath = fullfile(options.SaveDir, options.SaveName+ ".fig");
        savefig(figPath);
        print2pdf(figPath, 1);
    end
end
