function  pathCalcdir = prepareAlloPathCalc(simulation, receptorChain, settings, options)
    arguments
        simulation
        receptorChain
        settings

        options.Dir
        options.Name = "%s"
        options.LogPath

        options.ReceptorLigandResIds
        options.ReceptorGpResIds
        options.ReceptorResIds = receptorChain.resIds

        options.Traj = []

        options.DisCutoff = settings.disCutoff
        options.MIFractionCutoff = settings.miFractionCutoff
        options.NearCutoff = settings.nearCutoff
        options.OverlapCutoff = settings.overlapCutoff
    end    

    pathCalcdir = fullfile(options.Dir, sprintf(options.Name, 'alloPathCalc'));

%     if exist(pathCalcdir, 'dir')
%         error("Directory already exists in " + pathCalcdir);
%     end

    mkdir(pathCalcdir);

    % Renumbered residues

    renumber = @(resIds) arrayfun(@(x) find(options.ReceptorResIds == x), resIds);

    % Write protein only coordinates to pdb and renumber the residues from 1

     writepdbndx(fullfile(pathCalcdir, "protRenum.pdb"), simulation.entry.pdb, selectid(simulation.entry.pdb.resseq, options.ReceptorResIds), true, options.Traj);
%      writepdb(fullfile(pathCalcdir, "prot.pdb"), simulation.entry.pdb,options.Traj)

    add2log(options.LogPath, "Extracted receptor PDB, chain: " + receptorChain.name + ", renumbered from 1, and saved it to protRenum.pdb");
    add2log(options.LogPath, "All files after this point will be written to the path calculation directory: " + pathCalcdir);


    % Excluded residues

    if isempty(settings.excludedResidues)
%         excludedResidues = [options.ReceptorResIds(1) options.ReceptorResIds(4); (options.ReceptorResIds(end) - 3) (options.ReceptorResIds(end))];
        excludedResidues = [options.ReceptorResIds(1) options.ReceptorResIds(4); (options.ReceptorResIds(end) - 3) (options.ReceptorResIds(end))];
        add2log(options.LogPath, 'Added first and last 4 residues as default excluded residues in the pathway calculation');
    else
        excludedResidues = renumber(settings.excludedResidues);
%         excludedResidues = (settings.excludedResidues);
        add2log(options.LogPath, 'Using user defined excluded_residues variable');
    end

    writematrix(excludedResidues, fullfile(pathCalcdir, "excluded_res.txt"), 'Delimiter', 'space');


    % Helices

    if settings.isGPCR % GPCR specific option
        % Helices: extract secodary structure from VMD
        % run VMD from command line:  vmd -dispdev text -e
        if isempty(settings.helices)
            receptorChain.computeSecStruct();
            helices = receptorChain.helices;
            add2log(options.LogPath,['TM helices calculated from extracted secondary structure and written to ' pathCalcdir 'helices.txt'])
        else
            helices = renumber(settings.helices);
%             helices = (settings.helices);
            receptorChain.helices = helices;
            add2log(options.LogPath,['Using user defined helices variable that has been written to ' pathCalcdir 'helices.txt']);
        end

        writematrix(helices, fullfile(pathCalcdir, "helices.txt"), 'Delimiter', 'space');
    end


    % BAI residues

    if isempty(settings.baiResidues)
        baiResidues = [1]; % Placeholder since bai calculations won't matter
        add2log(options.LogPath,'Beta arrestin calculations are ignored');
    else
        baiResidues = renumber(settings.baiResidues);
%         baiResidues = (settings.baiResidues);
        add2log(options.LogPath,'Beta arrestin interacting residues saved to BAI_residues.txt');
    end

    writematrix(baiResidues', fullfile(pathCalcdir, "BAI_residues.txt"), 'Delimiter', 'space'); 


    % Ligand & G protein residues

    liResRenum = renumber(options.ReceptorLigandResIds);
%     liResRenum = (options.ReceptorLigandResIds);
    writematrix(liResRenum', fullfile(pathCalcdir, "BS_residues.txt"), 'Delimiter', 'space'); 
    add2log(options.LogPath,'Binding site residues renumbered and saved to BS_residues.txt');

    % GPI residues (we have them already just transform them to renum)
    gpiResRenum = renumber(options.ReceptorGpResIds);
%     gpiResRenum = (options.ReceptorGpResIds);
    writematrix(gpiResRenum', fullfile(pathCalcdir, "GPI_residues.txt"), 'Delimiter', 'space');
    add2log(options.LogPath,'G-protein interacting residues renumbered and saved to GPI_residues.txt');

    % Pathway analysis binding site residues and EC residues
    if settings.isGPCR
        helices = receptorChain.helices;

        all_BS = liResRenum; %to be filled in before running
        % Use the helices variable to calculate EC: EC is assume to be the top two
        % helical turns (~7 residues) of every helix and EC loops
        all_EC = [1:(helices(1,1)+7) (helices(2,2)-7):(helices(3,1)+7) ...
            (helices(4,2)-7):(helices(5,1)+7) (helices(6,2)-7):(helices(7,1)+7)];
    end


    % Renumber reSort

    % reSort: Make sure reSort is numbered according to protein-renum

    protAtoms = simulation.entry.pdb.serial(selectid(simulation.entry.pdb.resseq, options.ReceptorResIds));

    fileID = fopen(fullfile(pathCalcdir, 'reSort.txt'),'w');
    reSortRenum = zeros(length(simulation.reSort),6);

    for i = 1:length(simulation.reSort)
        atom1 = find(protAtoms == simulation.reSort(i,3));
        atom2 = find(protAtoms == simulation.reSort(i,4));
        atom3 = find(protAtoms == simulation.reSort(i,5));
        atom4 = find(protAtoms == simulation.reSort(i,6));
        reSortRenum(i,:) = [simulation.reSort(i,[1 2]) atom1 atom2 atom3 atom4];
        fprintf(fileID,'%d %d %d %d %d %d\r\n',reSortRenum(i,:));
    end

    fclose(fileID);
    add2log(options.LogPath,'reSort renumbered and saved to reSort.txt');


    % Rewrite MI

    miPath = fullfile(pathCalcdir, "MI_.mat");

    I = simulation.mi;
    save(miPath, 'I');


    % Rewrite dihedrals mat

    dihedralsPath = fullfile(pathCalcdir, "dihedrals.mat");

    dihedralsMat = simulation.dihedralsMat;
    save(dihedralsPath, 'dihedralsMat');


    % Write variables to workspace file:
    variables.dihedralsPath = dihedralsPath;
    variables.disCutOff = options.DisCutoff;
    variables.MIFractionCutoff = options.MIFractionCutoff;
    variables.miName = char(miPath);
    variables.nearcutoff = options.NearCutoff;
    variables.overlapcutoff = options.OverlapCutoff;
    variables.reSort = simulation.reSort;
    variables.includeLigPathwayCalc = settings.includeLigPathwayCalc;
    variables.isGPCR = settings.isGPCR;
    variables.diagnosticsOn = settings.diagnosticsOn;
    variables.MIWeightPaths = settings.MIWeightPaths;
    if settings.isGPCR
        variables.all_BS = all_BS;
        variables.all_EC = all_EC;
        variables.helices = helices;
    end

    save(fullfile(pathCalcdir, "workspace.mat"), '-struct', 'variables');

    add2log(options.LogPath, sprintf("All input is ready, running allosteric path calculation now with the following options:\n" ...
        + "  Distance cutoff for allosteric signal: %f Ang\n" ...
        + "  Near cutoff: %f Ang\n" ...
        + "  Pathway overlap cutoff: %f", options.DisCutoff, options.NearCutoff, options.OverlapCutoff));
end
