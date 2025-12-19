classdef Simulation < handle
    properties
        entry
        traj
        trajLength

        dihedrals
        dihedralsMat
        mi
        reSort
        klDivStuff
    end

    properties (Dependent)
        runCount
    end

    methods
        function obj = Simulation(entry)
            global settings;
            obj.entry = entry;
        end

        function value = get.runCount(obj)
            value = length(obj.traj);
        end

        function align(obj,options)
            arguments
                obj
                options.chain = 'all'
            end
            if strcmp(options.chain, 'all')
                index_CA = selectname(obj.entry.pdb.name, 'CA');
            else
                index_CA = selectname(obj.entry.pdb.name, 'CA') & selectname(obj.entry.pdb.chainid,options.chain);
            end
            traj_aligned = cell(max(obj.runCount),1);

            for runIndex = 1:obj.runCount
                [~, traj_aligned{runIndex}] = superimpose(obj.entry.crd, obj.traj{runIndex}, find(index_CA));
            end

            obj.traj = traj_aligned;
        end

        function traj = concatRuns(obj, options)
            arguments
                obj
                options.Atoms = 1:obj.entry.atomCount
                options.StartFrame = 1
            end

            runLengths = cellfun(@(x) size(x, 1), obj.traj) - options.StartFrame + 1;
            runStartFrames = [1; cumsum(runLengths) + 1];
            xyzIndices = to3(options.Atoms);
            traj = zeros(sum(runLengths), length(xyzIndices));

            for runIndex = 1:obj.runCount
                runStartFrame = runStartFrames(runIndex);
                traj(runStartFrame:(runStartFrame + runLengths(runIndex) - 1), :) = obj.traj{runIndex}(options.StartFrame:end, xyzIndices);
            end
        end


        function computeDihedrals(obj, chain, options)
            arguments
                obj
                chain = [] % Chain object
                options.HigherOrder = 'all'
                options.Path
                options.ReSortPath
                options.ResIds = chain.resIds
                options.StartFrame = 1
            end

            if isfield(options, 'Path') && exist(options.Path, 'file')
                load(options.Path, 'dihedrals', 'reSort');
            else
                concatTraj = obj.concatRuns('StartFrame', options.StartFrame);
                if ~isempty(chain) % Take dihedrals only from input chain
                    [indices, ~] =  obj.entry.getAtoms('Chain',chain.index);
                    pdbChain = substruct(obj.entry.pdb, indices);
                    [dihedrals, ~, reSort] = calcalldihedralsfromtrajs(pdbChain, concatTraj(:,to3(indices)), options.ResIds, 1, options.HigherOrder);
                else
                    [dihedrals, ~, reSort] = calcalldihedralsfromtrajs(obj.entry.pdb, concatTraj, options.ResIds, 1, options.HigherOrder);
                end

                if isfield(options, 'Path')
                    save(options.Path, 'dihedrals', 'reSort');
                end

                if isfield(options, 'ReSortPath')
                    fd = fopen(options.ReSortPath, 'w');

                    for i = 1: length(reSort)
                        fprintf(fd, '%d %d %d %d\r\n', reSort(i, [1 2 4 5]));
                    end

                    fclose(fd);
                end
            end

            obj.dihedrals = dihedrals;
            obj.reSort = reSort;


            % Compute the dihedrals matrix

            width = sum(cellfun(@(x) size(x, 2), obj.dihedrals));

            for res = 1:length(obj.dihedrals)
                resDihedrals = obj.dihedrals{res};

                if isnan(resDihedrals(1, 1)) || isnan(resDihedrals(1, 2))
                    width = width - 1;
                end
            end

            obj.dihedralsMat = zeros(length(obj.dihedrals{1}), width);

            dih = 1;

            for res = 1:length(obj.dihedrals)
                resDihedrals = obj.dihedrals{res};
                resDihedralsCount = size(resDihedrals, 2);

                offsetStart = isnan(resDihedrals(1, 1));
                offsetEnd = isnan(resDihedrals(1, 2));

                if offsetStart, resDihedrals = resDihedrals(:, 2:end); end
                if offsetEnd, resDihedrals = resDihedrals(:, [1 3:end]); end

                obj.dihedralsMat(:, dih:(dih + resDihedralsCount - offsetStart - offsetEnd - 1)) = resDihedrals;
                dih = dih + resDihedralsCount - offsetStart - offsetEnd;
            end

            obj.dihedralsMat = real(obj.dihedralsMat);
        end


        function computeMI(obj, options)
            arguments
                obj
                options.Path
            end

            if isfield(options, 'Path') && exist(options.Path, 'file')
                load(options.Path, 'I');
            else
                [~, ~, I] = computeEntropies(zeroStretchtotwopi(obj.dihedralsMat));
                I = I + I';

                save(options.Path, 'I');
            end

            obj.mi = I;
        end

        function [reSortCommon, reSortCommonRef] = reconcileDihedralList(obj,refSim,options) 
        % Find the union of the two dihedral lists to input into the KL 
        % divergence function
            arguments
                    obj
                    refSim
                    options.Path
                    options.refEntryNdx = 2
             end
            nRes = length(refSim.dihedrals); % Effective number of residues            

%             reSortndx = true(size(refSim.reSort,1),1);
%             reSortndxTest = true(size(obj.reSort,1),1);
              reSortndx = false(size(refSim.reSort,1),1);
              reSortndxTest = false(size(obj.reSort,1),1);

            for thisRes = 1:nRes

                % Find residue number in obj (test sim), assumes obj is
                % first column in database
                thisResTest = obj.entry.database.residues{1}.(1)(obj.entry.database.residues{1}.(options.refEntryNdx) == thisRes,1);
               
                % Indices for this residue in reSort
                tempRef = find(refSim.reSort(:,1)==thisRes);
                temp = find(obj.reSort(:,1)==thisResTest);
                
                % True or false for this dihedral
                tfRef = true(length(tempRef),1);
                tfTest = true(length(temp),1);
                
                % Dihedral types
                avDihsTest = obj.reSort(obj.reSort(:,1)==thisResTest,2);
                avDihsRef = refSim.reSort(refSim.reSort(:,1)==thisRes,2);
                
                % BB1
                if sum(avDihsTest==1) ~= sum(avDihsRef==1) 
                    tfTest(avDihsTest==1) = false;
                    tfRef(avDihsRef==1) = false;
                end
                % BB2
                if sum(avDihsTest==2) ~= sum(avDihsRef==2) 
                    tfTest(avDihsTest==2) = false;
                    tfRef(avDihsRef==2) = false;
                end
                % SC
                if sum(avDihsRef==0) > sum(avDihsTest==0)  % Ref structure has more SC dihedrals on this residue
                
                    diff = sum(avDihsRef==0) - sum(avDihsTest==0);
                    tfRef((end - diff + 1):end) = false;
                
                elseif sum(avDihsTest==0) > sum(avDihsRef==0) % Test structure has more SC dihedrals on this residue
                    diff = sum(avDihsTest==0) - sum(avDihsRef==0);
                    tfTest((end - diff + 1):end) = false;     
                end
                reSortndx(tempRef) = tfRef;
                reSortndxTest(temp) = tfTest;

%                 nDihHere = [ sum(refSim.reSort(:,1)==thisRes) sum(obj.reSort(:,1)==thisResTest)];         
%                 if nDihHere(2) > nDihHere(1) % Test structure has more dihedrals on this residue
% 
% 
%                     diffHere = nDihHere(2) - nDihHere(1);
%                     temp = find(obj.reSort(:,1)==thisResTest);
%                     reSortndxTest( temp((end - diffHere + 1):end)) = false;
%                 elseif nDihHere(1) > nDihHere(2) % Other way around, main needs to be cut
%                     diffHere = nDihHere(1) - nDihHere(2);
%                     temp = find(refSim.reSort(:,1)==thisRes);
%                     reSortndx( temp((end - diffHere + 1):end)) = false; 
% 
%                 end
            end
            
             assert( isequal(refSim.reSort(logical(reSortndx),2) ,obj.reSort(logical(reSortndxTest),2)), ...
     'Program was not able to reconcile the dihedral list from main and test structure!')

             reSortCommon = obj.reSort(logical(reSortndxTest),1:2);
             reSortCommonRef = refSim.reSort(logical(reSortndx),1:2);

             obj.klDivStuff.reSortndxHere = reSortndxTest;
             obj.klDivStuff.reSortndxRef = reSortndx;
             save(options.Path, 'reSortCommon', 'reSortCommonRef');
        end

        function [contactsPerRun, chainBAtomIndices, contactResAll] = computeContacts(obj, chainA, chainB, options)
            arguments
                obj
                chainA
                chainB
                options.ContactCut = 0.4
                options.RCut1 = 3.5
                options.RCut2 = 5
                options.StartFrame = 1
            end

            contactsPerRun = 0;
            contactResAll = cell(obj.runCount,1);

            for runIndex = 1:obj.runCount
                [~, contactRes, chainBAtomFilter, ~, runContacts] = proteinContactsDualCutoff(obj.entry.pdb, obj.traj{runIndex}(options.StartFrame:end, :), ...
                    chainB.resIds, char(chainB.name), char(chainA.name), options.RCut1, options.RCut2,options.ContactCut);
                runContacts = runContacts{1};

                if ~contactsPerRun
                    contactsSize = size(runContacts);
                    contactsPerRun = zeros(contactsSize(1), contactsSize(2), obj.runCount);
                end
                contactResAll{runIndex} = contactRes{1};
                contactsPerRun(:, :, runIndex) = runContacts;
            end

            chainBAtomIndices = find(chainBAtomFilter);
        end

        function rmsd = computeRmsd(obj, atomIndices)
            rmsd = zeros(obj.trajLength, obj.runCount);

            for runIndex = 1:obj.runCount
                rmsd(:, runIndex) = calcrmsd(obj.traj{runIndex}(1, :), obj.traj{runIndex}(1:obj.trajLength, :), atomIndices);
            end
        end
        
        function rmsd = computeRmsdAllFramePairs(obj, atomIndices, options)
            arguments
                obj
                atomIndices
                options.StartFrame = 1
            end
            
            C=[];
            nFramesEff = zeros(obj.runCount,1);
            for i = 1:obj.runCount
                % Take into consideration runs with different number of frames
                nFramesEff(i) = (size(obj.traj{i},1) - options.StartFrame)+1; % Frames used in calculations
                framesNdx = ((options.StartFrame):size(obj.traj{i},1))';
                Ctemp = [i*ones(nFramesEff(i),1) framesNdx];
                C=[C ;Ctemp];% Used for coloring and for labeling: [ run frameNdx]
            end
            
            
            trajTemp = obj.concatRuns('Atoms', atomIndices, 'StartFrame', options.StartFrame);
            rmsdMap = zeros(size(trajTemp,1));
            for runIndex = 1:obj.runCount
                rmsd(:, runIndex) = calcrmsd(obj.traj{runIndex}(1, :), obj.traj{runIndex}(1:obj.trajLength, :), atomIndices);
            end
        end

        function rmsf = computeRmsf(obj, atomIndicesInput, options)
            arguments
                obj
                atomIndicesInput
                options.Exact = false
                options.StartFrame = 1
            end

            if options.Exact
                atomIndicesCalc = 1:obj.entry.atomCount;
                atomIndicesFit = atomIndicesInput;
            else
                atomIndicesCalc = atomIndicesInput;
                atomIndicesFit = 1:length(atomIndicesInput);
            end

            rmsf = zeros(length(atomIndicesCalc), obj.runCount);

            for runIndex = 1:obj.runCount
                rmsf(:, runIndex) = calcrmsf(obj.traj{runIndex}(options.StartFrame:obj.trajLength, to3(atomIndicesCalc)), atomIndicesFit);
            end

            if options.Exact
                rmsf = rmsf(atomIndicesInput, :);
            end
        end

        function subset = createSubset(obj, selection, options)
            arguments
                obj
                selection
                options.DihedralMatIndices
            end

            subset = Simulation(obj.entry);

            subset.trajLength = length(selection);
            subset.traj = {zeros(subset.trajLength, obj.entry.atomCount * 3)};

            for frameSelectionIndex = 1:subset.trajLength
                runIndex = selection(frameSelectionIndex, 1);
                frameIndex = selection(frameSelectionIndex, 2);
                subset.traj{1}(frameSelectionIndex, :) = obj.traj{runIndex}(frameIndex, :);
            end

            subset.dihedrals = obj.dihedrals;
            subset.reSort = obj.reSort;

            if ~isempty(obj.dihedralsMat) && isfield(options, 'DihedralMatIndices')
                subset.dihedralsMat = obj.dihedralsMat(options.DihedralMatIndices, :);
            end
        end
    end

    methods(Static)
        function simulation = query(entry, query)
            global settings;
            runDirs = dir(query);
            runCount = length(runDirs);

            simulation = Simulation(entry);
            simulation.traj = cell(runCount, 1);
            maxTraj = 0;

            for runIndex = 1:runCount
                directory = runDirs(runIndex);
                runTraj = readdcdmat(fullfile(directory.folder, directory.name, settings.xtcName));

                if runIndex == 1
                    minTraj = size(runTraj,1);
                end

                minTraj = min(minTraj, size(runTraj, 1));
                maxTraj = max(maxTraj, size(runTraj, 1));

                simulation.traj{runIndex} = runTraj;
            end


            if (maxTraj - minTraj) / minTraj < 0.01
                for runIndex = 1:runCount
                    simulation.traj{runIndex} = simulation.traj{runIndex}(1:minTraj, :);
                end
            else
%                 add2log(md2pathdir,{'Different runs have more than 1% difference in length, may affect RMSD/RMSF calculations.',''});
            end

            simulation.trajLength = minTraj;

            % Make sure PDB and traj has the same number of atoms
            % assert(size(obj.traj{1},2)/3 == length(pdb.serial),'PDB and trajectory have different number of atoms!!!')
        end
    end
end
