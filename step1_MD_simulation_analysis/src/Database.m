classdef Database < handle
    properties
        dir
        entries
        entryState % defines state of entry: 0 1 2 3 --> (simulation, active, intive, other)
        residues
    end


    methods
        function obj = Database(dir)
            obj.dir = dir;
            obj.entries = cell(1, 0);
            obj.entryState = zeros(1,0);

            if ~exist(obj.dir, 'dir')
                mkdir(obj.dir)
            end
        end

        function fetch(obj, pdbCode, targetChainNames, stateHere)
            if nargin < 4
                stateHere = 3; % assume it's an "other" state
            end
            filename = fullfile(obj.dir, pdbCode + ".pdb");
            % add2log("Fetching PDB " + pdbCode + " to " + filename);

            if ~exist(filename, 'file')
                url = "https://files.rcsb.org/download/" + pdbCode + ".pdb";
                websave(filename, url);
            end

            obj.read(filename, targetChainNames, pdbCode, stateHere);
        end

        function read(obj, path, targetChainNames, name, stateHere)
            if nargin <5
                stateHere = 0; % assume we're inputing simulation data
            end
            entry = Entry(path, targetChainNames, name);
            obj.entries{end + 1} = entry;
            obj.entryState(end + 1) = stateHere;

            entry.database = obj;
            entry.databaseIndex = length(obj.entries);
        end

        function align(obj)
            chainCount = max(cellfun(@(entry) length(entry.chains), obj.entries));
            sequences = cell(chainCount, length(obj.entries));

            for chainIndex = 1:chainCount
                for entryIndex = 1:length(obj.entries)
                    entry = obj.entries{entryIndex};

                    if entry.hasChain(chainIndex)
                        sequences{chainIndex, entryIndex} = char(entry.seq{entry.chainIndices(chainIndex)});
                    end
                end
            end

            for chainIndex = 1:chainCount
                chainSequences = sequences(chainIndex, :);
                aligned = alignSequences(chainSequences);

                residueCount = size(aligned, 2);
                output = zeros(residueCount, length(obj.entries));
                names = repmat(blanks(length(obj.entries)), residueCount, 1);

                for entryIndex = 1:length(obj.entries)
                    entry = obj.entries{entryIndex};

                    if ~entry.hasChain(chainIndex)
                        continue;
                    end

                    chain = entry.chains{chainIndex};

                    entryResIndex = 1;

                    for resIndex = 1:residueCount
                        value = aligned(entryIndex, resIndex);

                        if value ~= '-'
                            output(resIndex, entryIndex) = chain.resIds(entryResIndex);
                            names(resIndex, entryIndex) = value;
                            entryResIndex = entryResIndex + 1;
                        end
                    end
                end

                obj.residues{chainIndex} = array2table(output);
                obj.residues{chainIndex}.Label = repmat("", residueCount, 1);
                obj.residues{chainIndex}.Name = names;
            end
        end

        function obj = label(obj, entryIndex, bwMapPath)
            bwMap = readtable(bwMapPath, 'MissingRule', 'omitrow', 'TextType', 'char', 'ReadVariableNames', false);
            bwMap = bwMap(2:end, 2:3);
            % entry = obj.entries{entryIndex};
            entryColumn = obj.residues{1}{:, entryIndex};
            % entryColumn = entryResIds(obj.residues{:, entryIndex});

            for mapRowIndex = 1:height(bwMap)
                label = string(bwMap{mapRowIndex, 1});
                refIndex = bwMap{mapRowIndex, 2};
                refIndex = refIndex{1};
                refIndex = str2double(refIndex(2:end));
                tableRowIndex = find(entryColumn == refIndex);
                % display(entryColumn(tableRowIndex));

                % resIndex = find(entry.resIds == refIndex);
                % tableRowIndex = find(entryColumn == resIndex);

                if tableRowIndex
                    obj.residues{1}.Label(tableRowIndex) = formatBw(label);
                end
            end
        end
        
        function [tempCell, mutPos, obj] = findMut(obj, entryIndex, chainIndex, fastaPath)
            % Read sequence from database and pdb

            [~, seqRef] = fastaread(fastaPath);
            seqTest = obj.entries{entryIndex}.seq{chainIndex};
             % align and compare
            [~,b] = nwalign(seqRef,seqTest);
            
            tempseq1 = b(1,1:end);
            tempseq2 = b(3,1:end);
            substitutions = tempseq1~=tempseq2 & tempseq1~='-' & tempseq2~='-';
            mutPos = find(substitutions);
            tempCell = cell(length(mutPos),1);
            for i=1:length(mutPos)
                
                tempCell{i} = [tempseq1(mutPos(i)) num2str(mutPos(i)) tempseq2(mutPos(i))];
            end
        end

        function result = findResidue(obj, name, options)
            arguments
                obj
                name
                options.ChainIndex = 1
            end

            rows = obj.residues{options.ChainIndex}(obj.residues{options.ChainIndex}.Label == name, 1:(end - 2));
            result = rows{1, :}';
        end

        function result = findFeature(obj, name, before, after, options)
            arguments
                obj
                name
                before
                after
                options.ChainIndex = 1
            end

            result = cell2mat(arrayfun(@(a) (a - before):(a + after), obj.findResidue(name, 'ChainIndex', options.ChainIndex), 'UniformOutput', false));
        end

        function alignStructures(obj)
            refEntry = obj.entries{1};

            for entryIndex = 2:length(obj.entries)
                entry = obj.entries{entryIndex};
                residueRows = obj.residues{1}((obj.residues{1}{:, 1} > 0) & (obj.residues{1}{:, entryIndex} > 0), :);

                entry.alignStructures(refEntry, residueRows{:, 1}, residueRows{:, entryIndex});
            end
        end

        function savePdbFiles(obj, targetPath)
            for i = 1:length(obj.entries)
                entry = obj.entries{i};

                if entry.name ~= ""
                    name = entry.name;
                else
                    name = "Entry " + i;
                end

                writepdb(fullfile(targetPath, name + ".pdb"), entry.pdb);
            end
        end

        function output = calcRmsd(obj, feature, refEntry)
            [~, refCrd] = refEntry.getAtoms('Backbone', 1, 'Chain', 1, 'Residues', feature(refEntry.databaseIndex, :));

            output = cell(1, length(obj.entries));

            for entryIndex = 1:length(obj.entries)
                entry = obj.entries{entryIndex};
                atomIndices = to3(entry.getAtoms('Backbone', 1, 'Chain', 1, 'Residues', feature(entryIndex, :)));

                if isempty(entry.simulation)
                    output{entryIndex} = calcrmsd(refCrd, entry.crd(atomIndices));
                else
%                     data = zeros(entry.simulation.runCount, size(entry.simulation.traj{1}, 1));
                    data = cell(entry.simulation.runCount,1);
                    for runIndex = 1:entry.simulation.runCount
                        traj = entry.simulation.traj{runIndex}(:, atomIndices);
                        data{runIndex} = calcrmsd(refCrd, traj);
%                         data(runIndex, :) = calcrmsd(refCrd, traj);
                    end

                    output{entryIndex} = data;
                end
            end
        end

        function output = calcDistance(obj, feature1, feature2)
            feature = [feature1 feature2];
            output = cell(1, length(obj.entries));

            for entryIndex = 1:length(obj.entries)
                entry = obj.entries{entryIndex};
                atomIndices = entry.getAtoms('Chain', 1, 'Name', 'CA', 'Residues', feature(entryIndex, :));
                
                if feature(entryIndex,1) == 0 || feature(entryIndex,2) == 0
                    output{entryIndex} = nan;
                else

                    if isempty(entry.simulation)
                        output{entryIndex} = calcbond(entry.crd, atomIndices');
                    else
    %                     data = zeros(entry.simulation.runCount, size(entry.simulation.traj{1}, 1));
                        data = cell(entry.simulation.runCount,1);
    
                        for runIndex = 1:entry.simulation.runCount
    %                         data(runIndex, :) = calcbond(entry.simulation.traj{runIndex}, atomIndices');
                            data{runIndex} = calcbond(entry.simulation.traj{runIndex}, atomIndices');
                        end
    
                        output{entryIndex} = data;
                    end
                end
            end
        end

        function output = calcDihedral(obj, feature, chainIndex)
            
            output =[];

            
            for entryIndex = 1:length(obj.entries)
                entry = obj.entries{entryIndex};
                % Take care of chains outside of calcdihedral function to avoid doubly numbered residues
                [indices, crd] =  entry.getAtoms('Chain',chainIndex);
                pdbChain = substruct(entry.pdb, indices);
                [dihedrals, ~, ~, ~] = calcalldihedralsfromtrajs(pdbChain, crd, feature(entryIndex, :), 1, 'all');
                output = [output dihedrals];
            end
            output = output';
        end
    end
end


function output = formatBw(bw)
   output = bw;
   fragments = bw.split(".");
   loop = fragments(1);

   if strlength(loop) == 2
      loop_str = char(loop);
      start_loop = str2double(loop_str(1));

      if mod(start_loop, 2) == 0
          output = string(['ECL' num2str(start_loop / 2)]);
      else
          output = string(['ICL' num2str((start_loop + 1) / 2)]);
      end
   end
end
