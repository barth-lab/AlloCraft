classdef Chain < handle
    properties
        entry
        index
        name

        atomIndices
        resIds

        helices
        secStruct
    end

    properties (Dependent)
        isSmall
        seq
    end

    methods
        function obj = Chain(entry, name) % Called in Entry.m
            obj.index = find(char(name) == entry.chainNames);
            obj.entry = entry;
            obj.name = name;

            obj.atomIndices = selectname(entry.pdb.chainid, name);
            obj.resIds = unique(entry.pdb.resseq(entry.pdb.chainid == char(name)),'stable');
        end

        function [pdb, crd] = export(obj)
            pdb = substruct(obj.entry.pdb, obj.atomIndices);
            crd = obj.entry.crd(to3(obj.atomIndices));
        end

        function value = get.seq(obj) 
            % Getting sequence depends on chainIndices of index, and not
            % the index itself
%             value = obj.entry.seq{obj.index};
            value = obj.entry.seq{obj.entry.chainIndices(obj.index)};
        end

        function value = get.isSmall(obj)
            value = (length(obj.resIds) <= 1);
        end

        function computeSecStruct(obj)
            outPath = tempname();
            system("vmd -dispdev text -e " + fullfile(pwd, "extract_ss.tcl") + " -args " + obj.entry.path + " " + outPath);
            data = importdata(outPath);

            if iscell(data)
                data = data{1};
            end

            data(isspace(data)) = [];

            obj.secStruct = data(obj.resIds);
            obj.helices = findHelices(obj.secStruct);
        end

        function output = getAtoms(obj)
            output = obj.entry.getAtoms('Chain', obj.index, 'Name', 'CA');
        end

        function output = getLigandAtoms(obj)
            if obj.isSmall
%                 output = obj.entry.getAtoms('Chain', obj.index, 'NoName', 'H*');
                indices = selectname(obj.entry.pdb.chainid, obj.name);
                indices = indices &  ~selectname(obj.entry.pdb.name, 'H*','LP*');
                output = find(indices);
            else
                output = obj.entry.getAtoms('Chain', obj.index, 'Name', 'CA');
            end
        end

        function output = formatName(obj)
            if obj.isSmall && (length(obj.resIds) == 1)
                output = strtrim(string(obj.entry.pdb.resname(find(obj.atomIndices, 1, 'first'), :)));
            else
                output = "Chain " + obj.name;
            end
        end

        function output = formatAtoms(obj, atomIndices)
            arguments
                obj
                atomIndices = find(obj.atomIndices)
            end

            names = obj.entry.pdb.name(atomIndices, :);
            names = strtrim(string(names));

            output = names;

            % output = strings(length(atomIndices), 1);

            % for index = 1:length(atomIndices)
            %     atomIndex = atomIndices(index);
            %     name = char(names(index));
            %     output(index) = sprintf("%c%s", name);
            % end
        end

        function output = formatResidues(obj, resIds, options)
            arguments
                obj
                resIds
                options.PrimaryEntry = obj.entry
                options.SecondaryEntry
                options.BWonly = 0
            end

            if ~isfield(options, 'SecondaryEntry')
                if length(obj.entry.database.entries) > 1
                    options.SecondaryEntry = obj.entry.database.entries{2};
                else
                    options.SecondaryEntry = options.PrimaryEntry;
                end
            end

            rowFilter = ismember(obj.entry.database.residues{obj.index}{:, obj.entry.databaseIndex}, resIds);
            rows = obj.entry.database.residues{obj.index}(rowFilter, :);

            output = strings(height(rows), 1);

            for rowIndex = 1:height(rows)
                row = rows(rowIndex, :);
                primaryId = row{1, options.PrimaryEntry.databaseIndex};
                primaryName = row.Name(options.PrimaryEntry.databaseIndex);
                secondaryId = row{1, options.SecondaryEntry.databaseIndex};

                if options.PrimaryEntry == options.SecondaryEntry
                    secondaryId = 0;
                end

                if options.BWonly

                    output(rowIndex) = sprintf("%c%u (%s)", primaryName, secondaryId, row.Label);
                    if isempty(output(rowIndex))
                        output(rowIndex) = "";
                    end
                else
                    if (primaryId > 0) && (secondaryId > 0) && (row.Label ~= "")
                        output(rowIndex) = sprintf("%c%u (%u, %s)", primaryName, primaryId, secondaryId, row.Label);
                    elseif (primaryId > 0) && (secondaryId > 0)
                        output(rowIndex) = sprintf("%c%u (%u)", primaryName, primaryId, secondaryId);
                    elseif (primaryId > 0) && (row.Label ~= "")
                        output(rowIndex) = sprintf("%c%u (%s)", primaryName, primaryId, row.Label);
                    elseif (primaryId > 0)
                        output(rowIndex) = sprintf("%c%u", primaryName, primaryId);
                    else
                        % error
                        output(rowIndex) = "";
                    end
                end
            end
        end

        function output = concatResIds(obj, otherChain)
            if obj.resIds(end) < otherChain.resIds(1)
                output = [obj.resIds; otherChain.resIds];
            else
                output = [otherChain.resIds; obj.resIds];
            end
        end
    end
end
