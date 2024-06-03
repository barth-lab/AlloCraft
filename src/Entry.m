classdef Entry < handle
    properties
        crd
        name
        path
        pdb
        seq

        chains
        chainIndices
        chainNames

        database
        databaseIndex

        simulation
    end

    properties (Dependent)
        atomCount
    end

    methods
        function obj = Entry(path, chainNames, name)
            if ~exist('name', 'var')
                name = "";
            end

            [obj.pdb, obj.crd] = readpdb(path);
            pdbChainNames = unique(obj.pdb.chainid, 'stable');

            obj.name = name;
            obj.path = path;
            obj.seq = pdb2seq(obj.pdb);

            chainCount = length(chainNames);

            obj.chains = cell(length(obj.chainNames), 1);
            obj.chainIndices = arrayfun(@(chainName) findIndex(chainName == pdbChainNames), chainNames); %%
            obj.chainNames = chainNames;

            for chainItemIndex = 1:chainCount
                chainIndex = obj.chainIndices(chainItemIndex);
                chainName = chainNames(chainItemIndex);

                if (chainName ~= ' ') && (chainName ~= '-') && (chainIndex < 1)
                    error("Missing chain '%c' in entry '%s'", chainName, obj.name);
                end

                if chainIndex > 0 % Call the Chain class!
                    obj.chains{chainItemIndex} = Chain(obj, chainName);
                else
                    obj.chains{chainItemIndex} = [];
                end
            end
        end

        function count = get.atomCount(obj)
            count = length(obj.pdb.name);
        end

        function alignStructures(obj, refEntry, refResidues, objResidues)
            refAtomIndices = find(ismember(refEntry.pdb.resseq, refResidues) & selectname(refEntry.pdb.name, 'CA') & refEntry.pdb.chainid == refEntry.chainNames(1));
            otherAtomIndices = find(ismember(obj.pdb.resseq, objResidues) & selectname(obj.pdb.name, 'CA') & obj.pdb.chainid == obj.chainNames(1));

            [~, obj.crd] = superimpose_double_index(refEntry.crd, obj.crd, otherAtomIndices, refAtomIndices);
            obj.pdb.xyz = reshape(obj.crd, 3, [])';
        end

        function [indices, crd] = getAtoms(obj, options)
            arguments
                obj
                options.Backbone
                options.Chain
                options.Name
                options.NoName
                options.Residues
            end

            indices = ones(obj.atomCount, 1);

            if isfield(options, "Backbone")
                backbone = {"CA" "C" "N" "O"};
                indices = indices & selectname(obj.pdb.name, backbone{:});
            end

            if isfield(options, "Chain")
                indices = indices & selectname(obj.pdb.chainid, obj.chainNames(options.Chain));
            end

            if isfield(options, "Name")
                indices = indices & selectname(obj.pdb.name, options.Name);
            end

            if isfield(options, "NoName")
                indices = indices & ~selectname(obj.pdb.name, options.NoName);
            end

            if isfield(options, "Residues")
                indices = indices & selectid(obj.pdb.resseq, options.Residues);
            end

            indices = find(indices);
            crd = obj.crd(to3(indices));
        end

        function simulation = addSimulation(obj, path, options)
            arguments
                obj
                path
                options.align2chain ='all'
            end
            obj.simulation = Simulation.query(obj, path);
            if strcmp(options.align2chain, 'all') % Take all CAs for alignment
                obj.simulation.align();
            else % Take only specified chain
                obj.simulation.align('chain' ,options.align2chain);
            end

            simulation = obj.simulation;
        end
        % This function will not use the dual cutoff scheme as it's
        % designed to compute contacts from static reference structure
        % (rather than dynamic trajectories)
        function contacts = computeContacts(obj, chainA, chainB, options)
            arguments
                obj
                chainA
                chainB
                options.ContactCut = 0.4
                options.ResIds = chainB.resIds
                options.RCut = 5
            end

            [~, ~, ~, ~, contacts] = proteinContacts(obj.pdb, obj.crd, options.ResIds, char(chainB.name), char(chainA.name), options.RCut, options.ContactCut);
            contacts = contacts{1};
        end

        function output = hasChain(obj, chainIndex)
            output = (chainIndex <= length(obj.chains)) && ~isempty(obj.chains{chainIndex});
        end
    end
end


function index = findIndex(value)
    index = find(value);

    if isempty(index)
        index = 0;
    end
end
