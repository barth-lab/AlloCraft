classdef HyperDataBase < handle
    %hyperDataBase Collects hyperEntries in one compact place
    %   Detailed explanation goes here

    properties
        hyperEntries
        entryName
        entryDataSize
        entryDataType
        alignment
    end

    methods
        function obj = HyperDataBase()
            %UNTITLED12 Construct an instance of this class
            %   Detailed explanation goes here
            obj.hyperEntries = cell(1, 0);
            obj.entryName = cell(1, 0);
        end

        function read(obj, rawData, labels, name, options)
            arguments
                obj
                rawData
                labels
                name
                options.VarNames = []
            end
            hyperEntry = HyperEntry(rawData,labels,name,'VarNames',options.VarNames);
            obj.hyperEntries{end + 1} = hyperEntry;
            obj.entryName{end + 1} = name;
            obj.entryDataSize{end+1} = size(rawData);
            obj.entryDataType{end+1} = class(rawData);

            hyperEntry.hyperDataBase = obj;
            hyperEntry.hyperDataBaseNdx = length(obj.hyperEntries);
        end

        function align(obj,options)
            % align: aligns the data according to similar labels, a "0"
            % means no alignment is found. Assumes first entry is reference
            % unless "refEntry" option is used
            arguments
            obj
            options.refEntry = 1 
            end
            
            refLabels = string(obj.hyperEntries{options.refEntry}.labels);
            ndxLabels = zeros(length(refLabels), length(obj.hyperEntries));
            % Option 1: look for the same labels in reference entry and align

%             (1:end ~= options.refEntry)

            for j = 1:length(obj.hyperEntries)
                if j == options.refEntry
                    ndxLabels(:,j) = 1:length(refLabels);               
                else
                    jLabels = string(obj.hyperEntries{j}.labels);
                    for i = 1:length(refLabels)
                    % strcmp(string(thisSysLabel),string(tabExp.Label(i)))
                    ndxHere = find(strcmp(jLabels,refLabels(i)));
                    if ndxHere
                        ndxLabels(i,j) = ndxHere;
                    end
                    end
                end
            end
            obj.alignment = [table(refLabels) array2table(ndxLabels,'VariableNames',obj.entryName)];
            
            % Option 2: make a special column with "union" that everything
            % aligns to
        end

        function dbFetchLabels(obj,hyperEntryNdx)

        end

        function dataTable = dbFetchData(obj,str,options)
            % dbFetchData: fetches data from labels matching "str" from 
            % selected hyperEntries that have a valid alignment. Defaults 
            % to attempt fetching from all hyperEntires
            arguments
            obj
            str
            options.hyperEntryNdx = 1:length(obj.hyperEntries)
            options.refEntry = 1 
            end
            assert(ismember(options.refEntry,options.hyperEntryNdx), ...
                'Make sure refEntry is included in hyperEntryNdx!!!!')
            data = cell(1,length(options.hyperEntryNdx));
            
            % Fetch labels and data from reference index:
            refLabelNdx = obj.hyperEntries{options.refEntry}.fetchLabels(str);
            % Fetch labels for full alignments for the chosen entries:
            A=obj.alignment{:,options.hyperEntryNdx+1}; 
            refLabelAligned = prod(A,2)~=0;
            
            dataTable = table(obj.alignment{refLabelNdx&refLabelAligned,1},'VariableNames',"Label");
            for i = options.hyperEntryNdx
                data{i} = obj.hyperEntries{i}.orderData(obj.alignment{refLabelNdx&refLabelAligned,i+1});
                % Make sure that the order of the data is right using the
                % alignment!
                dataTable = [dataTable table(data{i}','VariableNames',obj.entryName(i))];
            end
        end
    end
end