classdef HyperEntry < handle
    %hyperData Class to store labeled data 
    %   rawData can be of any form where the columns are variables and rows
    %   are observations. Size of labels should be the same as number of
    %   columns of input data

    properties
        rawData
        labels
        name
        VarNames

        hyperDataBase
        hyperDataBaseNdx
    end

    methods
        function obj = HyperEntry(rawData,labels,name,options)
            arguments
                rawData
                labels
                name
                options.VarNames = []
            end
            %hyperData Construct an instance of this class
            %   Detailed explanation goes here
            assert(size(rawData,2)== length(labels),"Columns of input data and length of labels are not consistent!")
            obj.rawData = rawData;
            if iscolumn(labels)
                obj.labels = labels;
            else
                obj.labels = labels';
            end
            obj.name = name;
            obj.VarNames = options.VarNames;
        end

        function strNdx = fetchLabels(obj,str)
            % fetchLabel: searches contents of labels for input patterns
            % str. str could be a string or a string array
            strNdx = ones(length(obj.labels),1);

            for i =1:length(str)
                strNdx = strNdx & contains(obj.labels,str(i));
            end
        end
        function data = fetchData(obj,str,options)
            arguments
                obj
                str
                options.Ndx = true(length(obj.labels),1)
            end
            % fetchData: searches contents of rawData for input pattern in
            % Labels
            data = obj.rawData(:,obj.fetchLabels(str) & options.Ndx);
        end
        function data = orderData(obj,Ndx)
            % orderData: orders rows of data according to input ndx
            arguments
                obj
                Ndx
            end
            % fetchData: searches contents of rawData for input pattern in
            % Labels
            data = obj.rawData(:,Ndx);
        end

        function h = plotDataPairs(obj,NdxRef,NdxTest,options)
            % plotDataPairs: plots pairs of rows rawData with difference
            % and labels
            arguments
                obj
                NdxRef
                NdxTest
                options.SaveName = "%s"
                options.SavePath
            end
            paperColors; % Loads paperColorPalette
            
            h = figure('Position',1.0e3*[0.3994    0.0762    1.1056    0.7000]);
            tiledlayout('flow')
            nexttile
            xData = 1:size(obj.rawData,1);
            s1 =scatter(xData, obj.rawData(:,NdxTest),25, ...
                           'LineWidth',1.5, ...
                           'MarkerEdgeColor',paperColorPalette.Hex(1), ...
                           'MarkerFaceColor',paperColorPalette.Hex(2));
            hold on
            s2 =scatter(xData, obj.rawData(:,NdxRef),25, ...
                           'LineWidth',1.5, ...
                           'MarkerEdgeColor',paperColorPalette.Hex(3), ...
                           'MarkerFaceColor',paperColorPalette.Hex(4));
            row = dataTipTextRow('VarNames',obj.varNames);
            s1.DataTipTemplate.DataTipRows(end+1) = row;
            s2.DataTipTemplate.DataTipRows(end+1) = row;
            
            % row = dataTipTextRow('Label',tabExp.Label(ligNdx&effectorNdx));
            % s.DataTipTemplate.DataTipRows(end+1) = row
%             drawTMhelices(sum(MIres),helices, mainChain.resIds)
            xlabel('Residue')
            ylabel(obj.name)
            legend([mainEntry.name ', \Sigma MI(i,j) = ' num2str(round(sum(obj.rawData(:,NdxTest))))], ...
                [refEntry.name ', \Sigma MI(i,j) = ' num2str(round(sum(obj.rawData(:,NdxRef))))],'Location','bestoutside')
            
            title(sprintf("%s, %s, %s",obj.name,obj.labels{NdxTest},obj.labels{NdxRef}))
            set(gca,'FontSize',16)
            
            % Difference plot:
            % figure; 
            nexttile
            deltaData = obj.rawData(:,NdxTest)-obj.rawData(:,NdxRef);
            p = plot(deltaData,'-o');
            hold on
%             drawTMhelices(deltaData,helices, mainChain.resIds)
%             row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds));
            p.DataTipTemplate.DataTipRows(end+1) = row;
            grid on
            
            legend_entries{1} = ['\Delta ' obj.name];
            hold on
            legend_count = 2;
            
            receptorResIdsNdx = zeros(max(receptorResIds),1); % Helpful for indexing
            receptorResIdsNdx(receptorResIds)=1:length(receptorResIds); % Hoping this would do the trick
            
            % Add mutations (from reference structure)
            mutPos = find(database.residues{1}.Name(:,1) ~= database.residues{1}.Name(:,2));
            % Remove any difference in length where alignment is not present, those are
            % probably not mutations
            mutPos(isspace(database.residues{1}.Name(mutPos,2)) | isspace(database.residues{1}.Name(mutPos,1))) = [];
            mutRes = table2array(database.residues{1}(mutPos,1)); % Mutated residues
            mutations = cellstr([database.residues{1}.Name(mutPos,2) num2str(mutRes) database.residues{1}.Name(mutPos,1)  ]);
            
            if settings.includeLigPathwayCalc % do same for peptide ligand
                mutPosLig = find(database.residues{3}.Name(:,1) ~= database.residues{3}.Name(:,2));
                mutations = [mutations ; cellstr( [database.residues{3}.Name(mutPosLig,2) num2str(table2array(database.residues{3}(mutPosLig,1))) database.residues{3}.Name(mutPosLig,1)])  ];
                mutRes = [mutRes ; table2array(database.residues{3}(mutPosLig,1))];
            end
            
            % Add Annotations to the figure
            if exist(fullfile(md2pathdir,'BS_residues.txt'),'file')
                bsRes.main = importdata(fullfile(md2pathdir,'BS_residues.txt'));
                scatter(bsRes.main,deltaData(receptorResIdsNdx(bsRes.main)), 'LineWidth',1.5)
                legend_entries{legend_count} = 'Binding residues';
                legend_count = legend_count + 1;
            end
            if exist(fullfile(settings.refdir,'md2path','BS_residues.txt'),'file')
                bsRes.ref = importdata(fullfile(settings.refdir,'md2path','BS_residues.txt'));
            
                % Translate bsRef to main numbering before plotting:
                ndxMainRef = [database.residues{1}.output1 database.residues{1}.output2];
                bsResRefRenum = zeros(size(bsRes.ref));
                for resi = 1:length(bsRes.ref)
                    bsResRefRenum(resi) = ndxMainRef(ndxMainRef(:,2)==bsRes.ref(resi),1);
                end
                scatter(bsResRefRenum,deltaData(bsResRefRenum), 25, 'filled')
                legend_entries{legend_count} = 'Ref binding residues';
                legend_count = legend_count + 1;
            end
            
            if ~isempty(mutations)
                s3 = scatter(mutRes, deltaData(receptorResIdsNdx(mutRes)),80, 'LineWidth',1.5);
                text(mutRes,deltaData(receptorResIdsNdx(mutRes)) +0.05*max(deltaData),mutations)
                row = dataTipTextRow('Mut',mutations);
                s3.DataTipTemplate.DataTipRows(end+1) = row;
            
                legend_entries{legend_count} = 'Mutation sites';
                legend_count = legend_count + 1;
            end
            
            if ~isempty(settings.highlightRes)
                scatter(settings.highlightRes, deltaData(receptorResIdsNdx(settings.highlightRes)),60,'x', 'LineWidth',1.5);
                legend_entries{legend_count} = settings.highlightText;
                legend_count = legend_count + 1;
            end
            
            legend(legend_entries,'Location','bestoutside')
            xlabel('Residue')
            ylabel(sprintf("delta %s, (%s - %s)",obj.name,obj.labels{NdxTest},obj.labels{NdxRef}))
            title([obj.name ', ' obj.labels{NdxTest} ' - ' obj.labels{NdxRef} ])
            set(gca,'FontSize',16)
            
            if isfield(options, 'SavePath')
            figPath = fullfile(options.SavePath,sprintf("delta%s_%s_%s",obj.name,obj.labels{NdxTest},obj.labels{NdxRef}));
            savefig(figPath + ".fig");
            print2pdf(figPath+ ".pdf");
            saveas(gcf,figPath + ".svg",'svg')
            end
        end
    end
end