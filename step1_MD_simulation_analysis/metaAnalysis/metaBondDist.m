%% MI per residue
metadir = 'D:\Telework_library\dopamine_phase_3';
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
name = 'DD2R_ligand_binding_distances';
pdbName = 'prot.pdb';

% Basic: DA + BRC + inactive
% foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT', '2-d2_dop_T174M-C220L', '4-d2_bromo_T174M-C220L', ...
%      '10-d2_dop_L214M','11-d2_brc_L214M','12-d2_dop_L92G', '13-d2_brc_L92G', ...
%       '14-d2_dop_F217M', '15-d2_brc_F217M',  '16-d2_dop_I125N', '17-d2_brc_I125N','18a-d2_dop_F217I', '18b-d2_brc_F217I', ...
%       '19a-d2_dop_L92H','19b-d2_brc_L92H','21-d2_ris_6cm4'};
% 
% thisSysLabel = {'DA-WT-Gi','BRC-WT-Gi','DA-T5.54M-C6.47L-Gi','BRC-T5.54M-C6.47L-Gi', ...
%     'DA-L6.41M-Gi','BRC-L6.41M-Gi','DA-L3.41G-Gi','BRC-L3.41G-Gi', ...
%     'DA-F6.44M-Gi','BRC-F6.44M-Gi','DA-I4.46N-Gi','BRC-I4.46N-Gi','DA-F6.44I-Gi','BRC-F6.44I-Gi', ...
%     'DA-L3.41H-Gi','BRC-L3.41H-Gi','RIS-inactive'};

% DA, BRC, AP comparison
% foldersToStudy = {'1-d2_dop_WT','3-d2_bromo_WT','61-D2_apom_GiH5_WT',...
%       '14-d2_dop_F217M', '15-d2_brc_F217M', '63-D2_apom_GiH5_F217M' ...
%       '16-d2_dop_I125N', '17-d2_brc_I125N', '62-D2_apom_GiH5_I125N'};
% 
% thisSysLabel = {'DA-WT-Gi','BRC-WT-Gi','APM-WT-Gi', ...
% 'DA-F6.44M-Gi','BRC-F6.44M-Gi', 'APM-F6.44M-Gi'...
% 'DA-I4.46N-Gi','BRC-I4.46N-Gi', 'APM-I4.46N-Gi'};

% DA + SNPs
% foldersToStudy = {'1-d2_dop_WT', '2-d2_dop_T174M-C220L','10-d2_dop_L214M', ...
%     '12-d2_dop_L92G','8-d2_dop_V130F','14-d2_dop_F217M', '16-d2_dop_I125N', '18a-d2_dop_F217I', '19a-d2_dop_L92H',  ...
%      '70b-D2_DA_V444I_VC_rescale','71b-D2_DA_I181F_BW561_VC_rescale'};
% 
% thisSysLabel = {'DA-WT-Gi','DA-T5.54M-C6.47L-Gi','DA-L6.41M-Gi','DA-L3.41G-Gi','DA-V4.51F-Gi', 'DA-F6.44M-Gi', ...
%     'DA-I4.46N-Gi','DA-F6.44I-Gi', 'DA-L3.41H-Gi','DA-V4.44I','DA-I5.61F'};

md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc_Culled_data';

% Generic GPCR numbering from https://gpcrdb.org/residue/residuetable
gpcrdbRefName = 'DD2R';

% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
chains = 'A-B';

% Reference PDB code that this simulation is based on
pdbCode = '6VMS';
pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
pdbCodeInactive = '6CM4';
pdbInactiveChains = 'A--';

isGPCR = true;
frames2skip = 500; % 50 ns
useDatabaseAlignment = true; % Attempts to align structures/simulations using database alignment
attempt2AlignTrajs = false;
paperColors;
%%

%% Modernized meta PCA, fit for AlloDy 1.0:

metaLoadAlign;
colormapName = 'parula';
colors = colormap(colormapName);

% This is for receptor chain: ONLY CA PCA FOR NOW
if attempt2AlignTrajs
    simResList =  database.residues{1}{:,1:length(foldersToStudy)}; % list from simulation set
    allAligned = sum(simResList==0,2)==0;
    alignmentVector = simResList(allAligned,:); % Attempt to line different systems
    % alignmentVector works if:
    % 1- residue numbers start from 1
    % 2- residue numbers of same chain are sorted in increasing order (no funny 
    % PDB files where residue 500 is before residue 1 in the same chain)
end

tempProtTraj = [];
tempLigTraj = [];
CSys = []; % Contains labels as: [ system run frameNdx ]
nFramesEff = cell(length(foldersToStudy),1); % I will use this to extract
clusterNdxPCA = cell(length(foldersToStudy),1);
% the pdb of the center

% Read and concatenate trajs
for thisSys = 1:length(foldersToStudy)

    myDir = fullfile(metadir,foldersToStudy{thisSys});
    entryHere  = database.entries{thisSys};
    chainHere =  entryHere.chains{Chains.receptor};
    simHere = entryHere.addSimulation(fullfile(myDir, "run*"),'align2chain',chains(Chains.receptor));
    
    clusterNdxPCA{thisSys} = load(fullfile(myDir,'md2pathdev','workspace.mat'),"-mat","indexOfCluster_pca");
%     % Protein stuff
%     protAtomIndices = chainHere.getAtoms(); % Grab CA atoms
%     traj = simHere.concatRuns('Atoms', protAtomIndices, 'StartFrame', frames2skip + 1);
%     if attempt2AlignTrajs
%         tempProtTraj = [tempProtTraj ;traj(:,to3(alignmentVector(:,thisSys)))]; % Matlab will complain about this :D :D 
%     else
%         tempProtTraj = [tempProtTraj ;traj];
%     end
% 
%     % Ligand stuff
%     ligAtomIndices = entryHere.chains{Chains.ligand}.getLigandAtoms; % Assumes same ligand among different simulations
%     trajLig = simHere.concatRuns('Atoms', ligAtomIndices, 'StartFrame', frames2skip + 1);
%     tempLigTraj = [tempLigTraj ;trajLig];
%     % Grab the index Csys
%     numRuns = simHere.runCount;
%     for runi = 1:numRuns
%         framesNdx = ((frames2skip + 1):size(simHere.traj{runi},1))';
%         nFramesEff{thisSys}(runi) = size(simHere.traj{runi}((frames2skip + 1):end,to3(protAtomIndices)),1);
%         Ctemp = [thisSys*ones(nFramesEff{thisSys}(runi),1)  ...
%            runi*ones(nFramesEff{thisSys}(runi),1) framesNdx];
%         CSys=[CSys ;Ctemp]; % Used for coloring and for labeling: [ system run frameNdx]
%     end
end


%%
% Distances:
% TM3 - TM7
% TM3 - TM6
% TM3 - TM5
% TM5 - TM7
res3BW = "3.32";
res5BW = "5.46";
res6BW = "6.55";
res7BW = "7.43";
res3BWalt = "3.40";
res6BWalt = "6.44";

res3 = database.findResidue(res3BW);
res5 = database.findResidue(res5BW);
res6 = database.findResidue(res6BW); 
res7 = database.findResidue(res7BW); 
res3alt = database.findResidue(res3BWalt);
res6alt =  database.findResidue(res6BWalt);

params.dis36.data = database.calcDistance(res3,res6);
params.dis36.label = sprintf("%s-%s Ca-Ca (Å)",res3BW,res6BW);

params.dis57.data = database.calcDistance(res5,res7);
params.dis57.label = sprintf("%s-%s Ca-Ca (Å)",res5BW,res7BW);

params.dis56.data = database.calcDistance(res5,res6);
params.dis56.label = sprintf("%s-%s Ca-Ca (Å)",res5BW,res6BW);

params.dis67.data = database.calcDistance(res6,res7);
params.dis67.label = sprintf("%s-%s Ca-Ca (Å)",res6BW,res7BW);

params.dis36alt.data = database.calcDistance(res3alt,res6alt);
params.dis36alt.label = sprintf("%s-%s Ca-Ca (Å)",res3BWalt,res6BWalt);

params.dis56alt.data = database.calcDistance(res5,res6alt);
params.dis56alt.label = sprintf("%s-%s Ca-Ca (Å)",res5BW,res6BWalt);
%% Do the plotting
systemsToPlot =[1 3];

figure; 
errStrideRMSD = 50;
fieldsHere = fields(params);
% factor to transform from frames to nanoseconds:
timeFactor = 10;
Ang = char(197);

tiledlayout('flow')
for j = 1:length(fieldsHere) % TM36 and TM37 distances
nexttile
colorCount = 1;
    for i =systemsToPlot
       
        RMSD_mean= mean([params.(fieldsHere{j}).data{i}{:}],2);
        RMSD_std= std([params.(fieldsHere{j}).data{i}{:}],0,2);
        runCount = length(params.(fieldsHere{j}).data{i});

        legendEntry = sprintf('%s:%0.2f%s [%0.2f]',database.entries{i}.name, mean([params.(fieldsHere{j}).data{i}{:}],'all'),Ang,std([params.(fieldsHere{j}).data{i}{:}],0,'all')/sqrt(runCount));
        % WARNING: using undocumented feature: 4th element in RGB triplet
        % is Alpha
        p = plot((1:length(RMSD_mean))/timeFactor,RMSD_mean,'LineWidth',1,'Color',[paperColorPalette.RGB(2*colorCount-1,:) 0.25],'DisplayName',legendEntry);

%         scatter((1:length(RMSD_mean))/timeFactor,RMSD_mean,5,'filled','MarkerEdgeAlpha',0.25,'MarkerFaceAlpha',0.25, ...
%             'MarkerEdgeColor',paperColorPalette.Hex(2*colorCount-1),'MarkerFaceColor',paperColorPalette.Hex(2*colorCount-1),'DisplayName',database.entries{i}.name);

        hold on
        errorbar((1:errStrideRMSD:length(RMSD_mean))/timeFactor,RMSD_mean(1:errStrideRMSD:length(RMSD_mean)), ...
            RMSD_std(1:errStrideRMSD:length(RMSD_mean))/sqrt(runCount),'o','MarkerSize',5,'Color',paperColorPalette.Hex(2*colorCount-1),'HandleVisibility','off')
        colorCount = colorCount + 1;
    end
%     yline(params.(fieldsHere{j}).data{end-1}, '--','Active reference','LineWidth',1.5,'Color',[0 0.4470 0.7410],'DisplayName',database.entries{end-1}.name);
%     yline(params.(fieldsHere{j}).data{end}, '--','Inactive reference','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980],'DisplayName',database.entries{end}.name);
    legend('-DynamicLegend');
    legend boxoff
    xlabel('Time (ns)'); ylabel(params.(fieldsHere{j}).label)
    title(params.(fieldsHere{j}).label)
    set(gca,'FontSize',16)
    
end

%% Calculate simple statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
systemsToPlot = 1:length(thisSysLabel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pm = char(177); % plus or minus char
stats = zeros(length(systemsToPlot),length(fieldsHere)*3);
% statsCell = cell(length(fieldsHere)*2,length(systemsToPlot));
 statsCell = cell(length(fieldsHere),length(systemsToPlot));
for j = 1:length(fieldsHere) % TM36 and TM37 distances
colorCount = 1;
    for i =systemsToPlot
        if contains(thisSysLabel{i},'DA')
            refsysID = find(contains(thisSysLabel,'DA-WT'));
            refsysID = refsysID(1);
        elseif contains(thisSysLabel{i},'BRC')
            refsysID = find(contains(thisSysLabel,'BRC-WT'));
        end
        [M,SEM,p] =  getStats(params,j,i,frames2skip,"refsysID",refsysID);
        stats(colorCount,(3*j-2):3*j) = [M SEM p];
%         statsCell{2*j-1,colorCount} = sprintf('%.2f %s %.2f',M,pm,SEM);
%         statsCell{2*j,colorCount} = p;
        if i == refsysID
            statsCell{j,colorCount} = sprintf('%.2f %s %.2f',M,pm,SEM);
        else
            statsCell{j,colorCount} = sprintf('%.2f %s %.2f [%.4f]',M,pm,SEM,p);
        end
        colorCount = colorCount + 1;
    end
end


varNames =[];
for j = 1:length(fieldsHere) % TM36 and TM37 distances
    name1 = [fieldsHere{j} '_mean'];
    name2 = [fieldsHere{j} '_SEM'];
    name3 = [fieldsHere{j} '_pVal'];
    varNames = [varNames string(name1) string(name2) string(name3)];
end

% Plot data in a fun way:
figure('Position',1e3*[0.1210    0.1138    1.3752    0.6482])
tiledlayout('flow')
colorCount = 1;
for j = 1:length(fieldsHere) % TM36 and TM37 distances
    nexttile
    errorbar(1:length(systemsToPlot),stats(:,3*j-2), ...
        stats(:,3*j-1),'o','MarkerSize',10,'Color',paperColorPalette.Hex(1), ...
        'MarkerEdgeColor',paperColorPalette.Hex(1), 'MarkerFaceColor',paperColorPalette.Hex(2),'LineWidth',1.5)
    hold on
    colorCount = colorCount + 1;
    xticks([1:length(systemsToPlot)])
    ylabel(params.(fieldsHere{j}).label)
    set(gca,'FontSize',16)
    xticklabels(thisSysLabel(systemsToPlot))
    ax = gca;
    ax.XAxis.FontSize=12;
    set(gca,'linewidth',1.5)
%     formatplot2
end



% % Write table (old table format)
% tabHere = [table(thisSysLabel(systemsToPlot)','VariableNames',"SysName" ) array2table(stats,'VariableNames',varNames)];
% tabFileName = fullfile(metadir,'md2pathMeta',name,sprintf('%s_distance_table_stats_paper_pval.xlsx',gpcrdbRefName));
% writetable(tabHere,tabFileName)

% Write table with cell

 tabHere = [ table(fieldsHere ,'VariableNames',"FieldName" ) cell2table(statsCell,'VariableNames',thisSysLabel(systemsToPlot))];
tabFileName = fullfile(metadir,'md2pathMeta',name,sprintf('%s_distance_table_stats_paper_pval_altFormat.xlsx',gpcrdbRefName));
% writetable(tabHere,tabFileName)



%% Plot as one system + WT all distances on one plot: no PCA version

xaxisLabel = [];
for j = 1:length(fieldsHere)
    xaxisLabel = [xaxisLabel extractBefore(params.(fieldsHere{j}).label,10)];
end

% plot
for i = systemsToPlot  
    
    colorCount = 1;
    if contains(thisSysLabel{i},'DA')
        refsysID = find(contains(thisSysLabel,'DA-WT'));
        refsysID = refsysID(1);
    elseif contains(thisSysLabel{i},'BRC')
        refsysID = find(contains(thisSysLabel,'BRC-WT'));
    end
    if i == refsysID % Skip WT
        continue
    end
    figure('Position',[683.4000  237.8000  664.0000  436.8000])
    errorbar(1:length(stats(:,1:3:end-2)),stats(refsysID,1:3:end-2),stats(refsysID,2:3:end-1), ...
        'o','MarkerSize',10,'Color',paperColorPalette.Hex(2*colorCount-1), ...
        'MarkerEdgeColor',paperColorPalette.Hex(2*colorCount-1), 'MarkerFaceColor',paperColorPalette.Hex(2*colorCount),'LineWidth',1.5, ...
        'DisplayName',thisSysLabel{refsysID})
    hold on
    colorCount = colorCount + 1;
    errorbar(1:length(stats(:,1:3:end-2)),stats(i,1:3:end-2),stats(i,2:3:end-1), ...
        'o','MarkerSize',10,'Color',paperColorPalette.Hex(2*colorCount-1), ...
        'MarkerEdgeColor',paperColorPalette.Hex(2*colorCount-1), 'MarkerFaceColor',paperColorPalette.Hex(2*colorCount),'LineWidth',1.5, ...
        'DisplayName',thisSysLabel{i})

    % Add pvalue as text above:
    pText = [];
    pValHere = stats(i,3:3:end);
    for pInd = 1:length(pValHere)
        if pValHere(pInd) < 0.0001
            pText = [pText "p<0.0001"];
        else
            pText = [pText sprintf("p=%.3f",pValHere(pInd))];
        end
    end
    text(1:length(pValHere), ...
        max([stats(refsysID,1:3:end-2)+stats(refsysID,2:3:end-1) ; stats(i,1:3:end-2)+stats(i,2:3:end-1)])+1, ...
        pText,'HorizontalAlignment','center','FontSize',16)


    xticks(1:length(stats(:,1:3:end-2)))
    ylabel("C\alpha-C\alpha distance (Å)")
    legend('-DynamicLegend','Location','best');
    legend boxoff
    set(gca,'FontSize',20)
    xticklabels(xaxisLabel)
    xlim([0 length(stats(:,1:3:end-2))+1])
%     ax = gca;
%     ax.XAxis.FontSize=12;
    set(gca,'linewidth',1.5)
end













%% Experiment with PCA clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
systemsToPlot = 1:length(thisSysLabel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clusterNdx = [1];


pm = char(177); % plus or minus char
stats = cell(length(clusterNdx),1);
for k = 1:length(clusterNdx)
    stats{k} = zeros(length(systemsToPlot),length(fieldsHere)*3);
end

% statsCell = cell(length(fieldsHere)*2,length(systemsToPlot));
 statsCell = cell(length(fieldsHere),length(systemsToPlot));
for j = 1:length(fieldsHere) % TM36 and TM37 distances
colorCount = 1;
    for i =systemsToPlot
        if contains(thisSysLabel{i},'DA')
            refsysID = find(contains(thisSysLabel,'DA-WT'));
            refsysID = refsysID(1);
        elseif contains(thisSysLabel{i},'BRC')
            refsysID = find(contains(thisSysLabel,'BRC-WT'));
        end
        if i ~= refsysID
            for k = 1:length(clusterNdx)
                if max(clusterNdxPCA{i}.indexOfCluster_pca) < k
                    stats{k}(colorCount,(3*j-2):3*j) = [nan nan nan];
                    continue
                end
                [M,SEM,p] =  getStats(params,j,i,frames2skip,"refsysID",refsysID,"clusterFrameNdx",clusterNdxPCA{i}.indexOfCluster_pca==clusterNdx(k));
                stats{k}(colorCount,(3*j-2):3*j) = [M SEM p];
            end
            
        else

            [M,SEM,p] =  getStats(params,j,i,frames2skip,"refsysID",refsysID);
            for k = 1:length(clusterNdx)
                stats{k}(colorCount,(3*j-2):3*j) = [M SEM p];
            end
        end
%          stats(colorCount,(3*j-2):3*j) = [M SEM p];
%         statsCell{2*j-1,colorCount} = sprintf('%.2f %s %.2f',M,pm,SEM);
%         statsCell{2*j,colorCount} = p;
        if i == refsysID
            statsCell{j,colorCount} = sprintf('%.2f %s %.2f',M,pm,SEM);
        else
            statsCell{j,colorCount} = sprintf('%.2f %s %.2f [%.4f]',M,pm,SEM,p);
        end
        colorCount = colorCount + 1;
    end
end

% Plot data in a fun way:
figure('Position',1e3*[0.1210    0.1138    1.3752    0.6482])
tiledlayout('flow')
% sgtitle(sprintf("Variant PCA cluster %d",clusterNdx))
for j = 1:length(fieldsHere) % TM36 and TM37 distances
    colorCount = 1;
    nexttile
    for k = 1:length(clusterNdx)
    errorbar(1:length(systemsToPlot),stats{k}(:,3*j-2), ...
        stats{k}(:,3*j-1),'o','MarkerSize',10,'Color',paperColorPalette.Hex(2*colorCount-1), ...
        'MarkerEdgeColor',paperColorPalette.Hex(2*colorCount-1), 'MarkerFaceColor',paperColorPalette.Hex(2*colorCount),'LineWidth',1.5, ...
        'DisplayName',['PCA C' num2str(clusterNdx(k))])
    hold on
    colorCount = colorCount + 1;
    end
    if j==1
        legend('-DynamicLegend','Location','best');
        legend boxoff
    end
    xticks([1:length(systemsToPlot)])
    ylabel(params.(fieldsHere{j}).label)
    set(gca,'FontSize',16)
    xticklabels(thisSysLabel(systemsToPlot))
    ax = gca;
    ax.XAxis.FontSize=12;
    set(gca,'linewidth',1.5)
%     formatplot2
end

%% Plot as one system + WT all distances on one plot:
% prep
k = 1; % Plot main PCA cluster
xaxisLabel = [];
for j = 1:length(fieldsHere)
    xaxisLabel = [xaxisLabel extractBefore(params.(fieldsHere{j}).label,10)];
end

% plot
for i = systemsToPlot  
    
    colorCount = 1;
    if contains(thisSysLabel{i},'DA')
        refsysID = find(contains(thisSysLabel,'DA-WT'));
        refsysID = refsysID(1);
    elseif contains(thisSysLabel{i},'BRC')
        refsysID = find(contains(thisSysLabel,'BRC-WT'));
    end
    if i == refsysID % Skip WT
        continue
    end
    figure('Position',[683.4000  237.8000  664.0000  436.8000])
    errorbar(1:length(stats{k}(:,1:3:end-2)),stats{k}(refsysID,1:3:end-2),stats{k}(refsysID,2:3:end-1), ...
        'o','MarkerSize',10,'Color',paperColorPalette.Hex(2*colorCount-1), ...
        'MarkerEdgeColor',paperColorPalette.Hex(2*colorCount-1), 'MarkerFaceColor',paperColorPalette.Hex(2*colorCount),'LineWidth',1.5, ...
        'DisplayName',thisSysLabel{refsysID})
    hold on
    colorCount = colorCount + 1;
    errorbar(1:length(stats{k}(:,1:3:end-2)),stats{k}(i,1:3:end-2),stats{k}(i,2:3:end-1), ...
        'o','MarkerSize',10,'Color',paperColorPalette.Hex(2*colorCount-1), ...
        'MarkerEdgeColor',paperColorPalette.Hex(2*colorCount-1), 'MarkerFaceColor',paperColorPalette.Hex(2*colorCount),'LineWidth',1.5, ...
        'DisplayName',thisSysLabel{i})

    % Add pvalue as text above:
    pText = [];
    pValHere = stats{k}(i,3:3:end);
    for pInd = 1:length(pValHere)
        if pValHere(pInd) < 0.0001
            pText = [pText "p<0.0001"];
        else
            pText = [pText sprintf("p=%.3f",pValHere(pInd))];
        end
    end
    text(1:length(pValHere), ...
        max([stats{k}(refsysID,1:3:end-2)+stats{k}(refsysID,2:3:end-1) ; stats{k}(i,1:3:end-2)+stats{k}(i,2:3:end-1)])+1, ...
        pText,'HorizontalAlignment','center','FontSize',16)


    xticks(1:length(stats{k}(:,1:3:end-2)))
    ylabel("C\alpha-C\alpha distance (Å)")
    legend('-DynamicLegend','Location','best');
    legend boxoff
    set(gca,'FontSize',20)
    xticklabels(xaxisLabel)
    xlim([0 length(stats{k}(:,1:3:end-2))+1])
%     ax = gca;
%     ax.XAxis.FontSize=12;
    set(gca,'linewidth',1.5)
end



%%
j=3
drawSwarm(params,j,[1 2 3 4],refsysID,frames2skip,'plotPval',1)

title(params.(fieldsHere{j}).label)
    set(gca,'FontSize',20)

%% Support functions

function [M,SEM,p,d] = getStats(params,field,sysID,frames2skip,options)
    arguments
    params
    field
    sysID
    frames2skip = 0
    options.refsysID = []
    options.calcPvalMeanPerRun = false
    options.clusterFrameNdx = [] % Logical array same size as data
    end
    fieldsHere = fields(params);
    runCount = length(params.(fieldsHere{field}).data{sysID});
    data = [];
    dataMeanPerRun = zeros(runCount,1);
    for i = 1:runCount
        data = [data; params.(fieldsHere{field}).data{sysID}{i}(frames2skip+1:end)];
        dataMeanPerRun(i) = mean(params.(fieldsHere{field}).data{sysID}{i}(frames2skip+1:end));
    end
    if ~isempty(options.clusterFrameNdx) 
        data = data(options.clusterFrameNdx);
    end

    if ~isempty(options.refsysID) % Calc pval with reference
        runCountRef = length(params.(fieldsHere{field}).data{options.refsysID});
        dataRef = [];
        dataMeanPerRunRef = zeros(runCountRef,1);
        for i = 1:runCountRef
            dataRef = [dataRef; params.(fieldsHere{field}).data{options.refsysID}{i}(frames2skip+1:end)];
            dataMeanPerRunRef(i) = mean(params.(fieldsHere{field}).data{options.refsysID}{i}(frames2skip+1:end));
        end
        if options.calcPvalMeanPerRun
            [~,p] = ttest2(dataMeanPerRunRef,dataMeanPerRun,'Vartype','unequal');
        else
            [~,p] = ttest2(dataRef,data,'Vartype','unequal');
        end
    end

    M = mean(data);
    SEM = std(data)/sqrt(runCount);
    
    if nargout>3
        d = cohensD(dataRef, data);
    end
end

function d = cohensD(x, y)
    % cohensD calculates the Cohen's d effect size between two vectors
    % Inputs:
    %   x - first vector of data
    %   y - second vector of data
    % Output:
    %   d - Cohen's d value
    
    % Check if the inputs are non-empty vectors
    if isempty(x) || isempty(y)
        error('Input vectors must not be empty');
    end
    
    % Calculate the means of both vectors
    mean_x = mean(x);
    mean_y = mean(y);
    
    % Calculate the standard deviations of both vectors
    std_x = std(x, 0);  % using N-1 (default) for the sample standard deviation
    std_y = std(y, 0);  % using N-1 (default) for the sample standard deviation
    
    % Calculate the pooled standard deviation
    n_x = length(x);
    n_y = length(y);
    pooled_std = sqrt(((n_x - 1) * std_x^2 + (n_y - 1) * std_y^2) / (n_x + n_y - 2));
    
    % Calculate Cohen's d
    d = (mean_x - mean_y) / pooled_std;
end


function  drawSwarm(params,field,sysID,refsysID,frames2skip,options)
    arguments
    params
    field
    sysID
    refsysID = 1
    frames2skip = 0
    options.calcPvalMeanPerRun = false
    options.plotPval = false
    options.clusterFrameNdx = [] % Logical array same size as data
    end
    if iscolumn(sysID)
        sysID = sysID';
    end
    figure
    fieldsHere = fields(params);
    
    % Get info on reference:
    runCountRef = length(params.(fieldsHere{field}).data{refsysID});
    dataRef = [];
    dataMeanPerRunRef = zeros(runCountRef,1);
    for i = 1:runCountRef
        dataRef = [dataRef; params.(fieldsHere{field}).data{refsysID}{i}(frames2skip+1:end)];
        dataMeanPerRunRef(i) = mean(params.(fieldsHere{field}).data{refsysID}{i}(frames2skip+1:end));
    end
    swarmchart(ones(size(dataRef)),dataRef,'.','MarkerEdgeAlpha',0.5)
    hold on
    counter = 1;
    p = zeros(size(sysID));
    for sysHere = sysID
        if sysHere == refsysID % Skip WT
            continue
        end
        runCount = length(params.(fieldsHere{field}).data{sysHere});
        data = [];
        dataMeanPerRun = zeros(runCount,1);
        for i = 1:runCount
            data = [data; params.(fieldsHere{field}).data{sysHere}{i}(frames2skip+1:end)];
            dataMeanPerRun(i) = mean(params.(fieldsHere{field}).data{sysHere}{i}(frames2skip+1:end));
        end
        if ~isempty(options.clusterFrameNdx) 
            data = data(options.clusterFrameNdx);
        end
    

        if options.calcPvalMeanPerRun
            [~,p(counter)] = ttest2(dataMeanPerRunRef,dataMeanPerRun,'Vartype','unequal');
        else
            [~,p(counter)] = ttest2(dataRef,data,'Vartype','unequal');
        end
        pValHere = p(counter);
        counter = counter + 1; % 1 is reserved for reference
        swarmchart(counter*ones(size(data)),data,'.','MarkerEdgeAlpha',0.5)
        hold on
        if options.plotPval
            pText = [];
            if pValHere < 0.0001
                pText = "p<0.0001";
            else
                pText = sprintf("p=%.3f",pValHere);
            end
            text(counter, ...
                max(data)+1, ...
                pText,'HorizontalAlignment','center','FontSize',16)
        end
    end

   
end