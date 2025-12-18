% A script to analyze md2path allosteric pathway output from several MD runs:

% This needs to be fixed for runs where the ligand is included 
% Goals: 
% 1- Compare top scoring residues between runs
% 2- Check positions that appear in top allosteric pipelines 
% 3- 

%% Prepare input files:

% Father directory for all the runs
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 
metadir = 'D:\Telework_library\dopamine_phase_3';
% foldersToStudy = {'1-d2_dop_WT','2-d2_dop_T174M-C220L','8-d2_dop_V130F','10-d2_dop_L214M','12-d2_dop_L92G','14-d2_dop_F217M'};

foldersToStudy = {'1-d2_dop_WT','2-d2_dop_T174M-C220L'};
% 
% foldersToStudy = {'3-d2_bromo_WT','4-d2_bromo_T174M-C220L','9-d2_bromo_V130F', ...
%     '11-d2_brc_L214M','13-d2_brc_L92G','15-d2_brc_F217M'};

% thisSysLabel = 'abcdefghijklmnop'; % Can you sing the abc?
% thisSysLabel = {'WT','T174M-C220L','V130F','L214M','L92G','F217M'};
thisSysLabel = {'WT','T174M-C220L'};
name = 'dd2_dpa_new';


% Chains in this order: [receptor, G protein, ligand]
% Missing chains can be set as '-'.
chains = 'ACB';

abc = char(65:65+25); % Can you sing the abc?

% Analysis parameters:
nPipes = 10; % number of pipelines to consider in analysis
topHubs = 25; % Number of top scoring global hubs to consider
isGPCR = true;
% Generic GPCR numbering from https://gpcrdb.org/residue/residuetable
gpcrdbRefName = 'DD2R';

% Reference PDB code that this simulation is based on
pdbCode = '6VMS';
pdbChains = 'RA-';

% Inactive state reference PDB, used for dihedral comparison
pdbCodeInactive = '6CM4';
pdbInactiveChains = 'A--';

md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc_Culled_data';
%% Make the meta-direcory
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end
md2pathdirMeta = [metadir slash 'md2pathMeta' slash];
% md2pathdir = [metadir slash 'md2pathMeta' slash name slash];
md2pathdir = [metadir slash 'md2pathMeta' slash];

% Create meta directory
if ~exist(md2pathdirMeta, 'dir')
    mkdir(md2pathdirMeta)
    add2log(md2pathdirMeta, ['Creating md2path directory in ' metadir]);
end

% Create subdirectories for every bundle of systems studied
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir)
    add2log(md2pathdir, ['Creating md2path directory in ' md2pathdirMeta]);
end

%% Load, fetch and align PDB files

metaLoadAlign;
%%
refentryNumber = length(foldersToStudy) + 1;

allhubRes = []; % Net hubscore for each allosteric hub (residues with non-zero hubscore); 
%Columns are: Residue no., Hubscore

allhubResGpBarr = []; % G protein and Arrestin hubscores for each residue
allhubEC2IC = []; % Hubscores for pathways connecting EC to IC (Gp or Barr)
channelinfo = cell(length(foldersToStudy),1); % Info on allosteric channels
gpiResCell = cell(length(foldersToStudy),1); % Gp interacting res
bsResCell = cell(length(foldersToStudy),1); % BS residues
helicesCell = cell(length(foldersToStudy),1); % Good old helices

% MI per residue: 
% For a given residue, MIperresidue is the average MI of other residues at 
% close and far distances
MIperresidue = cell(length(foldersToStudy),1);

for thisSys = 1:length(foldersToStudy)
    mydir = ([metadir slash foldersToStudy{thisSys} slash md2pathName slash alloPathCalcName]); 
%     mydir = ([metadir slash foldersToStudy{thisSys} slash]); 
    temp =  importdata([mydir slash 'allhubResidues.txt']);
    allhubRes = [allhubRes temp(1:topHubs,:)]; 
    allhubResGpBarr = [allhubResGpBarr importdata([mydir slash 'allhubResiduesAll.txt'])];
    allhubEC2IC = [allhubEC2IC importdata([mydir slash 'gpi_bai_hubs.txt'])];
%     allhubEC2IC = [allhubEC2IC];
    channelinfo{thisSys} = importdata([mydir slash 'hubinfo.dat']);
    gpiResCell{thisSys} = importdata([mydir slash 'GPI_residues.txt']);
    bsResCell{thisSys} = importdata([mydir slash 'BS_residues.txt']);
    MIperresidue{thisSys} = importdata([mydir slash 'MIperresidue.txt']);
    if isGPCR 
        helicesCell{thisSys} = importdata([mydir slash 'helices.txt']);
    end
end
resList = 1:size(allhubResGpBarr,1); % residue list, assumes all compared 
% receptors are the same!

xlsName = [md2pathdir 'hubResidues_' name '.xls',];
% Preallocate table names
% Globalhubs table
varTypesPerSys = ["string","double","double"];
varTypes = repmat(varTypesPerSys,1,length(foldersToStudy));
% GPI BAI
tempType = repmat("double",1,length(foldersToStudy)*2);
varTypes2 = ["double","string", tempType,"double" ,"string","double"];

% Prepare the labels:
tableLabels =[];
% tableSysLabel = cell(1,length(thisSysLabel)*2);
tableSysLabelSheet1 =[];
tableSysLabel =[];
tableLabelsGpBarr = {'Residue', 'Residue (Ref)'};
for thisSys = 1:length(foldersToStudy)
    tableLabels = [tableLabels "Residue (Ref)"+ num2str(thisSys) "Residue"+num2str(thisSys)  "Hubscore"+num2str(thisSys)];
    tableSysLabelSheet1 =  [tableSysLabelSheet1 thisSysLabel{thisSys} {''} {''}];
    tableSysLabel =  [tableSysLabel thisSysLabel{thisSys} {''}];
    tableLabelsGpBarr = [tableLabelsGpBarr {['Gpi score' num2str(thisSys)],['Barr score' num2str(thisSys)]}];
end
tableLabelsGpBarr = [tableLabelsGpBarr {'Sum Gpi score','Residue (ref)', 'residue'}];


sz = [topHubs length(foldersToStudy)*3];  % Global hubs table
sz2 = [length(resList) length(foldersToStudy)*2+5]; %GpBarr table
globalHubsTab = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableLabels);
gpibaiTab =  table('Size',sz2,'VariableTypes',varTypes2,'VariableNames',tableLabelsGpBarr);


% fill the tables with data!
for thisSys = 1:length(foldersToStudy)
    sysTblNdx = length(varTypesPerSys)*(thisSys-1)+1;
    sysArrNdx = 2*(thisSys-1)+1;

    globalHubsTab(:,sysTblNdx) =  table(arrayfun(@(x) refChain.formatResidues(x,'SecondaryEntry',database.entries{refentryNumber},'BWonly', true),allhubRes(:,sysArrNdx),'UniformOutput',false));
    globalHubsTab(:,sysTblNdx+1) = array2table(allhubRes(:,sysArrNdx));
    globalHubsTab(:,sysTblNdx+2) = array2table(allhubRes(:,sysArrNdx+1));

end

gpibaiTab(:,1) = array2table(resList');
gpibaiTab(:,2) = table(arrayfun(@(x) refChain.formatResidues(x,'SecondaryEntry',database.entries{refentryNumber},'BWonly', true),resList'));
gpibaiTab(:,end-1) = table(arrayfun(@(x) refChain.formatResidues(x,'SecondaryEntry',database.entries{refentryNumber},'BWonly', true),resList'));
gpibaiTab(:,end) = array2table(resList');

gpibaiTab(:,3:end-2) = array2table([allhubResGpBarr sum(allhubResGpBarr,2)]);
ec2gpTab =  gpibaiTab;
ec2gpTab(:,3:end-2) = array2table([ allhubEC2IC sum(allhubEC2IC,2) ]);


writecell(tableSysLabelSheet1,xlsName,'Sheet', 'Global Hubs','Range','A1')
writetable(globalHubsTab,xlsName,'Sheet', 'Global Hubs','Range','A2')
    
writecell(tableSysLabel,xlsName,'Sheet', 'GPI BAI Hubs','Range','C1')
writetable(gpibaiTab,xlsName,'Sheet', 'GPI BAI Hubs','Range','A2')

writecell(tableSysLabel,xlsName,'Sheet', 'EC2GP-BA Hubs','Range','C1')
writetable(ec2gpTab,xlsName,'Sheet', 'EC2GP-BA Hubs','Range','A2')

% % Save top global hubs:
% xlsName = [md2pathdir 'hubResidues_' name '.xls',];
% tableLabels =[];
% % tableSysLabel = cell(1,length(thisSysLabel)*2);
% tableSysLabelSheet1 =[];
% tableSysLabel =[];
% tableLabelsGpBarr = {'Residue'};
% % Prepare the labels:
% for thisSys = 1:length(foldersToStudy)
%     tableLabels = [tableLabels {'Residue','Residue (Ref)','Hubscore'}];
%     tableSysLabelSheet1 =  [tableSysLabelSheet1 thisSysLabel{thisSys} {''} {''}];
%     tableSysLabel =  [tableSysLabel thisSysLabel{thisSys} {''}];
%     tableLabelsGpBarr = [tableLabelsGpBarr {'Gpi score','Barr score'}];
% end
% tableLabelsGpBarr = [tableLabelsGpBarr {'Sum Gpi score','Residue'}];
% 
% if exist(xlsName,'file')
%     add2log(md2pathdir, ['Warning!! ' xlsName ' already exists! Will append to it' ]);
% end
% writecell(tableSysLabel,xlsName,'Sheet', 'Global Hubs','Range','A1')
% writematrix(allhubRes,xlsName,'Sheet', 'Global Hubs','Range','A3')
% writecell(tableLabels,xlsName,'Sheet', 'Global Hubs','Range','A2')
%     
% writecell(tableSysLabel,xlsName,'Sheet', 'GPI BAI Hubs','Range','B1')
% writematrix([resList' allhubResGpBarr sum(allhubResGpBarr,2) resList'],xlsName,'Sheet', 'GPI BAI Hubs','Range','A3')
% writecell(tableLabelsGpBarr,xlsName,'Sheet', 'GPI BAI Hubs','Range','A2')
% 
% writecell(tableSysLabel,xlsName,'Sheet', 'EC2GP-BA Hubs','Range','B1')
% writematrix([resList' allhubEC2IC sum(allhubEC2IC,2) resList'],xlsName,'Sheet', 'EC2GP-BA Hubs','Range','A3')
% writecell(tableLabelsGpBarr,xlsName,'Sheet', 'EC2GP-BA Hubs','Range','A2')
%% Divide allosteric pipelines to different parts of the protein and save 
% them to a spreadsheet

% In the pipe variable, columns re the following:
% channel number, residue number, hubscore, AvgMI
pipe_ligandCell = cell(length(foldersToStudy),1);
pipe_gpCell = cell(length(foldersToStudy),1);
pipe_tmCell = cell(length(foldersToStudy),1);
pipe_loopCell = cell(length(foldersToStudy),1);
    
xlsName = [md2pathdir 'pipes_BS_TM_Gp_' name '.xls',];
tableLabels = {'Channel';'Residue';'Residue (Ref)';'Hubscore';'AvgMI'};
varTypes = ["double","double","string","double","double"];

if exist(xlsName,'file')
    add2log(md2pathdir, ['Warning!! ' xlsName ' already exists! Will append to it' ]);
end
for thisSys = 1:length(foldersToStudy)
    pipe_ligand = [];
    pipe_gp = [];
    pipe_tm = [];
    pipe_loop =[];
    helicesList =[];
    if isGPCR % Construct full list of helix residues
        for thisHelix = 1:7
            helicesList = [helicesList [helicesCell{thisSys}(thisHelix,1):helicesCell{thisSys}(thisHelix,2)]];
        end
    end
    
    % Analyze pipelines in terms of ligand binding, TM or Gp contacting
    for pipe = 1:nPipes
       pipeRows = find(channelinfo{thisSys}(:,1)==pipe);

       for i =1:length(pipeRows)
           row = pipeRows(i);
           thisRes = channelinfo{thisSys}(row,5);

           if find(bsResCell{thisSys}==thisRes) % Is it ligand binding?
               temp = [channelinfo{thisSys}(row,[1 5 6]) MIperresidue{thisSys}(thisRes)];
               pipe_ligand = [ pipe_ligand ; temp]; 
           elseif find( gpiResCell{thisSys} ==thisRes) % Is it gp binding?
               temp = [channelinfo{thisSys}(row,[1 5 6]) MIperresidue{thisSys}(thisRes)];
               pipe_gp = [ pipe_gp ; temp];
           elseif ismember(thisRes,helicesList) && isGPCR % Is it in tm region? 
               % ONLY RELEVANT FOR GPCRS!
               temp = [channelinfo{thisSys}(row,[1 5 6]) MIperresidue{thisSys}(thisRes)];
               pipe_tm = [ pipe_tm ; temp];
               
           else % Other/loop regions
               temp = [channelinfo{thisSys}(row,[1 5 6]) MIperresidue{thisSys}(thisRes)];
               pipe_loop = [ pipe_loop ; temp];
           end
       end
    end
    pipe_ligandCell{thisSys} = pipe_ligand;
    pipe_gpCell{thisSys} = pipe_gp;
    pipe_tmCell{thisSys} = pipe_tm;
    pipe_loopCell{thisSys} = pipe_loop;

    % Write the arrays into tables for ease of storage:
    sz = [size(pipe_ligand,1) 5];
    pipesLigTab = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableLabels);
    pipesLigTab(:,[1 2 4 5]) = array2table(pipe_ligand);
    pipesLigTab(:,3) = table(arrayfun(@(x) refChain.formatResidues(x,'SecondaryEntry',database.entries{refentryNumber},'BWonly', true),pipe_ligand(:,2)));

    sz = [size(pipe_gp,1) 5];
    pipesGpTab = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableLabels);
    pipesGpTab(:,[1 2 4 5]) = array2table(pipe_gp);
    pipesGpTab(:,3) = table(arrayfun(@(x) refChain.formatResidues(x,'SecondaryEntry',database.entries{refentryNumber},'BWonly', true),pipe_gp(:,2)));
    
    sz = [size(pipe_tm,1) 5];
    pipesTMTab = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableLabels);
    pipesTMTab(:,[1 2 4 5]) = array2table(pipe_tm);
    pipesTMTab(:,3) = table(arrayfun(@(x) refChain.formatResidues(x,'SecondaryEntry',database.entries{refentryNumber},'BWonly', true),pipe_tm(:,2)));
    
    sz = [size(pipe_loop,1) 5];
    pipesLoopTab = table('Size',sz,'VariableTypes',varTypes,'VariableNames',tableLabels);
    pipesLoopTab(:,[1 2 4 5]) = array2table(pipe_loop);
    pipesLoopTab(:,3) = table(arrayfun(@(x) refChain.formatResidues(x,'SecondaryEntry',database.entries{refentryNumber},'BWonly', true),pipe_loop(:,2)));
    
    % Write the excel sheet!
%     xlsPosition = 1+ (size(pipe_ligand,2)+1)*((1:length(foldersToStudy))-1);
    xlsPosition = 1+ (size(pipe_ligand,2)+2)*(0:3);
    writecell({'Ligand binding'},xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(1)) '1'])
    writetable( pipesLigTab,xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(1)) '2'])
%     writematrix(pipe_ligand,xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(1)) '3'])
%     writecell(tableLabels',xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(1)) '2'])

    writecell({'Transmembrane'},xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(2)) '1'])
    if isGPCR
        writetable( pipesTMTab,xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(2)) '2'])
%         writematrix(pipe_tm,xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(2)) '3'])
    end
%     writecell(tableLabels',xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(2)) '2'])

    writecell({'G-protein binding'},xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(3)) '1'])
    writetable( pipesGpTab,xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(3)) '2'])
%     writematrix(pipe_gp,xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(3)) '3'])
%     writecell(tableLabels',xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(3)) '2'])
    
    writecell({'Other/loops'},xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(4)) '1'])
    writetable( pipesLoopTab,xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(4)) '2'])
%     writematrix(pipe_loop,xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(4)) '3'])
%     writecell(tableLabels',xlsName,'Sheet', thisSysLabel{thisSys},'Range',[abc(xlsPosition(4)) '2'])
end

add2log(md2pathdir, ['Wrote allosteric pipelines to pipes_BS_TM_Gp_' name '.xls' ]);

% a=table(pipe_tm(:,1), pipe_tm(:,2),pipe_tm(:,3),pipe_tm(:,4), 'VariableNames',tableLabels)

%% Extract the stuff that we need for comparison: ASSUMES 1ST SYSTEM IS REFERENCE
% Examples: 
% Mutation position hubscore
% Mutation neighborhood hubscores
% Conserved motifs score:
% CWxP 6.47-6.50
% NPxxY 7.49-7.53
% DRY 3.49-3.51
% PIF 3.40 5.50 6.44
%                                                                           

conserved = {["6.47" "6.48" "6.49" "6.50"], ["7.49" "7.50" "7.51" "7.52" "7.53"] ,...
    ["3.49" "3.50" "3.51"], ["3.40" "5.50" "6.44"]};
% Ligand binding res:  pipe_ligandCell
% GP binding scores: pipe_gpCell

resMotifs = cell(length(conserved),1);
for motif = 1:length(conserved)
    for resi = 1:length(conserved{motif})
        resMotifs{motif} = [ resMotifs{motif} database.findResidue(conserved{motif}(resi))];
    end
end
    
mutPosDB = cell(length(foldersToStudy),1);
mutations = cell(length(foldersToStudy),1);
% neighborhoodRes =  cell(length(foldersToStudy),1);
mutScore = [];
mutScoreWT = [];
mutNeighborScore = []; 
mutNeighborScoreWT = []; 
motifScores = zeros(length(conserved),length(foldersToStudy));

for thisSys = 1:length(foldersToStudy)
    % Mutations (either find automatically or enter manually)

    if exist(fullfile(databasePath, gpcrdbRefName + ".fasta"),'file') % Quick and dirty if statement, change later
        fastaPath = fullfile(databasePath, gpcrdbRefName + ".fasta");
        [mutations{thisSys}, mutPos] = database.findMut(thisSys, Chains.receptor, fastaPath);
        mutations{thisSys}
        mutPos
        
        mutPosDB{thisSys} = zeros(size(mutPos));
        neighborhoodRes =  cell(length(mutPosDB),1);

        if ~isempty(mutPosDB{thisSys}) % Mutant
            for j = 1:length(mutPosDB{thisSys}) % For every mutation
                mutPosDB{thisSys}(j) = find(database.residues{Chains.receptor}{:,length(foldersToStudy)+1}==mutPos(j));
    
                % Find mutation neighborhood, 5A heavy atom cutoff    
                [~, contact_res] = proteinContacts(database.entries{thisSys}.pdb,  ...
                    database.entries{thisSys}.crd , mutPosDB{thisSys}(j), chains(Chains.receptor), chains(Chains.receptor),4);
            
                neighborhoodRes{j} = find(sum(contact_res{1})); % Residue numbering may be problematic!
                mutNeighborScore = [mutNeighborScore sum(allhubResGpBarr(neighborhoodRes{j} ,2*thisSys-1))]; 
                mutNeighborScoreWT = [mutNeighborScoreWT sum(allhubResGpBarr(neighborhoodRes{j} ,1))]; 
            end
        else % no mutant!
            mutNeighborScore = [mutNeighborScore 0];
            mutNeighborScoreWT = [mutNeighborScoreWT 0];
        end
    end


    if ~isempty(mutPosDB{thisSys}) % Mutant
        mutScore = [mutScore allhubResGpBarr(mutPosDB{thisSys},2*thisSys-1)'];
        mutScoreWT = [mutScoreWT allhubResGpBarr(mutPosDB{thisSys},1)'];
    else % No mutations
        mutScore = [mutScore 0]; 
        mutScoreWT = [mutScoreWT 0]; 
    end

    % Lig and Gp binding:
    sum(pipe_ligandCell{thisSys}(:,3))
    sum(pipe_gpCell{thisSys}(:,3))
    
    % Conserved motifs
    for motif = 1:length(conserved)
        
        motifScores(motif,thisSys) = sum(allhubResGpBarr(resMotifs{motif}(thisSys,:),2*thisSys-1));
    end
end

% database.findResidue("3.50")


%% Write output to excel sheet:
xlsName = [md2pathdir 'conserved_motifs_' name '.xls',];
labels = {'Mut pos', 'Mut pos WT','Mut pos Neighborhood', 'Mut pos neighborhood WT','CWxP','NPxxY','DRY','PIF','Ligand binding','Gp binding'};
writecell(labels',xlsName,'Range','A2');


xlsXPos = 2;
% Write out mutation scores
writematrix( [mutScore; mutScoreWT; mutNeighborScore; mutNeighborScoreWT],xlsName,'Range',[abc(xlsXPos) '2'])
for thisSys = 1:length(foldersToStudy)
    % Label system
    writematrix(thisSysLabel{thisSys},xlsName,'Range',[abc(xlsXPos) '1'])
    % Write out mutation scores
%     writematrix( [mutScore(thisSys); mutScoreWT(thisSys); mutNeighborScore(thisSys); mutNeighborScoreWT(thisSys)],xlsName,'Range',[abc(xlsXPos) '2'])
    % Write out conserved position scores
    writematrix(motifScores(:,thisSys),xlsName,'Range',[abc(xlsXPos) '6'])
    % Write out Lig and Gp scores
    writematrix([sum(pipe_ligandCell{thisSys}(:,3)) ;sum(pipe_gpCell{thisSys}(:,3))],xlsName,'Range',[abc(xlsXPos) '10'])
    xlsXPos = xlsXPos + max(1,length(mutPosDB{thisSys})); % to account to more than 1 mutation
end
%% Compare what we have with experimental data (if available):

