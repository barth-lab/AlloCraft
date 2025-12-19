%% kldivMain: script that runs the KL divergence calculations and analysis

%% Open issues:
% 1- Benchmarking using experimental data or other software
% 2- Interpretability of the divergences, what does a KLDiv of a 0.5 mean?
% 3- Double check histogramming, do we need to zeroStretchtotwopi the
% dihedrals before inputting them?

%% Initialize variables:
global settings;

if ~isfield(settings,'frames2skipRef') 
    settings.frames2skipRef=settings.frames2skip;
end

if ~isfield(settings,'includeLigPathwayCalc') 
    settings.includeLigPathwayCalc = false;
end

if ~isfield(settings,'pdbCodeExtra') 
    settings.pdbCodeExtra = [];
end

klDivFileName = ['klDiv_' settings.refName  '.mat'];
%% Load, fetch and align PDB files and trajectories
database = Database(settings.databasePath);

% Chains: [receptor, G protein, ligand]
database.read(fullfile(settings.mydir, "prot_pymol.pdb"), settings.chains, settings.mainName);
% Load the KLDiv reference structure:
database.read(fullfile(settings.refdir, "prot.pdb"), settings.chains, settings.refName );

% Read reference structures
if ~isempty(settings.pdbCode)
    database.fetch(settings.pdbCode, settings.pdbChains);

    if ~isempty(settings.pdbCodeInactive)
        database.fetch(settings.pdbCodeInactive, settings.pdbInactiveChains);
    end
    for i = 1:length(settings.pdbCodeExtra)
        database.fetch(settings.pdbCodeExtra{i}, settings.pdbChainsExtra{i});
    end
end

% Align sequences and produce the database.residues table
database.align();

% Align structures to the first database entry
database.alignStructures();

% Store names for each chain
Chains.receptor = 1;
Chains.gprotein = 2;
Chains.ligand = 3;

% Keep references to entries and chains used frequently
mainEntry = database.entries{1};
refEntry = database.entries{2}; % This entry spot usually reserved for active ref pdb

mainChain = mainEntry.chains{Chains.receptor};
refChain = refEntry.chains{Chains.receptor};

if mainEntry.hasChain(Chains.ligand)
    ligandChain = mainEntry.chains{Chains.ligand};
end

% Transform xtcs to dcds if they're not there yet:
dirHere = settings.mydir;
for j=1:2
    areThereDCDs = dir(fullfile(dirHere, "run*", "traj.dcd"));
    
    if length(areThereDCDs) < settings.numRuns
        % run VMD from command line:  vmd -dispdev text -e
        [scriptPath, ~, ~] = fileparts(mfilename('fullpath'));
        pathToScript = fullfile(scriptPath, "load_save.tcl");  % assumes script is in the same directory
        cmdStr       = "vmd -dispdev text -e " + pathToScript + " -args " + dirHere + " " + num2str(settings.numRuns) + " " + num2str(settings.stride) + " " + database.entries{j}.path + " " + settings.xtcName;
        system(cmdStr);
%         add2log(md2pathdir, "Transformed .xtc files to .dcd files with the following command: """ + cmdStr + """");
    end
    dirHere = settings.refdir;
end
% Load trajs into sim class
mainSim = mainEntry.addSimulation(fullfile(settings.mydir, "run*"));
refSim = refEntry.addSimulation(fullfile(settings.refdir, "run*"));

md2pathdir = fullfile(settings.mydir, "md2pathdev");

% Create md2pathdir if it does not exist
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir);
end
%% Add labels

% Add labels with Ballesteros-Weinstein notation, aligned with the 2nd
% database entry
% Exported from GPCRdb (https://gpcrdb.org/residue/residuetable)
if length(database.entries) > 1 && isfield(settings, 'systemName') && settings.isGPCR
    residueTablePath = fullfile(database.dir, settings.systemName + "_residue_table.xlsx");
    database.label(3, residueTablePath); % Read 3rd entry (supposedly reference)
end
%% Calculate dihedrals from trajs, assumes same receptor IDs between entries!

% Test system
if settings.includeLigPathwayCalc
    receptorResIds = mainChain.concatResIds(ligandChain);
    receptorResIdsRef = refChain.concatResIds(ligandChain);
else
    receptorResIds = mainChain.resIds;
    receptorResIdsRef = refChain.resIds;
end

% if settings.alignTestRefDih % When Test and ref sims have different residue numbers/length
%     receptorResIds = database.residues{1}.output1(find(database.residues{1}.output2));
% end

mainSim.computeDihedrals( ...
    'Path', fullfile(md2pathdir, "dihedrals.mat"), ...
    'ReSortPath', fullfile(md2pathdir, "reSort.txt"), ...
    'ResIds', receptorResIds, ...
    'StartFrame', settings.frames2skip + 1 ...
);

% Reference system
refSim.computeDihedrals( 'Path', fullfile(settings.refdir, "md2pathdev","dihedrals.mat"),...
    'ReSortPath', fullfile(md2pathdir, "reSortRef.txt"), ...
    'ResIds', receptorResIdsRef, ...
    'StartFrame', settings.frames2skipRef  + 1 ...
);
%% Dihedral union:
% Now that we have dihedral time series, find the union of the two
% dihedral lists to input into the KL divergence function

 [reSortCommon, reSortCommonRef] = mainSim.reconcileDihedralList(refSim,'Path', fullfile(md2pathdir, "reSortKLDiv.mat"));
% Still a problem in deconcile dihedrals when resort is renumbered but the
% PDB is not

% Problem when one PDB starts at 1 but the other doesn't AND there is
% length difference
 %% Calculate KL divergence:
 
 % First order first
 [kl1] = computeKLDiv(refSim.dihedralsMat(:,mainSim.klDivStuff.reSortndxRef), ...
     mainSim.dihedralsMat(:,mainSim.klDivStuff.reSortndxHere), settings.nBlocks, ...
     'Grassberger', settings.Grassberger, 'SignificanceThreshold', settings.st);
 mainSim.klDivStuff.kl1 = kl1;
 save(fullfile(md2pathdir, klDivFileName),'kl1')
 %% 2nd order if you want: (Hint: I don't want)
%  
%  if exist(fullfile(md2pathdir, "klDiv.mat"),'file')
%     load(fullfile(md2pathdir, "klDiv.mat"));
%      else
%     [~, kl2, m2, pVal, pVal2, kl2H0] = computeKLDiv(refSim.dihedralsMat(:,mainSim.klDivStuff.reSortndxRef), ...
%      mainSim.dihedralsMat(:,mainSim.klDivStuff.reSortndxHere), settings.nBlocks, ...
%      'Grassberger', settings.Grassberger, 'SignificanceThreshold', settings.st);
%     save(fullfile(md2pathdir, "klDiv.mat"),'kl1','kl2','m2','pVal','pVal2','kl2H0')
%      end
% mainSim.klDivStuff.kl2 = kl2;
% mainSim.klDivStuff.m2 = m2;

%% Find top kl1 dihedrals and plot their histograms

% Temp dihedrals, please delete later
dihedralsRef = refSim.dihedralsMat(:,mainSim.klDivStuff.reSortndxRef);
dihedralsTest = mainSim.dihedralsMat(:,mainSim.klDivStuff.reSortndxHere);

dihCount = size(dihedralsRef, 2);

xDih = 1:dihCount;
figure 
plot(kl1)
if exist('pVal','var')
    hold on ; scatter(xDih(pVal<settings.st),kl1(pVal<settings.st),'filled')
end

[~,kl1ndx] = sort(kl1,'descend');

figure('Renderer', 'painters', 'Position', [10 10 1500 750]) ;
tiledlayout(3,3)
for i =1:9
    nexttile
    binLimitsHere = binLimits(dihedralsRef(:,kl1ndx(i)),dihedralsTest(:,kl1ndx(i)));
    histogram(dihedralsTest(:,kl1ndx(i)),50, 'BinLimits', binLimitsHere, 'Normalization', 'pdf')
    hold on 
    histogram(dihedralsRef(:,kl1ndx(i)),50, 'BinLimits', binLimitsHere, 'Normalization','pdf')

    % Find the resname and dihType:
    dihType = reSortCommon(kl1ndx(i),2);
    if dihType == 1
        dihText = '\phi';
    elseif dihType == 2
        dihText = '\psi';
    else
        dihText = '\chi';
    end

    resMain = (mainEntry.pdb.resname(mainEntry.pdb.resseq==reSortCommon(kl1ndx(i),1),:));
    resRef = (refEntry.pdb.resname(refEntry.pdb.resseq==reSortCommonRef(kl1ndx(i),1),:));
    textMain = ['Test: ' resMain(1,:) num2str(reSortCommon(kl1ndx(i),1)) ' ' dihText];
    textRef = ['Ref: ' resRef(1,:) num2str(reSortCommonRef(kl1ndx(i),1)) ' ' dihText];
    legend(textMain,textRef,'Location','best')
    legend boxoff
end

figPath = fullfile(md2pathdir,['kl1_dihedrals_' mainEntry.name '_' refEntry.name]);
savefig(figPath + ".fig");
print2pdf(figPath+ ".pdf");
% clear dihedralsRef dihedralsTest

%%  Find KLDiv residue-wise
% Visualize it on the structure

kl1Res = zeros(length(receptorResIds),1);
kl1BB_SC = zeros(length(receptorResIds),2); % Column 1 is BB, 2nd is SC KL1

receptorResIdsNdx = zeros(max(receptorResIds),1); % Helpful for indexing
receptorResIdsNdx(receptorResIds)=1:length(receptorResIds); % Hoping this would do the trick

isSC = reSortCommon(:,2)==0; % List of side chain dihedrals
isBB =  reSortCommon(:,2)~=0; % List of backbone dihedrals
for i = 1:length(receptorResIds)
    resHere = receptorResIds(i);
    kl1Res(i) = sum(kl1(reSortCommon(:,1) == resHere));

    kl1BB_SC(i,1) = sum(kl1(reSortCommon(:,1) == resHere & isBB)); % Fill in backbone KL1
    kl1BB_SC(i,2) = sum(kl1(reSortCommon(:,1) == resHere & isSC)); % Fill in sidechain KL1

end
save(fullfile(md2pathdir, klDivFileName),'kl1','kl1Res');



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

figure('Renderer', 'painters', 'Position', [10 10 1200 600]) 
% plot residue-wise KLDiv

% plot(receptorResIds,kl1Res,'LineWidth',1)
% legend_entries{1} = 'KLDiv';
% hold on
% legend_count = 2;

b = bar(receptorResIds, kl1BB_SC,'stacked'); % Use receptor IDs as X? 
b(1).FaceColor =  [0.6500 0.2250 0.0980]; %[0.6350 0.0780 0.1840]
b(2).FaceColor = [.3 .75 .6]; %[0.3010 0.7450 0.9330];

% Add residue name and numbers to the plot: 

row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds,'BWonly', true));
b(1).DataTipTemplate.DataTipRows(end+1) = row;
b(2).DataTipTemplate.DataTipRows(end+1) = row;

legend_entries{1} = 'Backbone KL1';
legend_entries{2} = 'Sidechain KL1';

hold on
legend_count = 3;

% Add ligand binding residues
if exist(fullfile(md2pathdir,'BS_residues.txt'),'file')
    bsRes.main = importdata(fullfile(md2pathdir,'BS_residues.txt'));
    scatter(bsRes.main,kl1Res(receptorResIdsNdx(bsRes.main)), 'LineWidth',1.5)
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
    scatter(bsResRefRenum,kl1Res(bsResRefRenum), 25, 'filled')
    legend_entries{legend_count} = 'Ref binding residues';
    legend_count = legend_count + 1;
end

if ~isempty(mutations)
    s3 = scatter(mutRes, kl1Res(receptorResIdsNdx(mutRes)),80, 'LineWidth',1.5);
    text(mutRes,kl1Res(receptorResIdsNdx(mutRes)) +0.05*max(kl1Res),mutations)
    row = dataTipTextRow('Mut',mutations);
    s3.DataTipTemplate.DataTipRows(end+1) = row;

    legend_entries{legend_count} = 'Mutation sites';
    legend_count = legend_count + 1;
end

if ~isempty(settings.highlightRes)
    scatter(settings.highlightRes, kl1Res(receptorResIdsNdx(settings.highlightRes)),60,'x', 'LineWidth',1.5);
    legend_entries{legend_count} = settings.highlightText;
    legend_count = legend_count + 1;
end

legend(legend_entries,'Location','best')
legend boxoff
title(['KL1 of ' mainEntry.name '|' refEntry.name])
% Draw helices if helices are defined
if ~isempty(settings.helices)
    drawTMhelices(kl1Res, settings.helices, receptorResIds)
end
xlabel('Residue'); ylabel('KL Divergence')
set(gca,'FontSize',16)

figPath = fullfile(md2pathdir,['kl1_' mainEntry.name '_' refEntry.name]);
savefig(figPath + ".fig");
print2pdf(figPath+ ".pdf");
%% Save the KLDiv as the B-factor
bfactor = zeros(1, mainEntry.atomCount);


if settings.includeLigPathwayCalc % Include ligand in b-factor array
    chainSelect = selectname(mainEntry.pdb.chainid, settings.chains(Chains.receptor)) | selectname(mainEntry.pdb.chainid, settings.chains(Chains.ligand));
    selection = selectname( mainEntry.pdb.name,'CA') & chainSelect;
else
    selection = selectname( mainEntry.pdb.name,'CA') & selectname(mainEntry.pdb.chainid, settings.chains(Chains.receptor));
end
% resFromReSort = selectid(mainEntry.pdb.resseq,unique(reSortCommon(:,1))); % Only residues contained in reSortCommon
% selection = selection & resFromReSort;

bfactor(selection) = kl1Res;
writepdb(fullfile(md2pathdir, ['prot_kl1div_' mainEntry.name '_' refEntry.name '.pdb']), mainEntry.pdb, [], 'default', bfactor);

% 
% bfactor(selection) = kl1BB_SC(:,1);
% writepdb(fullfile(md2pathdir, "prot_kl1div_BB.pdb"), mainEntry.pdb, [], 'default', bfactor);
% 
% bfactor(selection) =  kl1BB_SC(:,2);
% writepdb(fullfile(md2pathdir, "prot_kl1div_SC.pdb"), mainEntry.pdb, [], 'default', bfactor);




        
        


        
        
%% Support functions:

function output = binLimits(xref,x) 
    % Creates common bin limits for two dihderal series while keeping them
    % in the [ 0 2*pi] limit
    minCommon = min(min(x), min(xref));
    maxCommon = max(max(x), max(xref));
    output = [max(minCommon, 0), min(maxCommon, 2*pi)];

end




