% Utility script that only does the steps to calculate MI:

database = Database(settings.databasePath);
database.read(fullfile(settings.mydir, "prot.pdb"), settings.chains, "Protein");
md2pathdir = fullfile(settings.mydir, "md2pathdev");
% Create md2pathdir if it does not exist
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir);
end
% Store names for each chain
Chains.receptor = 1;
Chains.gprotein = 2;
Chains.ligand = 3;

% Keep references to entries and chains used frequently
mainEntry = database.entries{1};
mainChain = mainEntry.chains{Chains.receptor};

%% Transform .xtc files into .dcd

% Are the DCD files there? If yes, skip this step
% areThereDCDs = dir([settings.mydir '/run*/traj.dcd']);
areThereDCDs = dir(fullfile(settings.mydir, "run*", "traj.dcd"));

if length(areThereDCDs) < settings.numRuns
    % run VMD from command line:  vmd -dispdev text -e
    pathToScript = fullfile(pwd(), "load_save.tcl");  % assumes script is in current directory
    cmdStr       = "vmd -dispdev text -e " + pathToScript + " -args " + settings.mydir + " " + num2str(settings.numRuns) + " " + num2str(settings.stride) + " " + database.entries{1}.path + " " + settings.xtcName;
    system(cmdStr);
end

% Load sim
mainSim = mainEntry.addSimulation(fullfile(settings.mydir, "run*"),'align2chain',settings.chains(Chains.receptor));

%% Calculate dihedrals

resList = unique(mainEntry.pdb.resseq,'stable');

mainSim.computeDihedrals( ...
    'Path', fullfile(md2pathdir, "dihedrals.mat"), ...
    'ReSortPath', fullfile(md2pathdir, "reSort.txt"), ...
    'ResIds', resList, ...
    'StartFrame', settings.frames2skip + 1 ...
);

%% Calculate MI! (Including ALL frames)

mainSim.computeMI('Path', fullfile(md2pathdir, "MI.mat"));


%% Assess convergence of entropies (Including ALL frames)

assessEntropyConvergence(mainSim, 'SavePath', fullfile(md2pathdir, "convergence"));
add2log(md2pathdir, "S2 convergence calculated and plotted");


%% Filter MI
MI = filterCorrectMI(mainSim.dihedralsMat, mainSim.mi, 'SaveDir', md2pathdir,'SignificanceThreshold',0.05);

%%

resind = mainSim.reSort(:,1);
dihtype = mainSim.reSort(:,2);
dihedral = mainSim.reSort(:,4:5);

Nres = length(resList);

MIres=zeros(Nres);

for i=1:Nres-1
        for j=i+1:Nres
            temp = MI(resind==(i),resind==(j));
            temp = temp(:);
            MIres(i,j) = sum(temp);
            MIres(j,i) = MIres(i,j);
        end
end


figure;imagesc(MIres) 
axis xy; axis square;
xlabel('residue','fontsize',40);
ylabel('residue','fontsize',40);
colorbar; colormap turbo
formatplot2

if ~isempty(settings.helices)
    drawTMhelices(1:Nres,settings.helices,1:Nres)
    drawTMhelices(1:Nres,settings.helices,1:Nres,'Y')
end

savefig(fullfile(md2pathdir,"MIres"))

save(fullfile(md2pathdir,"workspace.mat"),'MIres','-append');

MIresPrint = [resList MIres ];
MIresPrint = [[0 resList'] ; MIresPrint];

writematrix(round(MIresPrint,3), fullfile(md2pathdir,"MIres.dat"),'Delimiter','tab');
%% Write out summed MI as b-factors

MIresSum = sum((MIres)); % Sum the MI over every residue with given residue

bfactor = zeros(1, mainEntry.atomCount);
bfactor(mainEntry.getAtoms('Name','CA')) = sqrt(MIresSum);
writepdb(fullfile(md2pathdir, "prot_MIres.pdb"), mainEntry.pdb, [], 'default', bfactor);
