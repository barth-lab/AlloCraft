% Trajectory distribution differences (trajDD)
% Comparison of two trajectories residue wise to figure out difference
% in motion



simResList =  database.residues{1}{:,1:2}; % list from simulation set
allAligned = sum(simResList==0,2)==0;
alignmentVector = simResList(allAligned,:); % alignmentVector works if:
% 1- residue numbers start from 1
% 2- residue numbers of same chain are sorted in increasing order (no funny 
% PDB files where residue 500 is before residue 1 in the same chain)

% Concatenate trajs:
% test traj
atomIndices = mainChain.getAtoms(); % Grab CA
traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
[~, trajMean] = meanstructure(traj);

% Ref traj 
atomIndicesRef = refChain.getAtoms(); % Grab CA
trajRef = refSim.concatRuns('Atoms', atomIndicesRef, 'StartFrame', settings.frames2skip + 1);
[~, trajMeanRef] = meanstructure(trajRef);

clear traj trajRef
%% Naive comparison: difference of means, good for unimodal distributions
% Here's my quick attempt to fix alignment issues:
meanDiff = mean(trajMean(:,to3(alignmentVector(:,1))))  - mean(trajMeanRef(:,to3(alignmentVector(:,2))));

figure;
% Initial structure, plot reference structure to be source of vectors

p3 = plot3(refEntry.pdb.xyz(atomIndicesRef,1),refEntry.pdb.xyz(atomIndicesRef,2),refEntry.pdb.xyz(atomIndicesRef,3), ...
    'Color','k','LineWidth',1);
row = dataTipTextRow('Residue',refChain.formatResidues(refChain.resIds));
p3.DataTipTemplate.DataTipRows(end+1) = row;

hold on
s3 = scatter3(refEntry.pdb.xyz(atomIndicesRef,1),refEntry.pdb.xyz(atomIndicesRef,2),refEntry.pdb.xyz(atomIndicesRef,3), ...
    15,refChain.resIds, 'o','filled');
row = dataTipTextRow('Residue',refChain.formatResidues(refChain.resIds));
s3.DataTipTemplate.DataTipRows(end+1) = row;

% quiver3(mainEntry.pdb.xyz(atomIndices,1),mainEntry.pdb.xyz(atomIndices,2),mainEntry.pdb.xyz(atomIndices,3), ...
%     meanDiff(1:3:end-2)', meanDiff(2:3:end-1)', meanDiff(3:3:end)')
quiver3(refEntry.pdb.xyz(atomIndicesRef,1),refEntry.pdb.xyz(atomIndicesRef,2),refEntry.pdb.xyz(atomIndicesRef,3), ...
    meanDiff(1:3:end-2)', meanDiff(2:3:end-1)', meanDiff(3:3:end)','r')
axis equal
xlabel('X [A]');
ylabel('Y [A]');
zlabel('Z [A]');
title(['Trajectory difference: ' mainEntry.name ' - '  refEntry.name  ])

figPath = "trajDiffMean_" +  mainEntry.name  + '_' + refEntry.name ;
savefig(figPath)
% Reverse Z?
% set(gca,'zdir','reverse')

% Plot test structure too?
% p3 = plot3(mainEntry.pdb.xyz(atomIndices,1),mainEntry.pdb.xyz(atomIndices,2),mainEntry.pdb.xyz(atomIndices,3), ...
%     'Color',[ 0.9290 0.6940 0.1250 0.5],'LineWidth',1);
% row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds));
% p3.DataTipTemplate.DataTipRows(end+1) = row;

%% Visualize motion of residues of interest as scatters:
roi = [40 152 178 205 261]-2; % roi depends on alignmentVector!!!

stride = 50;
hold on
% scatter3(mainEntry.pdb.xyz(atomIndices,1),mainEntry.pdb.xyz(atomIndices,2),mainEntry.pdb.xyz(atomIndices,3), ...
%     15, 'o','MarkerFaceColor','K','MarkerFaceAlpha',0.5)
defaultColors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];

for i = roi
scatter3(trajMean(1:stride:end,3*alignmentVector(i,1)-2), trajMean(1:stride:end,3*alignmentVector(i,1)-1), trajMean(1:stride:end,3*alignmentVector(i,1)),25,defaultColors(1,:),'.')
hold on
scatter3(trajMeanRef(1:stride:end,3*alignmentVector(i,2)-2), trajMeanRef(1:stride:end,3*alignmentVector(i,2)-1), trajMeanRef(1:stride:end,3*alignmentVector(i,2)),25,defaultColors(2,:),'.')
end

axis equal
legend([refEntry.name ' struct'],'', 'Mean \Delta traj', 'Test traj', 'Ref traj ')

figPath = fullfile(md2pathdir,['trajDD_' mainEntry.name '_' refEntry.name]);
savefig(figPath + ".fig");
% clear trajMean trajMeanRef
%% Additional visualization options, specific atom choices other than CA
targetRes = 125;
targetAtom = 'CB';

ndx = selectid(mainEntry.pdb.resseq,targetRes) & selectname(mainEntry.pdb.name,targetAtom);
ndxRef = selectid(refEntry.pdb.resseq,targetRes) & selectname(refEntry.pdb.name,targetAtom);


traj = mainSim.concatRuns('Atoms', find(ndx), 'StartFrame', settings.frames2skip + 1);
trajRef = refSim.concatRuns('Atoms', find(ndxRef), 'StartFrame', settings.frames2skip + 1);

for i = 1:length(targetRes)
scatter3(traj(1:stride:end,3*i-2), traj(1:stride:end,3*i-1), traj(1:stride:end,3*i),25,defaultColors(1,:),'.')
hold on
scatter3(trajRef(1:stride:end,3*i-2), trajRef(1:stride:end,3*i-1), trajRef(1:stride:end,3*i),25,defaultColors(2,:),'.')
end
%% Greedy comparison: 3D binning and density comparison
% Options: compare density with KL or Wasserstein?

%% Smarter comparison, fit Gaussians to every residue coordinate and then what?
