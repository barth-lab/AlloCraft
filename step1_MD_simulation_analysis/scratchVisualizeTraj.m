% Two main problems to solve here:
% 1- visualizing trajectory in a 3D plot if possible 
% 2- comparison of two trajectories residue wise to figure out difference
% in motion

settings.frames2skip = 100;
atomIndices = mainChain.getAtoms(); % Grab CA
traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
[~, trajMean] = meanstructure(traj);




%% Visualize some trajectory elements


stride = 50;

figure
% Initial structure
% scatter3(trajMean(1,3*(1:length(atomIndices))-2), trajMean(1:stride:end,3*(1:length(atomIndices))-1), trajMean(1:stride:end,3*(1:length(atomIndices))),5,'filled')
plot3(mainEntry.pdb.xyz(atomIndices,1),mainEntry.pdb.xyz(atomIndices,2),mainEntry.pdb.xyz(atomIndices,3), ...
    'Color','k','LineWidth',1)
hold on
scatter3(mainEntry.pdb.xyz(atomIndices,1),mainEntry.pdb.xyz(atomIndices,2),mainEntry.pdb.xyz(atomIndices,3), ...
    15, 'o','MarkerFaceColor','K','MarkerFaceAlpha',0.5)

for i = 200:215%length(atomIndices)
scatter3(trajMean(1:stride:end,3*i-2), trajMean(1:stride:end,3*i-1), trajMean(1:stride:end,3*i),25,'.')
hold on
end

axis equal


%% Bubble chart maybe?

% Here's a good post about visualizing scatters:
% https://www.mathworks.com/matlabcentral/answers/802966-binning-a-3d-scatter-plot

% Put points into 3D bins; xyzBinNum is an nx3 matrix containing
% the bin ID for n values in xyz for the [x,y,z] axes.
xyz = [trajMean(1:stride:end,3*i-2) trajMean(1:stride:end,3*i-1) trajMean(1:stride:end,3*i)];
nBins = 8;  % number of bins
xbins = linspace(min(xyz(:,1)),max(xyz(:,1))*1,nBins+1);
ybins = linspace(min(xyz(:,2)),max(xyz(:,2))*1,nBins+1);
zbins = linspace(min(xyz(:,3)),max(xyz(:,3))*1,nBins+1);
xyzBinNum = [...
    discretize(xyz(:,1),xbins), ...
    discretize(xyz(:,2),ybins), ...
    discretize(xyz(:,3),zbins), ...
    ];
% bin3D is a mx3 matrix of m unique 3D bins that appear 
% in xyzBinNum, sorted.  binNum is a nx1 vector of bin
% numbers identifying the bin for each xyz point. For example,
% b=xyz(j,:) belongs to bins3D(b,:).
[bins3D, ~, binNum] = unique(xyzBinNum, 'rows');
% density is a mx1 vector of integers showing the number of 
% xyz points in each of the bins3D. To see the number of points
% in bins3D(k,:), density(k).  
density = histcounts(binNum,[1:size(bins3D,1),inf])'; 
% Compute bin centers
xbinCnt = xbins(2:end)-diff(xbins)/2;
ybinCnt = ybins(2:end)-diff(ybins)/2;
zbinCnt = zbins(2:end)-diff(zbins)/2;

% Plot bubblechart3
fig = figure();
bubblechart3(...
    xbinCnt(bins3D(:,1)), ...
    ybinCnt(bins3D(:,2)), ...
    zbinCnt(bins3D(:,3)), ...
    density, ...
    density, ...
    'MarkerFaceAlpha', .3, ...
    'MarkerEdgeAlpha', .3)
title('bubblechart3')
%% Compare test and ref simulations coordinate wise?

% Ref traj 
atomIndices = refChain.getAtoms(); % Grab CA
trajRef = refSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
[~, trajMeanRef] = meanstructure(trajRef);

%% Visualize what we're comparing
figure
i=1;
stride = 5;

scatter3(trajMean(1:stride:end,3*i-2), trajMean(1:stride:end,3*i-1), trajMean(1:stride:end,3*i),25,'.')
hold on
scatter3(trajMeanRef(1:stride:end,3*i-2), trajMeanRef(1:stride:end,3*i-1), trajMeanRef(1:stride:end,3*i),25,'.')


%% Naive comparison: difference of means, good for unimodal distributions


meanDiff = mean(trajMean)  - mean(trajMeanRef);


figure;
% Initial structure, plot reference structure to be source of vectors
% p3 = plot3(mainEntry.pdb.xyz(atomIndices,1),mainEntry.pdb.xyz(atomIndices,2),mainEntry.pdb.xyz(atomIndices,3), ...
%     'Color','k','LineWidth',1);
% 
% row = dataTipTextRow('Residue',mainChain.formatResidues(mainChain.resIds));
p3 = plot3(refEntry.pdb.xyz(atomIndices,1),refEntry.pdb.xyz(atomIndices,2),refEntry.pdb.xyz(atomIndices,3), ...
    'Color','k','LineWidth',1);

row = dataTipTextRow('Residue',refChain.formatResidues(refChain.resIds));

p3.DataTipTemplate.DataTipRows(end+1) = row;

hold on

% quiver3(mainEntry.pdb.xyz(atomIndices,1),mainEntry.pdb.xyz(atomIndices,2),mainEntry.pdb.xyz(atomIndices,3), ...
%     meanDiff(1:3:end-2)', meanDiff(2:3:end-1)', meanDiff(3:3:end)')
quiver3(refEntry.pdb.xyz(atomIndices,1),refEntry.pdb.xyz(atomIndices,2),refEntry.pdb.xyz(atomIndices,3), ...
    meanDiff(1:3:end-2)', meanDiff(2:3:end-1)', meanDiff(3:3:end)')
axis equal
xlabel('X [A]');
ylabel('Y [A]');
zlabel('Z [A]');
title(['Trajectory difference: ' mainEntry.name ' - '  refEntry.name  ])
%% Greedy comparison: 3D binning and density comparison
% Options: compare density with KL or Wasserstein?

%% Smarter comparison, fit Gaussians to every residue coordinate and then what?
