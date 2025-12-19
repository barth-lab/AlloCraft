% Dimensionality reduction of protein coordinates scratch 

%% tICA
% tICA finds high-autocorrelation linear combinations of the input degrees
% of freedom, which are the slowest relaxing ones in a sense

% tICA is usally combined with markov state models to extract relaxation
% times, such as in: http://docs.markovmodel.org/lecture_tica.html#PerezEtAl-JCP13-TICA
% Author of MDtoolbox https://doi.org/10.1063/1.5019750
% First application of the method: https://doi.org/10.1063/1.4811489




atomIndices = mainChain.getAtoms();
traj = mainSim.concatRuns('Atoms', atomIndices, 'StartFrame', settings.frames2skip + 1);
[~, trajMean] = meanstructure(traj);

timeLag = 10;
 [p, mode, lambda] = calctica(trajMean, timeLag);

 figure
 tiledlayout('flow')
 nexttile
scatter(p(:, 1), p(:, 2),5,C(:,1),'filled');
 xlabel('tICA 1', 'fontsize', 25);
 ylabel('tICA 2', 'fontsize', 25);
 nexttile
scatter(mode(:, 1), mode(:, 2),5,'filled');
 xlabel('Mode 1', 'fontsize', 25);
 ylabel('Mode 2', 'fontsize', 25);
clear trajMean
 %% PCA:
 % PCA finds high-variance linear combinations of the input degrees of 
 % freedom

 % PCA on cartesian coords
 pc = calcpca(traj);

 % PCA on pairwise distances
 % Step 1: calculate pairwise distances for desired stuff:
 resDisCalc = 1:4:size(traj,2)/3; % Take every Nth residue
 distMat = zeros(length(resDisCalc),length(resDisCalc),size(traj,1));

 for frameHere = 1:size(traj,1)
    distMat(:,:,frameHere) = calcdistancematrix(traj(frameHere,to3(resDisCalc)));
 end
 % Try PCA on dihedrals
 
 % Backbone
%  pc = calcpca(mainSim.dihedralsMat(:,(mainSim.reSort(:,2)==1 | mainSim.reSort(:,2)==2)));
% 
%  % Sidechain
%  pc = calcpca(mainSim.dihedralsMat(:,(mainSim.reSort(:,2)==0)));
% 
%  % All
%  pc = calcpca(mainSim.dihedralsMat);


%  pc = calcpca(temp');
temp = reshape(distMat,length(resDisCalc)*length(resDisCalc),size(traj,1));
 stdPairs = std(temp');
 [a,b] = sort(stdPairs,'descend');                                          
[coeff,pc,latent,tsquared,explained,mu] = pca(temp(b(1:200),:)');

 figure
scatter(pc(:, 1), pc(:, 2),5,C(:,1),'filled');
 xlabel('PC 1', 'fontsize', 25);
 ylabel('PC 2', 'fontsize', 25);

%% tsne:
% Do pca before?
% perplexity? it's a guess about the number of close neighbors each point has
% distance metric?
% Explanation of tsne: https://distill.pub/2016/misread-tsne/

% [Y,loss] = tsne(trajMean(1:1000,:),'Algorithm','exact','Distance','cosine');
 
 [Y,loss] = tsne(trajMean,'NumPCAComponents',50,'Perplexity',100);
figure
scatter(Y(:, 1), Y(:, 2),5,C(:,1),'filled');
% scatter(Y(:, 1), Y(:, 2),5,'filled');
colormap turbo

%% Diffusion maps?

 %% CV

 path = [0.25 0.3; 0.5 0.6; 0.75 0.3];
data = rand(100000, 2);
[progress, distance] = calcpathcv(path, data);
figure
subplot(1, 2, 1);
scatter(data(:, 1), data(:, 2), 20, progress, 'filled');
hold on;
scatter(path(:, 1), path(:, 2), 200, 'k', 'filled');
 axis square; 
 subplot(1, 2, 2);
 scatter(data(:, 1), data(:, 2), 20, distance, 'filled');
 hold on;
 scatter(path(:, 1), path(:, 2), 200, 'k', 'filled');
 axis square;