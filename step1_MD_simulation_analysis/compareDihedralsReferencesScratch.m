
feature(1,:) % residue numbers for desired conserved residues
dihedralsDatabase = database.calcDihedral(feature, Chains.receptor); % Dihedrals

mainSim.dihedrals

% What to do with preiodicity of dihedrals? 
% 1- zeroStretchtotwopi the whole thing (including ref. dihedrals)
% 2- compare cos/sin rather than raw values
% 3- calcDistMatPBC2D BAAAAM

% How to assign states?
% 1- simple "tol" and standard dev based method (calcState.m)
% 2- cluster the states with dbscan or sth?

% How to visualize? 
% 1- time series 
% 2- scatter (ramachandran plot, chi1-chi2 plots, etc ...)
%%
feature = zeros(length(database.entries),length(resImp));
for i=1:length(resImp)
    featureHere = database.findResidue(resImp(i));
    feature(:,i) = featureHere;
end
dihedralsDatabase = database.calcDihedral(feature, Chains.receptor);


options.ColorMap = [];

dihLabels = ["\phi","\psi","\chi_1","\chi_2","\chi_3","\chi_4","\chi_5"];
temp = mainSim.dihedrals(feature(1,:));
resFormatted = mainChain.formatResidues(feature(1,:));
dihi = 3;
dihj = 4;

figure
tiledlayout('flow')
for resNdx = 1:length(resImp)
    nexttile

    if size(temp{resNdx},2) < max(dihi,dihj) % skip if dihedral does not exist
        continue
    end
    Xdat = temp{resNdx}(:,dihi);
    Ydat = temp{resNdx}(:,dihj);
    if ~isempty(options.ColorMap)
        scatter(Xdat,Ydat,5,options.ColorMap,'filled')
    else
        scatter(Xdat,Ydat,5,'filled')
    end
    hold on;
    xlabel(dihLabels(dihi))
    ylabel(dihLabels(dihj))
    handles = zeros(1, length(database.entries));
    
    for entryIndex = 1:length(database.entries)
    if length(dihedralsDatabase{entryIndex,resNdx}) < max(dihi,dihj) % skip if dihedral does not exist
        continue
    end
    if ~isnan(dihedralsDatabase{entryIndex,resNdx}) % nan when alignment is missing from this entry
        h = scatter(dihedralsDatabase{entryIndex,resNdx}(dihi),dihedralsDatabase{entryIndex,resNdx}(dihj), ...
              100,'o','MarkerFaceAlpha',0.75 , 'DisplayName', database.entries{entryIndex}.name);
    else
        h = scatter(nan,nan, ...
              100,'o','MarkerFaceAlpha',0.75 , 'DisplayName', database.entries{entryIndex}.name);
    end

          set(h, 'MarkerFaceColor', get(h, 'MarkerEdgeColor'));
          set(h, 'MarkerEdgeColor', 'w');
        
          handles(entryIndex ) = h;
    end
    handles(handles==0)=[]
    legend(handles,'Location','best');
    legend boxoff;
    title(resFormatted(resNdx))
    
    axis([0 2*pi 0 2*pi])
end
sgtitle(name)
%% Clustering angles with periodicity:
resNdx = 5;
Xdat = temp{resNdx}(:,dihi);
Ydat = temp{resNdx}(:,dihj);


% figure % Visualize data on a unit circle
% scatter(cos(Xdat), sin(Xdat),5,'filled')
% 
% hold on 
% 
% scatter(cos(Ydat), sin(Ydat),5,'filled')
% legend('dih1','dih2')
% axis equal
tic
% DistMat with periodic boundary conditions
distMat = zeros(length(Xdat));
for pt1 = 1:length(Xdat)-1

    for pt2 = pt1+1:length(Xdat)
        dx = abs(Xdat(pt2) - Xdat(pt1));
        if dx > pi
            dx = 2*pi - dx;
        end
        dy = abs(Ydat(pt2) - Ydat(pt1));
        if dy > pi
            dy = 2*pi - dy;
        end
        distMat(pt1,pt2) = sqrt(dx^2 + dy^2);
        distMat(pt2,pt1) = distMat(pt1,pt2) ;
    end 
end
toc

tic
[distMat] = calcDistMatPBC2D(Xdat,Ydat,2*pi,2*pi);
toc
%% Try GMMs:

% Preallocation
k = 1:5;
nK = numel(k);
gm = cell(nK,1);         
aic = zeros(nK,1);
bic = zeros(nK,1);
converged = false(nK,1);

clusterDat = [cos(Xdat), sin(Xdat) ,cos(Ydat), sin(Ydat)];
% Find optimal number of clusters:
for k = 1:nK
    gm{k} = fitgmdist( clusterDat,k);
    aic(k) = gm{k}.AIC;
    bic(k) = gm{k}.BIC;
    converged(k) = gm{k}.Converged;
end
allConverge = (sum(converged(:)) == nK)


gm = fitgmdist( clusterDat,3);
idx = cluster(gm,clusterDat);
cluster1 = (idx == 1); % |1| for cluster 1 membership
cluster2 = (idx == 2); % |2| for cluster 2 membership

figure
gscatter(Xdat,Ydat,idx)
% legend('Cluster 1','Cluster 2','Location','best')

%% DBscan?
epsilon = 0.25;
minpts = 50;
% [idx,corepts] = dbscan(clusterDat,epsilon,minpts)
[idx,corepts]= dbscan(distMat,epsilon,minpts,'Distance','precomputed')

figure
gscatter(Xdat,Ydat,idx)
title(['\epsilon = ' num2str(epsilon) ', minpts = ' num2str(minpts)])