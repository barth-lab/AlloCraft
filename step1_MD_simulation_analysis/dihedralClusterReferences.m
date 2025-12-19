function [idxCell] = dihedralClusterReferences(database,resImp,options)
% dihedralClusterReferences Clusters dihedral pairs in a 2D map using
% DBscan and assigns "states" using reference structures in the database
%
%% Usage:
% dihedralClusterReferences(database,resImp)
% [idxCell] = dihedralClusterReferences(database,resImp,options)
% * database: database class from AlloDy containing:
%   1- a main entry with dihedrals calculated using 
% simulation.computeDihedrals 
%   2- (Optional) inactive and active reference pdbs for assignment of
%   states
%
% * resImp: list of important residues to choose. Right now the function
% only supports BW numbering for GPCRs. Ex.
% resImp = ["3.50","5.54","6.47","6.48","7.49" ,"7.53"];
% 
% * idxCell: cell structure containing cluster indices for every resImp 
%
% * options: 
% 'dihij': which dihedral pair to plot? must be input in the format: 
% [dihi dihj]. Dihedrals are numbered from 1 to 7, with [1 2] being
% backbone dihedrals and [3:7] being side chain dihedrals. Defaults to 
% [1 2]
% 
% For other options, check the code, don't be lazy


arguments
    database
    resImp
    options.mainEntryNdx = 1; % First entry is used as reference
    options.receptorChain = 1; % Assumes receptor is 1st chain
    options.dihij = [1 2]; % Dihedrals to plot
    options.epsilon = 0.25; % For DBscan calculation
    options.minpts = 50; % For DBscan
end

% define simulation, chain, and reference states
mainSim = database.entries{options.mainEntryNdx}.simulation;
mainChain = database.entries{options.mainEntryNdx}.chains{options.receptorChain};

actRefNdx = find(database.entryState==1);
actRefNdx = actRefNdx(1); % Does not support mutliple active/inactive states
inactRefNdx = find(database.entryState==2);
inactRefNdx = inactRefNdx(1);

% Get features from list of important residues and dihedral list
feature = zeros(length(database.entries),length(resImp));
for i=1:length(resImp)
    featureHere = database.findResidue(resImp(i));
    feature(:,i) = featureHere;
end

dihedralsDatabase = database.calcDihedral(feature, options.receptorChain);
dihLabels = ["\phi","\psi","\chi_1","\chi_2","\chi_3","\chi_4","\chi_5"];

assert(~isempty(mainSim.dihedralsMat),'Simulation does not contain dihedrals! Calculate dihedrals first!!!')
temp = mainSim.dihedrals(feature(1,:));
resFormatted = mainChain.formatResidues(feature(1,:));
idxCell = cell(length(resImp),1);
% Pick dihedrals to plot:
dihi = options.dihij(1);
dihj = options.dihij(2);

% Plot desired dihedrals with their time series
figure
tiledlayout('flow')
for resNdx = 1:length(resImp)
    nexttile

    if size(temp{resNdx},2) < max(dihi,dihj) % skip if dihedral does not exist
        continue
    end
    Xdat = temp{resNdx}(:,dihi);
    Ydat = temp{resNdx}(:,dihj);
    
    if actRefNdx && inactRefNdx % Add active and inactive references to time series
        Xdat = [ Xdat; dihedralsDatabase{actRefNdx,resNdx}(dihi) ;dihedralsDatabase{inactRefNdx,resNdx}(dihi)];
        Ydat = [ Ydat; dihedralsDatabase{actRefNdx,resNdx}(dihj) ;dihedralsDatabase{inactRefNdx,resNdx}(dihj)];
    end
    % Get distance matrix between dihedrals with periodic boundary
    % conditions
    [distMat] = calcDistMatPBC2D(Xdat,Ydat,2*pi,2*pi);

    % Cluster elements
    [idx,~]= dbscan(distMat,options.epsilon,options.minpts,'Distance','precomputed');
    idxCell{resNdx} = idx;
    
    if actRefNdx && inactRefNdx % Find how many elements in active and inactive like cluster:
        if idx(end-1) == idx(end) % Active and inactve are in the same clusters
            stateLabel = "Ambiguous state";
        else
            activeLike = sum(idx(1:end-3) == idx(end-1))/length(temp{resNdx}(:,dihi));
            inactiveLike = sum(idx(1:end-3) == idx(end))/length(temp{resNdx}(:,dihi));
            stateLabel = "Act=" + num2str(activeLike,3) + ", Ina=" + num2str(inactiveLike,3) ;
        end
    else
        stateLabel = "";
        warning("Could not find active and inactive reference states in database!!!")
    end

    gscatter(Xdat,Ydat,idx)
    hold on;
    xlabel(dihLabels(dihi))
    ylabel(dihLabels(dihj))
    handles = zeros(1, length(database.entries));
    
    for entryIndex = 1:length(database.entries)
    if length(dihedralsDatabase{entryIndex,resNdx}) < max(dihi,dihj) % skip if dihedral does not exist
        continue
    end
    if entryIndex == actRefNdx % Active state ref, mention that
        dispName = [database.entries{entryIndex}.name ':Active'];
        mrkr = "^";
    elseif entryIndex == inactRefNdx % Inactive state ref, mention that
        dispName = [database.entries{entryIndex}.name ':Inactive'];
        mrkr = "v";
    else
        dispName = database.entries{entryIndex}.name;
        mrkr = 'o';
    end
    if ~isnan(dihedralsDatabase{entryIndex,resNdx}) % nan when alignment is missing from this entry
        h = scatter(dihedralsDatabase{entryIndex,resNdx}(dihi),dihedralsDatabase{entryIndex,resNdx}(dihj), ...
              100,mrkr,'MarkerFaceAlpha',0.75 , 'DisplayName', dispName);
    else
        h = scatter(nan,nan, ...
              100,mrkr,'MarkerFaceAlpha',0.75 , 'DisplayName', dispName);
    end

          set(h, 'MarkerFaceColor', get(h, 'MarkerEdgeColor'));
          set(h, 'MarkerEdgeColor', 'w');
        
          handles(entryIndex ) = h;
    end
    handles(handles==0)=[];
    legend(handles,'Location','best');
    legend boxoff;
    title(resFormatted(resNdx) + ";" + stateLabel)
    
    axis([0 2*pi 0 2*pi])
end
sgtitle([database.entries{options.mainEntryNdx}.name ', \epsilon = ' ...
    num2str(options.epsilon) ', minpts = ' num2str(options.minpts)])


end