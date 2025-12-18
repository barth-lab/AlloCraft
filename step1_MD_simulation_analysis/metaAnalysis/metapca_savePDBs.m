%% Script to save the pdbs of the centers of the clusters of separate systems:
% this is a separate script since it was too awkward to make into a
% function

% CSys and PCs of the ligand cluster centers
CSysClusterCenters = [CSys(ind_centersHere,:)  (1:length(ind_centersHere))'];

 pcaTraj = cell(length(foldersToStudy),1);
 pcaTrajHighestDen = cell(length(foldersToStudy),1);
 pcaClusterTraj = cell(length(ind_centersHere),1);
  for thisSys = 1:length(foldersToStudy)

    myDir = fullfile(metadir,foldersToStudy{thisSys});
    entryHere  = database.entries{thisSys};
    chainHere =  entryHere.chains{Chains.receptor};
    atomIndices = chainHere.getAtoms(); % Grab CA atoms
    simHere = entryHere.addSimulation(fullfile(myDir, "run*"),'align2chain',chains(Chains.receptor));

    % Save pdbs of centers of systems
    for modelN = 1:nCenters
        pcaTraj{thisSys} = [ pcaTraj{thisSys} ;simHere.traj{run_FrameNdxHere{thisSys}(modelN,1)}(run_FrameNdxHere{thisSys}(modelN,2),:)];
        pcaTrajHighestDen{thisSys} = [ pcaTrajHighestDen{thisSys} ;simHere.traj{run_FrameHighestDenNdxHere{thisSys}(modelN,1)}(run_FrameHighestDenNdxHere{thisSys}(modelN,2),:)];
    end
    writepdbndx([md2pathdir pcaName '_' thisSysLabel{thisSys} '.pdb'], entryHere.pdb, [], 0, pcaTraj{thisSys});

    % Save pdbs of highest density points of every system
    writepdbndx([md2pathdir pcaName 'HighestDen_' thisSysLabel{thisSys} '.pdb'], entryHere.pdb, [], 0, pcaTrajHighestDen{thisSys});

    % Save pdbs of PCA cluster centers:

    for thisCenter = 1:length(ind_centersHere)
        if CSysClusterCenters(thisCenter,1) == thisSys % Cluster is in this system, save PDB!
            pcaClusterTraj{thisCenter} = simHere.traj{CSysClusterCenters(thisCenter,2)}(CSysClusterCenters(thisCenter,3),:);
            writepdbndx([md2pathdir pcaName 'Cluster_' thisSysLabel{thisSys} '_C' num2str(thisCenter) '.pdb'], entryHere.pdb, [], 0, pcaClusterTraj{thisCenter});
        end
    end

  end