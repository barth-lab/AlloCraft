%% Import the data:
% THIS IS A LEGACY SCRIPT TO BE DELETED
% Choice of data to import will be made before calling the script. 

% Initialize traj cells:
traj = cell(max(numRuns),1);
    
% Read the PDBs (useful for selections)

[pdb, crd ] = readpdb([mydir '/' pdbName]); 
% readdcdmat is from the MATDCD package, and is more reliable than readdcd
% from the mdtoolbox package

minTraj = 0;
maxTraj = 0;
for runi=1:numRuns
    % Read the trajectories:
    traj{runi} = readdcdmat([mydir '/run' mat2str(runi) '/traj.dcd']);
    temp = decenter(traj{runi});
    traj{runi} = temp;
    if runi == 1
        minTraj = size(traj{runi},1);
    end
    minTraj = min(minTraj,size( traj{runi},1));
    maxTraj = max(maxTraj,size( traj{runi},1));
end
clear temp

% If trajectories have different lengths and the difference in lengths is
% small, just make them all same length

if (maxTraj-minTraj)/minTraj < 0.01 % Less than 1% difference
   for runi=1:numRuns
        traj{runi} = traj{runi}(1:minTraj,:);
   end
   
else
    add2log(md2pathdir,{'Different runs have more than 1% difference in length, may affect RMSD/RMSF calculations.',''});
end
% Make sure PDB and traj has the same number of atoms
assert(size(traj{1},2)/3 == length(pdb.serial),'PDB and trajectory have different number of atoms!!!')
%% Align/superimpose the structures:

% SIDE NOTE: I tested aligning with and without helices and the difference
% is minimal, so I'll drop it for now
% Choose helices for superimposing (alignment): 
% D2 active state model
% helices_to_align = [2:30 36:65 77:106 117:142 151:192 199:235 240:278];

% D1 model:
% helices_to_align = [4:30 39:67 72:107 119:142 159:202 210:245 258:291];
traj_aligned = cell(max(numRuns),1);

    
% Get indices for CAs and CBs, regardless of chain
index_CA_all = selectname(pdb.name, 'CA');
index_CB_all = selectname(pdb.name, 'CB');
% chain?
index_chain = selectname(pdb.chainid,chains(1));

for runi=1:numRuns
    [~, traj_aligned{runi}] = superimpose(crd, traj{runi}, find(index_CA_all & index_chain));
end

traj = traj_aligned; % Shorter var name
clear traj_aligned; % Clean memory from unused trajectory/ies

% Don't forget to log!
add2log(md2pathdir,{'Imported and aligned trajectories successfully',''});

%% Useful selections to be used later:


% Get some diagnostics on the PDB and trajs
[nFrames,nAtoms] =  size(traj{1});
nAtoms = nAtoms/3; % Traj has every atom three times for x,y, and z

Ures = unique(pdb.resseq,'stable'); % List of ALL residues in the pdb (all chains)
noH = ~selectname(pdb.name, 'H*'); % Heavy atoms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using Ures when different chains have common residue IDs is a BAD idea
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NresAll = length(Ures); % Number of ALL residues in the pdb file (all chains)

% Divide them by chain
chains = unique(pdb.chainid); % 1st chain is receptor, 2nd is ligand, 3rd is Gp
protChain = chains(1);
protRes = unique(pdb.resseq(pdb.chainid==protChain)); % Gives protein residues

if length(chains)>1 % Make ligand (and Gp) selections
    ligandChain = chains(2);
%     gpChain = chains(3);

    ligandRes = unique(pdb.resseq(pdb.chainid==ligandChain)); % Gives ligand residues
%     gpRes = unique(pdb.resseq(pdb.chainid==gpChain)); % Gives protein residues
    index_chainL = selectname(pdb.chainid,ligandChain);
    % Is ligand small molecule or peptide?
    if length(ligandRes)>1
      isLigSmall = 0;  
      index_L = index_chainL & index_CA_all;
      add2log(md2pathdir, 'Ligand is a peptide!')
    else
      isLigSmall = 1; 
      index_L = index_chainL & noH ;
      add2log(md2pathdir, 'Ligand is a small molecule!')
    end
end

% Grab some useful selections:
index_chainR =  selectname(pdb.chainid, protChain);
index_CA = index_chainR & index_CA_all; % CAs of receptor only 
index_CB = index_chainR & index_CB_all; % CBs of receptor only 