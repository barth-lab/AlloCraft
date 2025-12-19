function [x,y,np] = orderparamsgpcr(pdb,traj,res3,res6,res7,numRuns,pdb_ref,npxxy,npxxy_ref)
%orderparamsgpcr Calculates TM3-6 and TM3-7 distances as well as RMSD for
% select residues from a reference structure.
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
%
%% Usage
% [x,y] = orderparamsgpcr(pdb,traj,res3,res6,res7)
% [x,y] = orderparamsgpcr(pdb,traj,res3,res6,res7,numRuns)
% [x,y,np] = orderparamsgpcr(pdb,traj,res3,res6,res7,numRuns,pdb_ref,npxxy,npxxy_ref)
%
%% Description:
% * x is the TM3-6 between res3 and res6, it is numRuns x 1 cell structure.
%
% * y is the TM3-7 between res3 and res7, it is numRuns x 1 cell structure.
%
% * np is the RMSD between residues npxxy in traj (using pdb to grab residues) 
%   using pdb_ref and npxxy_ref as reference. If npxxy_ref is not given,
%   it defaults to npxxy. Make sure that trajectories in traj are aligned
%   to the reference pdb structure pdb_ref!
%
% * pdb is the pdb structure obtained by pdb = readpdb('pdb.pdb').
%
% * traj is the trajectory in the form of a numRuns x 1 cell structure if 
%   numRuns > 1, where each cell contains Nframes x Natoms coordinates
%   or just a matrix with Nframes x Natoms if numRuns = 1.
%   Default of numRuns = 1, and in that case the function takes traj as a
%   matrix and not a cell structure.
%
% * res3, res6, and res7, are the numbers of the residues used for
%   calculation of the distances res3 - res 6 and res3 - res7.
%
% * numRuns is the number of cells to consider in traj (if traj is a cell
%   structure).
%
% * pdb_ref is the reference pdb structure to measure RMSD against.
%
% * npxxy and npxxy_ref are the residues to be used for RMSD calculation, 
%   if npxxy is not given, it defaults to npxxy. npxxy_ref references back
%   to pdb_ref and is to be used if pdb_ref has a different residue
%   numbering than pdb.

% Set the default value for numRuns
if nargin<6
  numRuns = 1;
end

% Cell structures to hold the data
x = cell(numRuns,1);
y = cell(numRuns,1);
np = cell(numRuns,1);

% Pairs won't change between runs, so calculate them once outside the loop
index_tm36 = selectid(pdb.resseq, [res3 res6]); % Logical array for TM3-6
index_tm37 = selectid(pdb.resseq, [res3 res7]); % Logical array for TM3-7
index_CA = selectname(pdb.name, 'CA'); % Logical array for all CA's
pair_tm36 =  find(index_tm36.*index_CA); % Indices for the TM3-6 pairs
pair_tm37 =  find(index_tm37.*index_CA); % Indices for the TM3-7 pairs

for runi=1:numRuns   
% Calculate the bond lengths
    if numRuns == 1
        x{runi} = calcbond(traj,pair_tm36');
        y{runi} = calcbond(traj,pair_tm37'); 
    else
        x{runi} = calcbond(traj{runi},pair_tm36');
        y{runi} = calcbond(traj{runi},pair_tm37');
    end
end

if nargin >= 7 && nargout >= 3 % If user asked for np RMSD variable
    % Calculate the NPxxY rmsd from reference structure:
    % calcrmsd takes atom number as input for index
    
    % If npxxy for reference is not given, assume it is the same as the trajectory pdb
    if ~exist('npxxy_ref', 'var') || isempty(npxxy_ref)
        npxxy_ref = npxxy;
    end
    
    index_npxxy = selectid(pdb.resseq, npxxy); % Logical array for NPxxY
    index_npxxy_ref = selectid(pdb_ref.resseq, npxxy_ref); 
    % Choose backbone atoms:
    % Surveyed structure:
    index_CA = selectname(pdb.name, 'CA');
    index_C = selectname(pdb.name, 'C');
    index_N = selectname(pdb.name, 'N');
    index_O = selectname(pdb.name, 'O');
    index_backbone = index_CA | index_C | index_N | index_O;
    % Crystal/reference structure:
    index_CA_ref = selectname(pdb_ref.name, 'CA');
    index_C_ref = selectname(pdb_ref.name, 'C');
    index_N_ref = selectname(pdb_ref.name, 'N');
    index_O_ref = selectname(pdb_ref.name, 'O');
    index_backbone_ref= index_CA_ref | index_C_ref | index_N_ref | index_O_ref;
    
    % Calculate the cooridnates in [1 x 3Natoms] format from xyz
    % [Natom x 3] format
    crd_ref = zeros(1,length(pdb_ref.xyz)*3);
    counter = 1;
    for i=1:length(pdb_ref.xyz)
        crd_ref(counter:counter+2) = pdb_ref.xyz(i,:);
        counter = counter + 3;
    end

    % Grab the relevant part of the reference coordinates
    traj_ref = crd_ref(1,to3(index_backbone_ref & index_npxxy_ref));
    
    for runi=1:numRuns
        % Calculate the NPxxY RMSD:
        if numRuns == 1
            traj_rmsd = traj(:,to3(index_backbone & index_npxxy));
        else
            traj_rmsd = traj{runi}(:,to3(index_backbone & index_npxxy));
        end
        np{runi} = calcrmsd(traj_ref,traj_rmsd);
    end
end

end

