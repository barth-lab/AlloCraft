function [nHBonds,pair_cell,pair_res_cell,hbond_per_res, angles_cell, is_bb_cell] = hydrogenBondContacts(pdb, traj, target ,prot, targetChain, protChain, cut_off_max, cut_off_min, cut_off_angle, tol_res)
%hydrogenBondAnalysis Performs Hydrogen bond analysis of the PDB for every
%frame in traj.
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
%
%% Usage:
% nHBonds = hydrogenBondAnalysis(pdb, traj)
% nHBonds = hydrogenBondAnalysis(pdb, traj, cut_off_max, cut_off_min, cut_off_angle, tol_res)
% [nHBonds,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] = hydrogenBondAnalysis(pdb, traj)
% [nHBonds,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] = hydrogenBondAnalysis(pdb, traj, cut_off_max, cut_off_min, cut_off_angle, tol_res))
%
%% Description:
% * nHBonds are the number of hydrogen bonds detected for every frame of
% the trajectory. [nFrames x 1] array
%
% * pdb is the pdb structure obtained by pdb = readpdb('pdb.pdb'). Note
% that other structure files (.gro for example) will not work, as the way
% mdtoolbox names the structure elements is different.
%
% * traj is the trajectory, as obtained by traj = readdcdmat(traj.dcd) for
% example. [Nframes x 3*Natoms]
%
% * cut_off_max is the maximum cut-off for the pair-list calculations. This
% is the distance between H and the acceptor atom. Defaults to 2.5 Ang.
%
% * cut_off_min is the minimum considered distance for H-bonds, any
% distance less that this cut-off is not considered. Defaults to 1.5 Ang.
%
% * cut_off_angle is the cut-off value for the angle formed between the
% Donor-H--Acceptor, criterion for acceptance is calculated as follows: 
% abs(angle - 180) <= cut_off_angle. Defaults to 30 degrees.
% 
% * tol_res, pairs in the same residue OR within "tol_res" residues are not
% considered as H-bond forming pairs. Defaults to 1, meaning that atoms on
% the following, same, and previous residues are not considered for bonding 
% by default.
%
% * pair_cell, pair_res_cell, pair_names_cell, angles_cell, and is_bb_cell 
% contain the atom pairs, residue pairs, atom names, D-H--A angles, and 1/0 
% for backbone Hbond respectively, as numbered and named in PDB. 
% nFrame x 1 cell structure, where each cell is nHBonds(frame) long.
%
% See also hydrogenBondPeaks, hydrogenBondManipulate

% Set the default values:
if ~exist('cut_off_max','var')
    cut_off_max = 2.5; %  Angstrom
end

if ~exist('cut_off_min','var')
    cut_off_min = 0; % Angstrom
end

if ~exist('cut_off_angle','var')
    cut_off_angle = 30; % Degrees
end

if ~exist('tol_res','var')
   tol_res = 1;
end

chains = unique(pdb.chainid); % Get all chains in PDB
Nchains = length(chains);

if ~exist('protChain', 'var') || isempty(protChain) % Default protChain to 
    % only chain in PDB, if PDB has more than 1 chain, assert an error!   
    assert(Nchains == 1, 'PDB has more than 1 chain, please specify a chain to be the receptor/protein')    
    protChain = unique(pdb.chainid);
end

if ~exist('targetChain', 'var') || isempty(targetChain) % Default targetChain to 
    % only chain in PDB, if PDB has more than 1 chain, assert an error!   
    assert(Nchains == 1, 'PDB has more than 1 chain, please specify a chain to be the target/ligand')    
    targetChain = unique(pdb.chainid);
end 

protChainNum = find(chains == protChain);
Ures = cell(Nchains,1);
NresAll = 0;

for chain = 1:Nchains % List of residues in the pdb by chain
    Ures{chain} = unique(pdb.resseq(selectid(pdb.chainid,chains(chain))),'stable');
    NresAll = NresAll +  length(Ures{chain}); % Number of residues in the protein
end
Nres = length(Ures{protChainNum});

nFrames =  size(traj,1);
nHBonds = zeros(nFrames, 1);

if nargout>1 % Did the user ask for the detailed pair info?
    pair_cell = cell(nFrames, 1);
    pair_res_cell = cell(nFrames, 1);
    pair_names_cell = cell(nFrames, 1);
    angles_cell = cell(nFrames, 1);
    is_bb_cell = cell(nFrames, 1);
end

% Backbone - backbone interactions: Choose C, CA, N, O, HN, and HA only
index_ca = selectname(pdb.name, 'CA');
index_c = selectname(pdb.name, 'C');
index_n = selectname(pdb.name, 'N');
index_o = selectname(pdb.name, 'O');
index_hn = selectname(pdb.name, 'HN');
index_ha = selectname(pdb.name, 'HA');
% For general calculations, choose ALL H, N, and O's
index_Hall = selectname(pdb.name, 'H*');
index_Nall = selectname(pdb.name, 'N*');
index_Oall = selectname(pdb.name, 'O*');

index_bb = index_ca | index_c | index_n | index_o | index_hn | index_ha; % 1 for backbone, 0 for not

% Grab index for target and protein/receptor atoms and their placement in the pdb 
targetIndex = selectid(pdb.resseq,target); % 1's for target, 0's for others
temp =  selectid(pdb.chainid,targetChain); % Grab the proper chain
targetIndex = targetIndex.*temp; % choose target from a specific chain only
targetIndexFind = find(targetIndex);
pdb_resseq_target = pdb.resseq(find(targetIndex)); % contains resseq 
% only for target

% protein
protIndex = selectid(pdb.resseq,prot);
temp = selectid(pdb.chainid,protChain);
protIndex = protIndex.*temp;
protIndexFind = find(protIndex);

hbond_res = zeros(length(target), Nres, nFrames);

for frame = 1:nFrames
    
        % Debugging
    if mod(frame,100) ==0

        frame

    end
    
    crd_target = traj(frame,to3(logical(targetIndex)));
    crd_prot = traj(frame,to3(logical(protIndex)));

     % Calculate all atoms of target within a cutoff distance of protein
     [pair, dists] = searchrange(crd_prot, crd_target, cut_off_max);


    % Remove any bonds that are shorter than cut_off_short (1.5 A), those are not Hbonds:
    pair1 = pair(:,1);
    pair2 = pair(:,2);
    pair1(dists<cut_off_min)=[];
    pair2(dists<cut_off_min)=[];
    pair = [pair1 pair2];

    
    % Transform pair index from crd_prot and crd_target index to pdb index:
    pair1 = targetIndexFind(pair(:,1)); % target is in column 1
    pair2 = protIndexFind(pair(:,2)); % protein is in column 2
    pair = [pair1 pair2];
    
    
    % the section below isn't needed in this framework
    
%     % Remove pairs within the same or within tol_res (def = 1) residues
%     same_res = zeros(size(pair,1),1);
%     for row = 1:size(pair,1)
%         if abs(pdb.resseq(pair(row,1)) - pdb.resseq(pair(row,2))) <= tol_res
%             same_res(row) = 1;
%         end
%     end
%     pair = pair.*not(same_res); % Zeroes out pairs on same or within "tol_res" residues
%     pair1 = pair(:,1);
%     pair2 = pair(:,2);
%     pair1 = nonzeros(pair1);
%     pair2 = nonzeros(pair2);
%     pair = [pair1 pair2];



    % Set rules for accepted donor - acceptor pairs: 
    % Backbone: HN -- O, HA -- O 
    % All: H* -- O*, H* -- N*
    acc_pair = zeros(size(pair,1),1); % Will contain the pairs to accept

    for row = 1:size(pair,1)
        acc_pair(row) = index_Hall(pair(row,1))*index_Oall(pair(row,2)) + ...
            index_Hall(pair(row,2))*index_Oall(pair(row,1)) + ...
            index_Hall(pair(row,1))*index_Nall(pair(row,2)) + ...
            index_Hall(pair(row,2))*index_Nall(pair(row,1));
    end

    pair = pair.*acc_pair; % Zeroes out pairs that do not form H-bonds
    pair1 = pair(:,1);
    pair2 = pair(:,2);
    pair1 = nonzeros(pair1);
    pair2 = nonzeros(pair2);
    pair = [pair1 pair2];

    pair_res = zeros(size(pair,1),2);
    pair_names = repmat('12345678',size(pair,1),1);

    for row = 1:size(pair,1)
        pair_res(row,:) = [pdb.resseq(pair(row,1)) pdb.resseq(pair(row,2))];
        pair_names(row,:) = [pdb.name(pair(row,1),:) pdb.name(pair(row,2),:)];
    end
    
    % Last criterion: angle
    triplet = zeros(size(pair,1),1); % Triplet D-H -- A
    % Hpos is either 2 or 6
    %     ' O   HN '
    %     ' HA  O  '
    for row = 1:size(pair,1) % Rearrange "pair" variable so H is in the second column
        atom_name = pdb.name(pair(row,1),:);
        atom_name = atom_name(~isspace(pdb.name(pair(row,1),:))); % Remove white spaces
        if(atom_name(1)=='H') % H is in the 1st column 
            triplet(row,1) = pair(row,2); % Switch so H is in the 2nd column
            triplet(row,2) = pair(row,1);
        else % H is in the 2nd column
            triplet(row,1) = pair(row,1);
            triplet(row,2) = pair(row,2); 
        end
        % The atom H is attached to is always before it in
        % the PDB, either -1, -2, or -3: so
        test_temp = 1;
        donor = 1;
        atom_name_test = pdb.name(triplet(row,2) - donor,:);
        atom_name_test = atom_name_test(~isspace(atom_name_test)); % Remove white spaces
        while test_temp==1
            if atom_name_test(1) ~= 'H' % it's not a Hydrogen, we found the donor!!!
                break
            end
            donor = donor + 1;
            atom_name_test = pdb.name(triplet(row,2) - donor,:);
            atom_name_test = atom_name_test(~isspace(atom_name_test)); % Remove white spaces
        end
         triplet(row,3) = triplet(row,2) - donor; % Which will give us the donor atom!
    end

    % FINALLY CALCULATE THE FREAKING ANGLES
    angles = calcangle(traj(frame,:),triplet).*180./pi; % In degrees
    angles_acc = angles;
    % Now remove bonds beyond the angle cutoff: cut_off_angle
    accept_angle = zeros(size(pair,1),1);
    for row = 1:size(pair,1)
        if abs(angles(row)-180) <= cut_off_angle
            accept_angle(row) = 1;
        end
    end
    pair = pair.*accept_angle; % Zeroes out pairs beyond the angle cutoff
    angles_acc = angles_acc'.*accept_angle;
    angles_acc = nonzeros(angles_acc);
    pair1 = pair(:,1);
    pair2 = pair(:,2);
    pair1 = nonzeros(pair1);
    pair2 = nonzeros(pair2);
    pair = [pair1 pair2];
    dists_acc = calcbond(traj(frame,:),pair); % Calculate distances for references

    % Calculate residue numbers again
    pair_res = zeros(size(pair,1),2);
    pair_names = repmat('12345678',size(pair,1),1);

    % Fill the occupancies into occupancy matrix
    for pairi = 1:size(pair,1)
        
        pair_res(pairi,:) = [pdb.resseq(pair(pairi,1)) pdb.resseq(pair(pairi,2))];
        pair_names(pairi,:) = [pdb.name(pair(pairi,1),:) pdb.name(pair(pairi,2),:)];
        
        % use protIndex to figure out where the atom is in the PDB
        res = pair_res(pairi,2);% This residue is found in this interaction
        resNdx_prot = find(Ures{protChainNum}==res);

        % Fill in the residue contact matrix:

        res_tar = pair_res(pairi,1);% This residue is found in this interaction
        resNdx_tar = find(target==res_tar);
        hbond_res(resNdx_tar,resNdx_prot,frame) = 1;

    end 
        
    nHBonds(frame) = size(pair,1); % Number of H-bonds per frame
    
    if nargout>1 % Did the user ask for the detailed pair info?
    pair_cell{frame} = pair;
    pair_res_cell{frame} = pair_res;
    pair_names_cell{frame} = pair_names;
    angles_cell{frame} = angles_acc;
%     is_bb_cell{frame} = is_bb;

    % calculate percentage of hbond occupancy:
    hbond_per_res = zeros(length(target), Nres); 
        for res_tar = 1:length(target) % fill contact map residue-wise for target
           for res = 1:Nres % now do residues for protein
                   dist_here = hbond_res(res_tar,res,:);
                   dist_here = reshape(dist_here,[nFrames,1]);
                   % fill in the contact percentage
                   hbond_per_res(res_tar,res) = 1-length(find(dist_here==0))/length(dist_here);
            end 
        end
    end
end
end
