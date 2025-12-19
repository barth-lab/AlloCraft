function [dist, contact_res, targetindexHeavy, notHindex, contact_per_res, contact_per] = proteinContactsDualCutoff(pdb, traj, target, targetChain, protChain, rCut1, rCut2, contactCut, allAtom)
%proteinContacts Calculates the contacts between the target and surrounding 
% residues using a dual cutoff scheme. A contact "connects" once 
% dist(i,j) < rCut1 but only disconnects when dist(i,j) > rCut2. Dual
% cutoff works at the atom level rather than at the residue level.
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
% 
%% Usage:
%  [dist] = proteinContacts(pdb, traj, target)
%  [dist, contact_res, targetindexHeavy, notHindex, contact_per_res, contact_per] = proteinContacts(pdb, traj, target, targetChain, protChain, rCut, contactCut, allAtom)
%
%% Description:
% OUTPUT:
% * dist: cell structure containing distance between chosen target atoms and  
%   protein atoms over the trajectory
%
% * contact_res: contacts between target atoms and residues for every frame
% of the simualtion, 1 for contact, 0 for no contact
%
% * targetindexHeavy: index of heavy atoms of the target
%
% * notHindex: index of heavy atoms of the protein
%
% * contact_per_res: % Percentage of simulation in which target and protein 
% RESIDUES are within cutoff
%
% * contact_per:  % Percentage of simulation in which target and protein 
% ATOMS are within cutoff
%
% INPUT:
% * pdb is the pdb structure obtained by pdb = readpdb('pdb.pdb')
%
% * traj is the trajectory in the form of a numRuns x 1 cell structure if 
%   numRuns > 1, where each cell contains Nframes x Natoms coordinates
%   or just a matrix with Nframes x Natoms if numRuns = 1
%
% * target is the residue number of the target in the pdb
%   IF the target is a single residue: contact_res outputs the target atom
%   by atom, while IF target is more than one, contact_res outputs the
%   target residue by residue
%
% * targetChain complements target by specifying a chain, helpful for pdbs
%   with multiple chains and overlapping residue numbers. You can ignore 
%   this option by inputting []  
%
% * protChain specifies which chain in the PDB is to be considered the
% receptor or protein to use as a reference for calculating contacts. this
% one defaults to the only chain if PDB has only 1 chain. This input is
% NECESSARY if PDB has more than one chain
%
% * rCut1 is the cutoff distance in A below which a pair is considered in
% contact
%
% * rCut2 is the cutoff distance in A above which a pair is considered out
% contact
%
% * contactCut: If two residues are below rCut for contactCut fraction 
 % of the simulation or more, they're considered in contact
%
% * allAtom: consider all atoms or just heavy ones? Defaults to false (only
% consider heavy atoms)

%% Define default values:
if ~exist('rCut1', 'var') || isempty(rCut1)
        rCut1 = 3.5; % defaults to 3.5 A
end

if ~exist('rCut2', 'var') || isempty(rCut2)
        rCut2 = 4.5; % defaults to 4.5 A
end

assert(rCut2>=rCut1,'rCut2 is smaller than rCut1, make sure that rCut2 is equal or larger than rCut1');

if ~exist('allAtom', 'var') || isempty(allAtom)
        allAtom = 0; % defaults to considering H atoms
end

if ~exist('contactCut', 'var') || isempty(contactCut)
        contactCut = 0.4; % defaults to 40% of frames within rcut to be considered
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


if iscell(traj)
    numRuns = length(traj);
else % just put traj in a cell of 1
    numRuns = 1;
    temp = cell(1);
    temp{1} = traj;
    traj = temp;
end

%%%%%%%%%%%%%%%%%%%%%

Ures = cell(Nchains,1);
NresAll = 0;

for chain = 1:Nchains % List of residues in the pdb by chain
    Ures{chain} = unique(pdb.resseq(selectid(pdb.chainid,chains(chain))),'stable');
    NresAll = NresAll +  length(Ures{chain}); % Number of residues in the protein
end
Nres = length(Ures{protChainNum});
%%%%%%%%%%%%%%%%%%%%%%%%%%%


targetindex = selectid(pdb.resseq,target); % 1's for target, 0's for others
temp =  selectid(pdb.chainid,targetChain);
targetindex = targetindex.*temp; % choose target from a specific chain only


notHindex = ~selectname(pdb.name, 'H*'); % Remove H's from calculation
targetindexHeavy = targetindex.*notHindex; 
pdb_resseq_target = pdb.resseq(find(targetindexHeavy)); % contains resseq 
% only for target

% Grab index for protein/receptor atoms and their placement in the pdb
protIndex = selectid(pdb.chainid,protChain);

if allAtom == 1
    Natoms = sum(protIndex); % Natoms of the protein/receptor
    protIndexFind = find(protIndex);
%     targetindexHeavy = targetindex; % Bad naming I know
else
    Natoms = sum(protIndex.*notHindex); % Take not H atoms
    protIndexFind = find(protIndex.*notHindex);
end

Nframes = zeros(numRuns,1);
dist = cell(numRuns,1); % each cell contains distance between chosen target 
% atoms and protein atoms over the trajectory
contact_res = cell(numRuns,1); % contacts between target atoms and residues
contact_per = cell(numRuns,1); % Percentage in which target and protein 
%atoms are within cutoff
contact_per_res = cell(numRuns,1); % Percentage in which target and protein 
% residues are within cutoff

for runi=1:numRuns 
    contact_last_frame = zeros(sum(targetindexHeavy), Natoms); % Was there contact in the last frame?
    contact_this_frame = zeros(sum(targetindexHeavy), Natoms); % What about this frame?

    Nframes(runi) = size(traj{runi},1);
    % initialize dist to be a 3D: Nheavy x Natoms_protein x Nframes matrix
    
    % THESE 3D matrices TAKE TOO MUCH memory

    if length(target) == 1 % If target is just one residue, calculate atom-wise,
        contact_res{runi} = zeros(sum(targetindexHeavy), Nres, Nframes(runi));
        dist{runi} = zeros(sum(targetindexHeavy), Natoms, Nframes(runi));
    else % Otherwise, calculate residue wise, and only save dist in this frame
        contact_res{runi} = zeros(length(target), Nres, Nframes(runi));
        dist{runi} = zeros(sum(targetindexHeavy), Natoms);
    end
    for frame = 1:Nframes(runi)
        
        % Debugging
        if mod(frame,100) ==0
            
            frame
            
        end
        
        % Grab the crd of the target/ligand and protein of this run in this frame
         crd_target = traj{runi}(frame,to3(logical(targetindexHeavy)));
         if allAtom == 1
            crd_prot = traj{runi}(frame,to3(logical(protIndex)));
         else
             crd_prot = traj{runi}(frame,to3(logical(protIndex.*notHindex)));
         end
         
         % Search range gives you all pairs of traj that are within cutoff
         % of crd_target
         [pair, dist_list] = searchrange(crd_prot, crd_target, rCut2);


         % Fill the distances into dist matrix
        for pairi = 1:length(dist_list)
            if dist_list(pairi) <= rCut1 || contact_last_frame(pair(pairi,1),pair(pairi,2)) == 1
                if length(target) == 1
                    dist{runi}(pair(pairi,1),pair(pairi,2),frame) = dist_list(pairi);  
                else
                    dist{runi}(pair(pairi,1),pair(pairi,2)) = dist_list(pairi);  
                end
                % use protIndex to figure out where the atom is in the PDB
                res = pdb.resseq(protIndexFind(pair(pairi,2)));% This residue is found in this interaction
                resNdx_prot = find(Ures{protChainNum}==res);
                
                % Fill in the residue contact matrix:
                if length(target) == 1 % fill contact map atom-wise for target
                    contact_res{runi}(pair(pairi,1),resNdx_prot,frame) = 1;
                else % Fill contact map residue-wise for target
                     res_tar = pdb_resseq_target(pair(pairi,1));% This residue is found in this interaction
    %                 res_tar = pdb.resseq(pair(pairi,1));
                    resNdx_tar = find(target==res_tar);
                    contact_res{runi}(resNdx_tar,resNdx_prot,frame) = 1;
                end

                contact_this_frame(pair(pairi,1),pair(pairi,2)) = 1;
                
            end
        end 
        % Now fix the contact_last_frame: pairs that don't show up in the
        % contacts in this frame get reset!
        contact_this_frame(contact_this_frame - contact_last_frame < 0) = 0;
        contact_last_frame = contact_this_frame;
    end
    % If asked for contact percentage, calculate it:
    if nargout > 4
        
        contact_per{runi} = zeros(sum(targetindexHeavy), Natoms);
         if length(target) == 1 % fill contact map atom-wise for target
            contact_per_res{runi} = zeros(sum(targetindexHeavy), Nres);
         else % fill contact map residue-wise for target
            contact_per_res{runi} = zeros(length(target), Nres); 
         end
         
        for target_atom = 1:sum(targetindexHeavy)
           if nargout > 5 % did user ask for contact_per atom-wise?
               for atom = 1:Natoms % Go over atoms first

                   dist_here = dist{runi}(target_atom,atom,:);
                   dist_here = reshape(dist_here,[Nframes(runi),1]);
                   % fill in the contact percentage
                   contact_per{runi}(target_atom,atom) = 1-length(find(dist_here==0))/length(dist_here);
                   % If the percentage is less than cutoff, set it to 0
                   if contact_per{runi}(target_atom,atom) < contactCut
                       contact_per{runi}(target_atom,atom) = 0;
                   end
               end
           end
           
           if length(target) == 1 % fill contact map atom-wise for target
               for res = 1:Nres % now do residues

                   dist_here = contact_res{runi}(target_atom,res,:);
                   dist_here = reshape(dist_here,[Nframes(runi),1]);
                   % fill in the contact percentage
                   contact_per_res{runi}(target_atom,res) = 1-length(find(dist_here==0))/length(dist_here);
                   % If the percentage is less than cutoff, set it to 0
                   if contact_per_res{runi}(target_atom,res) < contactCut
                       contact_per_res{runi}(target_atom,res) = 0;
                   end
               end
           end
        end
        
        if length(target) > 1 % fill contact map residue-wise for target
        for res_tar = 1:length(target)
           for res = 1:Nres % now do residues
                   dist_here = contact_res{runi}(res_tar,res,:);
                   dist_here = reshape(dist_here,[Nframes(runi),1]);
                   % fill in the contact percentage
                   contact_per_res{runi}(res_tar,res) = 1-length(find(dist_here==0))/length(dist_here);
                   % If the percentage is less than cutoff, set it to 0
                   if contact_per_res{runi}(res_tar,res) < contactCut
                       contact_per_res{runi}(res_tar,res) = 0;
                   end
            end 
        end
        end
    end
end
end

