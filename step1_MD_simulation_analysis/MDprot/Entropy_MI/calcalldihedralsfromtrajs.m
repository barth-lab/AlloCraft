function [dihedrals, dihIndex, reSort, resname_cell] = calcalldihedralsfromtrajs(pdb,traj,rotamers,numRuns,higherOrder,transform)
%% calcdihedralsfromtrajs: 
% This function takes a set of trajectories and calculates the backbone and
% rotamer dihedrals (Psi), (Phi), and (Chi) of  chosen residues.
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
% Dihedrals based on: http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
%
%% Usage:
% dihedrals = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns);
% dihedrals = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns, 'all',transform);
% dihedrals = calcdihedralsfromtrajs(pdb,traj,rotamers);
% [dihedrals, reSort, resname_cell] = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns);
% [dihedrals, reSort, resname_cell] = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns, 'all',transform);
% [dihedrals, reSort, resname_cell] = calcdihedralsfromtrajs(pdb,traj,rotamers);
%
%% Description
% * pdb is the pdb structure obtained by pdb = readpdb('pdb.pdb'). Note
% that other structure files (.gro for example) will not work, as the way
% mdtoolbox names the structure elements is different.
%
% * traj is the trajectory in the form of a numRuns x 1 cell structure if 
% numRuns > 1, where each cell contains Nframes x Natoms coordinates
% or just a matrix with Nframes x Natoms if numRuns = 1.
% Default of numRuns = 1, and in that case the function takes traj as a
% matrix and not a cell structure.
%
% * rotamers is the list of residues to calculate psi, phi, and chiN for.
%
% * numRuns number of cells to consider in traj (if traj is a cell
% structure).
%
% * dihedrals is a length(rotamers) x numRuns cell structure containing the
% dihedral angles, where each cell has Nframes x N angle. where N = 2 BB +1
% chi or 2 BB + chiN if the option 'all' is given. ( a NaN in the BB slots 
% means the code encountered a C or N term or the dihedral is not defined)
%
% * resname_cell is length(rotamers) x numRuns cell structure containing the
% names of the chosen residues, useful when dealing with homologues that
% may have different residues in the same spot.
%
% * dihIndex is a matrix that indexes present dihedrals (1 for present, 0 for
% not) and it has the following format: 
%    phi psi chi1 chi2 chi3 chi4 chi5
% 
% * 'all' instructs the function to calculate higher order dihedrals chiN
% for residues that have them.
%
% * transform instructs the function to transform the dihedral output from
% from [0,360] to [-180,180] if equal to 1.

% Set the default value for numRuns
if nargin<4
  numRuns = 1;
end

% If higherOrder is not given, do not calculate them
if nargin<5
  higherOrder = 'Non';
end

% If transform is not given, default it to 0
if nargin<6
    transform = 0;
end



dihedrals = cell(length(rotamers),numRuns);
resname_cell = cell(length(rotamers),numRuns);
dihIndex = zeros(length(rotamers),7); % indexes present dihedrals, 7 for 2 BB + 5 Chi
reSort = [];
% number of frames
if numRuns == 1 %traj is a matrix
    nframes = size(traj,1);
    assert(size(pdb.xyz,1) == size(traj,2)/3, 'PDB and trajectory have different number of atoms!')
else            % traj is a cell
    nframes = size(traj{1},1);
    assert(size(pdb.xyz,1) == size(traj{1},2)/3, 'PDB and trajectory have different number of atoms!')
end

counter = 1; % Counter for number of residues
    
% Choose the atoms for the dihedral calculation:
% Backbone atoms are always the same:
% phi: C(-1)-N-Ca-C
% psi: N-Ca-C-N(+1)

atoms_bb{1} = selectname(pdb.name, 'N');
atoms_bb{2} = selectname(pdb.name, 'CA');
atoms_bb{3} = selectname(pdb.name, 'C');

atoms_bb_all = atoms_bb{1} | atoms_bb{2} | atoms_bb{3};
 
% Now chi1:
% The first three atoms are always the same (so we keep them out of the
% loop)
atoms{1} = selectname(pdb.name, 'N');
atoms{2} = selectname(pdb.name, 'CA');
atoms{3} = selectname(pdb.name, 'CB');

% Here comes the rub: (GLY and ALA don't have a chi angle, and 
% their values will be replaced by a NaN)
res_CG = ['ARG '; 'ASN '; 'ASP '; 'GLN '; 'GLU '; 'HIS '; 'HSD '; 'HSE ';...
    'LEU '; 'LYS '; 'MET '; 'PHE '; 'PRO '; 'TRP '; 'TYR '];
res_CG1 = ['ILE '; 'VAL '];
res_OG1 = ['THR '];
res_OG = ['SER '];
res_SG = ['CYS '];
res_break = ['ALA ';'GLY '];
% 
    
if  iscolumn(rotamers) % Make rotamers a row vector so for loop can take it
    % element by element
    rotamers = rotamers';
end

% Do some termini recognition using chainids:
chainChange = find(diff(pdb.chainid)); % find when chain ID changes

% chainChange reports the index of the element BEFORE the chain break
isCterm =  pdb.resseq(chainChange);
isNterm =  pdb.resseq(chainChange+1);

for resNum=rotamers    % Choose a residue from the list of rotamers
    if resNum == 0 % Fake input (useful in tandem with database from md2path)
        dihedrals{counter} = nan;
        counter = counter + 1;
        continue
    end
    index_res = selectid(pdb.resseq, resNum);
    index_resNext = selectid(pdb.resseq, resNum + 1);
    index_resBefore = selectid(pdb.resseq, resNum - 1);
    
    % Deal with backbone dihedrals first:
    % remember: N-term residues do not have phi, C-term res do not have psi
    index_psi = (index_res & atoms_bb_all) | (atoms_bb{1} & index_resNext);
    index_phi =  (atoms_bb{3} & index_resBefore) | (index_res & atoms_bb_all);
    calc_bb = [1 1]; % Defaults to calculate both phi and psi respectively
    resNum
    
    % Recognize termini:
    if resNum == pdb.resseq(1) | ismember(resNum,isNterm) | resNum == 1 |... 
            (abs(resNum - pdb.resseq(min(find(index_res)) - 1)) > 1) % This is most likely an N-term
           calc_bb(1) = 0;    
         
    elseif resNum == pdb.resseq(end) | ismember(resNum,isCterm) | ...
            (abs(pdb.resseq(max(find(index_res)) + 1) - resNum) > 1) % This is most likely an C-term
          calc_bb(2) = 0; 
    end
    
    
    % reSort section!!!!!
    if calc_bb(1)
        temp = find(index_phi);
%         reSort = [reSort; counter 1 temp(2) temp(3)]; % Phi -> 1
    reSort = [reSort; counter 1 temp']; % Phi -> 1
    end
    if calc_bb(2)
        temp = find(index_psi);
%         reSort= [reSort; counter 2 temp(2) temp(3)]; % Psi -> 2
        reSort= [reSort; counter 2 temp']; % Psi -> 2
    end
    
    % CHI1 index-searching:
    % Now depending on the residue: decide what you're gonna do:
    calc_chi1 = 1; % Defaults to calculating (except GLY and ALA)
    resname = pdb.resname(pdb.resseq == resNum,:);
    resname = resname(1,:);
    resname = standardizeProtonatedStateName(resname); % Change protonated state names to standard
    % Standardize resname (Including protonation states)

    if sum(ismember(res_CG,resname,'rows')) == 1
        atoms{4} = selectname(pdb.name, 'CG');
    elseif sum(ismember(res_CG1,resname,'rows')) == 1
        atoms{4} = selectname(pdb.name, 'CG1');
    elseif sum(ismember(res_OG1,resname,'rows')) == 1
        atoms{4} = selectname(pdb.name, 'OG1');
    elseif sum(ismember(res_OG,resname,'rows')) == 1
        atoms{4} = selectname(pdb.name, 'OG');
    elseif sum(ismember(res_SG,resname,'rows')) == 1
        atoms{4} = selectname(pdb.name, 'SG');

    % If it's Ala or Gly just calculate backbone dihedrals
    elseif sum(ismember(res_break,resname,'rows')) == 1
        calc_chi1 = 0;
    end
    
    % reSort section!!!!!
    if calc_chi1
        % Activate an OR to get all the indices of the required atoms
        atoms_all = atoms{1} | atoms{2} | atoms{3} | atoms{4};
        % The dihedral for this residue:
        index_dihedral = index_res & atoms_all; % This will be used in the loop
        % where the dihedral calculations are done
        
        temp = find(index_dihedral);
%         reSort = [reSort; counter 0 temp(2) temp(3)]; % chi1 -> 0
        reSort = [reSort; counter 0 temp']; % chi1 -> 0
    end
    
    % Calculate the higher order chi terms (ch2 - chi5) if the option
    % is given AND the amino acid actually has higher order terms
    res_chi2 = ['ARG '; 'ASN '; 'ASP '; 'GLN '; 'GLU '; 'HIS '; 'ILE '; ...
        'LEU '; 'LYS '; 'MET '; 'PHE '; 'PRO '; 'TRP '; 'TYR '];
    if (higherOrder == 'all' & sum(ismember(res_chi2,resname,'rows')) == 1)

        % CHI2:
        atoms_c2{1} = selectname(pdb.name, 'CA');
        atoms_c2{2} = selectname(pdb.name, 'CB');
        atoms_c2{3} = selectname(pdb.name, 'CG'); % Except for Isoleucine,
        % which will be overwritten later

        % CHI3: defined for ARG, GLN, GLU, LYS, MET
        atoms_c3{1} = selectname(pdb.name, 'CB');
        atoms_c3{2} = selectname(pdb.name, 'CG');
        atoms_c3{3} = selectname(pdb.name, 'CD'); % Except for Methionine,
        % which will be overwritten later
        calc_chi3 = 0; % Becomes 1 if a residue with chi3 is encountered
        calc_chi4 = 0;
        calc_chi5 = 0;
        % Atoms for CHI4 (ARG and LYS) and CHI5 (ARG) will be directly 
        % defined in the loop

        % Here comes the rub PART 2:
        res_CD_OE1 = ['GLN '; 'GLU '];
        res_OD1 = ['ASP' ; 'ASN'];
        res_CD1 = ['LEU' ; 'TRP' ; 'TYR' ; 'PHE'];
        % 
        if sum(ismember(res_CD_OE1,resname,'rows')) == 1
        atoms_c2{4} = selectname(pdb.name, 'CD');
        atoms_c3{4} = selectname(pdb.name, 'OE1');
        calc_chi3 = 1;
        elseif sum(ismember(res_OD1,resname,'rows')) == 1
        atoms_c2{4} = selectname(pdb.name, 'OD1');  
        elseif sum(ismember(res_CD1,resname,'rows')) == 1
        atoms_c2{4} = selectname(pdb.name, 'CD1');    
        elseif resname == 'PRO '
        atoms_c2{4} = selectname(pdb.name, 'CD');    
        elseif resname == 'HIS '
        atoms_c2{4} = selectname(pdb.name, 'ND1'); 
        elseif resname == 'ILE '
        atoms_c2{3} = selectname(pdb.name, 'CG1');  
        atoms_c2{4} = selectname(pdb.name, 'CD', 'CD1'); % Either called CD or CD1
        elseif resname == 'MET '
        atoms_c2{4} = selectname(pdb.name, 'SD'); 
        atoms_c3{3} = selectname(pdb.name, 'SD');
        atoms_c3{4} = selectname(pdb.name, 'CE');
        calc_chi3 = 1;
        elseif resname == 'LYS '
        atoms_c2{4} = selectname(pdb.name, 'CD');
        atoms_c3{4} = selectname(pdb.name, 'CE');
        atoms_c4{1} = selectname(pdb.name, 'CG');
        atoms_c4{2} = selectname(pdb.name, 'CD');
        atoms_c4{3} = selectname(pdb.name, 'CE');
        atoms_c4{4} = selectname(pdb.name, 'NZ');
        calc_chi3 = 1;
        calc_chi4 = 1;
        elseif resname == 'ARG '
        atoms_c2{4} = selectname(pdb.name, 'CD');
        atoms_c3{4} = selectname(pdb.name, 'NE');
        atoms_c4{1} = selectname(pdb.name, 'CG');
        atoms_c4{2} = selectname(pdb.name, 'CD');
        atoms_c4{3} = selectname(pdb.name, 'NE');
        atoms_c4{4} = selectname(pdb.name, 'CZ');
        atoms_c5{1} = selectname(pdb.name, 'CD');
        atoms_c5{2} = selectname(pdb.name, 'NE');
        atoms_c5{3} = selectname(pdb.name, 'CZ');
        atoms_c5{4} = selectname(pdb.name, 'NH1');
        calc_chi3 = 1;
        calc_chi4 = 1;    
        calc_chi5 = 1; 
        end

        % Do we calculate chiN? calc_chi2 = 1 by default
        % atoms_higher indices will be used in the calculation in the
        % numRuns for loop
        calc_chi = [calc_chi1 1 calc_chi3 calc_chi4 calc_chi5];
        Nchis = 2 + calc_chi3 + calc_chi4 + calc_chi5;
        atoms_higher = cell(Nchis-1,1);
        atoms_higher{1} = atoms_c2{1} | atoms_c2{2} | atoms_c2{3} | atoms_c2{4};
        
        % reSort!!!!!
        index_dihedral_higher = index_res & atoms_higher{1};
        temp = find(index_dihedral_higher);
%         reSort = [reSort; counter 0 temp(2) temp(3)]; % chi2 -> 0
        reSort = [reSort; counter 0 temp']; % chi2 -> 0
        
        if calc_chi3 == 1
        atoms_higher{2} = atoms_c3{1} | atoms_c3{2} | atoms_c3{3} | atoms_c3{4};
        
        % reSort!!!!!
        index_dihedral_higher = index_res & atoms_higher{2};
        temp = find(index_dihedral_higher);
%         reSort = [reSort; counter 0 temp(2) temp(3)]; % chi3 -> 0
        reSort = [reSort; counter 0 temp']; % chi3 -> 0
        end
        if calc_chi4 == 1
        atoms_higher{3} = atoms_c4{1} | atoms_c4{2} | atoms_c4{3} | atoms_c4{4};
        
        % reSort!!!!!
        index_dihedral_higher = index_res & atoms_higher{3};
        temp = find(index_dihedral_higher);
%         reSort = [reSort; counter 0 temp(2) temp(3)]; % chi4 -> 0
        reSort = [reSort; counter 0 temp']; % chi4 -> 0
        end
        if calc_chi5 == 1
        atoms_higher{4} = atoms_c5{1} | atoms_c5{2} | atoms_c5{3} | atoms_c5{4};
        
        % reSort!!!!!
        index_dihedral_higher = index_res & atoms_higher{4};
        temp = find(index_dihedral_higher);
%         reSort = [reSort; counter 0 temp(2) temp(3)]; % chi5 -> 0
        reSort = [reSort; counter 0 temp']; % chi5 -> 0
        end
    else
        calc_chi = [calc_chi1 0 0 0 0]; % We're not calculating any higher order chis
    end   
    dihIndex(resNum,:) = [calc_bb calc_chi];
    % Now loop over the runs and execute the calculations
    for runi = 1:numRuns 
        % Length of each element in dihedral: 2 + sum(calc_chi);
        dihedrals{counter, runi} = zeros(nframes,(2+sum(calc_chi)));
        resname_cell{counter, runi} = resname;
        
        % PHI-PSI calculations:
        if calc_bb(1) % calculate phi
            % Extract trajectory of dihedral that residue ONLY
            if numRuns == 1 %traj is a matrix
                traj_dihedral = traj(:,to3(index_phi));
            else            % traj is a cell
                traj_dihedral = traj{runi}(:,to3(index_phi));
            end
%             phi = calcdihedral(traj_dihedral).*180./pi;
            phi = calcdihedral(traj_dihedral);
        else % phi not defined, set it to NaN
            phi = nan(nframes,1);
        end
        dihedrals{counter, runi}(:,1) = phi;
        
        if calc_bb(2) % calculate psi
            % Extract trajectory of dihedral that residue ONLY
            if numRuns == 1 %traj is a matrix
                traj_dihedral = traj(:,to3(index_psi));
            else            % traj is a cell
                traj_dihedral = traj{runi}(:,to3(index_psi));
            end
%             psi = calcdihedral(traj_dihedral).*180./pi;
            psi = calcdihedral(traj_dihedral);
            else % psi not defined, set it to NaN
            psi = nan(nframes,1);
        end
        dihedrals{counter, runi}(:,2) = psi;
        
        % CHI1 calculations   
        if calc_chi1 % Not ALA or GLY
            % Extract trajectory of dihedral that residue ONLY
            if numRuns == 1 %traj is a matrix
                traj_dihedral = traj(:,to3(index_dihedral));
            else            % traj is a cell
                traj_dihedral = traj{runi}(:,to3(index_dihedral));
            end
%             chi1 = calcdihedral(traj_dihedral).*180./pi;
            chi1 = calcdihedral(traj_dihedral);
            dihedrals{counter, runi}(:,3) = chi1;
        end
        
        % Execute higher order calculations?
        if (higherOrder == 'all' & sum(ismember(res_chi2,resname,'rows')) == 1)
            % Calculate chi3/4/5 if calc_chiN = 1
            for chiN=2:5
            if calc_chi(chiN) == 1 
                index_dihedral_higher = index_res & atoms_higher{chiN-1};
                if numRuns == 1 %traj is a matrix
                    traj_dihedral = traj(:,to3(index_dihedral_higher));
                else            % traj is a cell
                    traj_dihedral = traj{runi}(:,to3(index_dihedral_higher));
                end
%                 chiHigher = calcdihedral(traj_dihedral).*180./pi;
                chiHigher = calcdihedral(traj_dihedral);
                dihedrals{counter, runi}(:,chiN + 2) = chiHigher;
            end
            end
        end
        if transform~=1 % transform from [-180,180] to [0,360]
            for dih = 1:(2+sum(calc_chi))
                for k=1:nframes
                if dihedrals{counter, runi}(k,dih) < 0
%                     dihedrals{counter, runi}(k,dih) = dihedrals{counter, runi}(k,dih) + 360;
                    dihedrals{counter, runi}(k,dih) = dihedrals{counter, runi}(k,dih) + 2*pi;
                end   
                end
            end
        end
    end
    counter = counter + 1;
end
end

