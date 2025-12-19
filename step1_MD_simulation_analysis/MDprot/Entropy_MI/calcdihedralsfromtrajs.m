function [dihedrals, resname_cell] = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns,higherOrder,transform)
%% calcdihedralsfromtrajs: 
% This function takes a set of trajectories and calculates the rotamer dihedral (Chi) of  chosen residues.
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
% Dihedrals based on: http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
%
%% Usage:
% dihedrals = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns);
% dihedrals = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns, 'all',transform);
% dihedrals = calcdihedralsfromtrajs(pdb,traj,rotamers);
% [dihedrals, resname_cell] = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns);
% [dihedrals, resname_cell] = calcdihedralsfromtrajs(pdb,traj,rotamers,numRuns, 'all',transform);
% [dihedrals, resname_cell] = calcdihedralsfromtrajs(pdb,traj,rotamers);
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
% * rotamers is the list of residues to calculate chiN for.
%
% * numRuns number of cells to consider in traj (if traj is a cell
% structure).
%
% * dihedrals is a length(rotamers) x numRuns cell structure containing the
% dihedral angles, where each cell has Nframes x 1 angles. It outputs a NaN 
% for Alanines and Glycines. If the option 'all' is given, dihedrals will
% also contain higher order dihedrals for residues that have them, and each
% cell will have a size of Nframes x chiN, where chiN is the number of
% available dihedrals for this residue
%
% * resname_cell is length(rotamers) x numRuns cell structure containing the
% names of the chosen residues, useful when dealing with homologues that
% may have different residues in the same spot.
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
counter = 1; % Counter for number of residues
    
for resNum=rotamers    % Choose a residue from the list of rotamers
    index_res = selectid(pdb.resseq, resNum);

    % CHI1 index-searching:
    % Now depending on the residue: decide what you're gonna do:
    resname = pdb.resname(pdb.resseq == resNum,:);
    resname = resname(1,:);
    resname = standardizeProtonatedStateName(resname); % Change protonated state names to standard

    % Choose the atoms for the dihedral calculation:
    % The first three atoms are always the same
    atoms{1} = selectname(pdb.name, 'N');
    atoms{2} = selectname(pdb.name, 'CA');
    atoms{3} = selectname(pdb.name, 'CB');

    % Here comes the rub: (GLY and ALA don't have a chi angle, and 
    % their values will be replaced by a NaN)
    res_CG = ['ARG '; 'ASN '; 'ASP '; 'GLN '; 'GLU '; 'HIS '; 'HSD ';'HSE '; ...
        'LEU '; 'LYS '; 'MET '; 'PHE '; 'PRO '; 'TRP '; 'TYR '];
    res_CG1 = ['ILE '; 'VAL '];
    res_OG1 = ['THR '];
    res_OG = ['SER '];
    res_SG = ['CYS '];
    res_break = ['ALA ';'GLY '];
    % 
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

    % If it's Ala or Gly just give NaN to ALL RUNS, end this iteratrion and continue
    elseif sum(ismember(res_break,resname,'rows')) == 1
        for runi = 1:numRuns
            dihedrals{counter, runi} = NaN;  
            resname_cell{counter, runi} = resname;
        end
        counter = counter + 1;
        continue
    end
    % Activate an OR to get all the indices of the required atoms
    atoms_all = atoms{1} | atoms{2} | atoms{3} | atoms{4};
    % The dihedral for this residue:
    index_dihedral = index_res & atoms_all; % This will be used in the loop
    % where the dihedral calculations are done

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
        atoms_c2{4} = selectname(pdb.name, 'CD');
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
        calc_chi = [1 calc_chi3 calc_chi4 calc_chi5];
        Nchis = 2 + calc_chi3 + calc_chi4 + calc_chi5;
        atoms_higher = cell(Nchis-1,1);
        atoms_higher{1} = atoms_c2{1} | atoms_c2{2} | atoms_c2{3} | atoms_c2{4};
        if calc_chi3 == 1
        atoms_higher{2} = atoms_c3{1} | atoms_c3{2} | atoms_c3{3} | atoms_c3{4};
        end
        if calc_chi4 == 1
        atoms_higher{3} = atoms_c4{1} | atoms_c4{2} | atoms_c4{3} | atoms_c4{4};
        end
        if calc_chi5 == 1
        atoms_higher{4} = atoms_c5{1} | atoms_c5{2} | atoms_c5{3} | atoms_c5{4};
        end
    end   
    
    % Now loop over the runs and execute the calculations
    for runi = 1:numRuns 
        % CHI1 calculations    
        % Extract trajectory of dihedral that residue ONLY
        if numRuns == 1 %traj is a matrix
            traj_dihedral = traj(:,to3(index_dihedral));
        else            % traj is a cell
            traj_dihedral = traj{runi}(:,to3(index_dihedral));
        end
        chi1 = calcdihedral(traj_dihedral).*180./pi;
        
        if transform~=1 % transform from [-180,180] to [0,360]
            for k=1:length(chi1)
            if chi1(k) < 0
                chi1(k) = chi1(k) +360;
            end 
            end
        end
        
        dihedrals{counter, runi} = chi1;
        resname_cell{counter, runi} = resname;
        
        % Execute higher order calculations?
        if (higherOrder == 'all' & sum(ismember(res_chi2,resname,'rows')) == 1)
            % Expand the dihedral cell so it can hold as many chis as needed
            dihedrals{counter, runi} = zeros(length(chi1),Nchis);
            % Return chi1 to dihedrals variable:
            dihedrals{counter, runi}(:,1) = chi1;
            % Calculate chi3/4/5 if calc_chiN = 1
            for chiN=2:5
            if calc_chi(chiN-1) == 1 
                index_dihedral_higher = index_res & atoms_higher{chiN-1};
                if numRuns == 1 %traj is a matrix
                    traj_dihedral = traj(:,to3(index_dihedral_higher));
                else            % traj is a cell
                    traj_dihedral = traj{runi}(:,to3(index_dihedral_higher));
                end
                chiHigher = calcdihedral(traj_dihedral).*180./pi;
                
                if transform~=1 % transform from [-180,180] to [0,360]
                    for k=1:length(chiHigher)
                    if chiHigher(k) < 0
                        chiHigher(k) = chiHigher(k) +360;
                    end   
                    end
                end
                
                dihedrals{counter, runi}(:,chiN) = chiHigher;
            end
            end
        end
    end
    counter = counter + 1;
end
end

