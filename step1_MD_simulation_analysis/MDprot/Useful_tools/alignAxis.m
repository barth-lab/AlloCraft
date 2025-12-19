function [crd_rot, rot_mat] = alignAxis(pdb, res, u, name_pdb)
% alignAxis Aligns the structure from pdb by aligning the axis defined by
% the two residues in res to the vector u
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
%
% Syntax:
% [crd_rot] = alignAxis(pdb, res, u)
% [crd_rot, rot_mat] = alignAxis(pdb, res, u)
% [crd_rot, rot_mat] = alignAxis(pdb, res, u, name_pdb)
% [crd_rot, rot_mat] = alignAxis(crd, atom_pair, u)
%
% Description:
% * pdb is the pdb structure obtained by pdb = readpdb('pdb.pdb'). Note
% that other structure files (.gro for example) will not work, as the way
% mdtoolbox names the structure elements is different.
% * res is a pair of RESIDUES that form a vector v, which will be aligned
% to u. [1 x 2] vector containing two residue numbers (if the first input 
% is a pdb structure)
% * crd is a matrix of coordinates in the form of [1 x 3*natom], which can be
% obtained by [~,crd] = readpdb('pdb.pdb').
% * atom_pair is a pair of ATOMS that form a vector v, which will be aligned to u, if
% the first input is a [1 x 3*natom] matrix of coordinates.
% * u is the target vector to which v will be aligned to [1 x 3] vector.
% * name_pdb is an optional parameter, if given, the function will write
% out a PDB file with the name "name_pdb". 
% * crd_rot are the coordinates of the aligned structure. [nframe x
% 3*natom].
% * rot_mat is the rotation matrix used [3 x 3].
%
% Examples:
% pdb = readpdb('pdb.pdb')
% residue_pair = [1 50]
% [crd_rot, rot_mat] = alignAxis(pdb, residue_pair, [0 0 1], 'pdb_aligned.pdb')
%
% [~,crd] = readpdb('pdb.pdb')
% atom_pair = [1 2300]
% [crd_rot, rot_mat] = alignAxis(crd, atom_pair, [0 0 1])

if isnumeric(pdb) == 1 % Input is coordinates, not a pdb file
    Natoms = length(pdb)/3;
    crd = decenter(pdb); % Decenter the coordinates
    res1 = res(1);
    res2 = res(2);
    % Extract coordinates
    crd_CA1 = crd(:, 3*res1-2:3*res1);
    crd_CA2 = crd(:, 3*res2-2:3*res2);
else
    Natoms =  length(pdb.serial);
    % Extract and format cooridnates from PDB
    crd = pdb.xyz;
    crd_new = zeros(1,3*Natoms);
    for i=1:Natoms
       crd_new(3*i-2) =  crd(i,1);
       crd_new(3*i-1) = crd(i,2);
       crd_new(3*i) = crd(i,3);
    end
    crd = crd_new;
    crd = decenter(crd); % Decenter the coordinates

    res1 = res(1);
    res2 = res(2);
    % Grab the indices of the residues
    index1 = selectname(pdb.name, 'CA') & selectid(pdb.resseq,res1);
    index2 = selectname(pdb.name, 'CA') & selectid(pdb.resseq,res2);

    % Extract coordinates
    crd_CA1 = crd(:, to3(index1));
    crd_CA2 = crd(:, to3(index2));    
end

v = (crd_CA2 - crd_CA1)/norm(crd_CA2 - crd_CA1); % unit vector defining the line
u = u/norm(u); % Make sure u is a unit vector too

% Find the angle between v and u:
ang = acosd(dot(u,v));

% Find a vector perpendicular to both v and the u (axis you want to align
% to): (to use as an axis of rotation)
w = cross(v, u);
w = w/norm(w); % Normalize w

% Calculate rotation matrix around w by ang
rot_mat = rot(w,ang);

% Apply the rot mat on the coordinates
crd_rot = zeros(1,length(crd));
for atom=1:Natoms
   crd_rot(3*atom-2:3*atom)=rot_mat*(crd(3*atom-2:3*atom)');
end

% Write the PDB (if the user asked for it)
if exist('name_pdb', 'var') & isnumeric(pdb) == 0
writepdb(name_pdb,pdb,crd_rot);
end

end