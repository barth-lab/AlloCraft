function writepdbndx(filename, pdb, ndx, renum, trj, format_type, bfactor)
%% writepdb
% write Protein Data Bank (PDB) file
%
%% Syntax
%# writepdb(filename, pdb);
%# writepdb(filename, pdb, ndx, renum, crd);
%# writepdb(filename, pdb, ndx, renum, trj);
%# writepdb(filename, pdb, ndx, renum, crd, format_type);
%# writepdb(filename, pdb, ndx, renum, [], format_type);
%# writepdb(filename, pdb, ndx, renum, [], format_type, bfactor);
%
%% Description
% This code writes only just ATOM or HETATM records. 
% Note that icode fields (code for insetion of residues, see
% "References" below) are ignored in this routine. 
%
% * filename  - input filename of PDB [char]
% * pdb       - structure data
%          record: 'ATOM  ' or 'HETATM' [natomx6 char]
%          serial: Atom serial number. [natomx1 double]
%            name: Atom name. [natomx4 char]
%          altloc: Alternate location indicator. [natomx1 char]
%         resname: Residue name. [natomx3 char]
%         chainid: Chain identifier. [natomx1 char]
%          resseq: Residue sequence number. [natomx1 double]
%           icode: Code for insertion of residues. [natomx1 char] ()
%             xyz: Cartesian coordinate of atom in Angstrom [natomx3 double]
%       occupancy: Occupancy. [natomx1 double]
%      tempfactor: Temperature factor. [natomx1 double]
%         element: Element symbol, right-justified. [natomx2 char]
%          charge: Charge on the atom. [natomx2 char]
% * ndx         -  index of the atoms that you want to save to the pdb file
%                  [1 x natom] logical with 1 for the atoms you want to
%                  save or simply the number of atoms in the serial.
% * renum       -  1 if you want the residue and atom numbers to start from
%                  1. Defaults to 0.
% * crd, trj    -  coordinates, or trajecotry of the molecule
%                  if given, the coordiates of pdb is replaced with
%                  this data. Trajectory is written as mutiple
%                  models in pdb file. 
%                  [1 x 3natom double] or [nframe x 3natom double]
% * format_type -  format type [chars. only 'vmd' can be available,
%                  otherwise default(standard?) format is used]
%
%% Example
%# pdb = readpsf('jac.pdb');
%# pdb.xyz(:, 1) = pdb.xyz(:, 1) + 2.5; %translate in x-axis by 2.5 Angstrom
%# writepdb('jac_x.pdb', pdb);
%# writepdb('jac_x_vmd.pdb', pdb, [], 'vmd');
%
%% See also
% readpdb
% 
%% References
% http://www.wwpdb.org/documentation/format33/sect9.html
% ATOM Record Format
% COLUMNS        DATA  TYPE    FIELD        DEFINITION
% -------------------------------------------------------------------------------------
%  1 -  6        Record name   "ATOM  " or "HETATM"
%  7 - 11        Integer       serial       Atom serial number.
% 13 - 16        Atom          name         Atom name.
% 17             Character     altLoc       Alternate location indicator.
% 18 - 20        Residue name  resName      Residue name.
% 22             Character     chainID      Chain identifier.
% 23 - 26        Integer       resSeq       Residue sequence number.
% 27             AChar         iCode        Code for insertion of residues.
% 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
% 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
% 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
% 55 - 60        Real(6.2)     occupancy    Occupancy.
% 61 - 66        Real(6.2)     tempFactor   Temperature factor.
% 77 - 78        LString(2)    element      Element symbol, right-justified.
% 79 - 80        LString(2)    charge       Charge on the atom.

%% check existing file
if exist(filename, 'file')
  filename_old = sprintf('%s.old', filename);
  display(sprintf('existing file %s is moved to %s', filename, filename_old));
  movefile(filename, filename_old);
end

%% preparation
natom = size(pdb.record, 1);

if (nargin < 3) || isempty(ndx) % initialize atom selection to all
    ndx = pdb.serial';
else
    if  islogical(ndx)
        ndx = pdb.serial(ndx)';
    end
end

if (nargin < 4) % Renumber or not?
    renum = 0; % Default is do not renumber
end

if (nargin < 5) || (numel(trj) == 0)
  trj = pdb.xyz';
  trj = trj(:)';
end
nframe = size(trj, 1);
  
if nargin < 6
  format_type = 'default';
end

is_bfactor = false;
if (nargin > 6) && (numel(trj) ~= 0)
  is_bfactor = true;
end

%% open file
fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open file.');
cleaner = onCleanup(@() fclose(fid));

if strncmpi(format_type, 'namd', numel('namd'))
  fprintf(fid, 'CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n');
end

%% write file
ndxRes = pdb.resseq(ndx); % Sequence of original pdb
for iframe = 1:nframe

  if(nframe > 1)
    fprintf(fid, 'MODEL %8d\n', iframe);
  end
  
  resCounter = 1; % Residue counter for renumbering
  i = 1; % counter loops
  for iatom = ndx
      if renum == 1
            if i==1
                thisRes = 1;
            elseif ndxRes(i) == ndxRes(i-1) % Are we in the same residue still?
                thisRes = resCounter;
            else % we moved to next residue, increase counter!
                resCounter = resCounter + 1;
                thisRes = resCounter;
            end
            thisAtom = i; % Renumber atoms starting from 1 too!
            i = i + 1;
      else
          thisRes = pdb.resseq(iatom);
          thisAtom = pdb.serial(iatom);
      end
    if strncmpi(format_type, 'vmd', numel('vmd'))
      % VMD format
      fprintf(fid, '%6s', pdb.record(iatom, :));
      fprintf(fid, '%5d', mod(thisAtom, 100000));
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%4s', pdb.name(iatom, :));
      fprintf(fid, '%1s', pdb.altloc(iatom, :));
      %fprintf(fid, '%3s', pdb.resname(iatom, :));
      %fprintf(fid, '%1s', ' ');
      fprintf(fid, '%4s', pdb.resname(iatom, :));
      fprintf(fid, '%1s', pdb.chainid(iatom, :));
      fprintf(fid, '%4d', mod(thisRes, 10000));
      %fprintf(fid, '%1s', pdb.icode(iatom, :));
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 1));
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 2));
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 3));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+1));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+2));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+3));
      fprintf(fid, '%6.2f', pdb.occupancy(iatom));
      if is_bfactor
        fprintf(fid, '%6.2f', bfactor(iframe, iatom));
      else
        fprintf(fid, '%6.2f', pdb.tempfactor(iatom));
      end
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%2s', pdb.element(iatom, :));
      fprintf(fid, '%2s', pdb.charge(iatom, :));
      fprintf(fid, '\n');
      
    elseif strncmpi(format_type, 'namd', numel('namd'))
      % NAMD format
      fprintf(fid, '%6s', pdb.record(iatom, :));
      if thisAtom < 100000
        fprintf(fid, '%5d', thisAtom);
      else
        fprintf(fid, '%5s', lower(dec2hex(thisAtom)));
      end
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%4s', pdb.name(iatom, :));
      fprintf(fid, '%1s', pdb.altloc(iatom, :));
      %fprintf(fid, '%3s', pdb.resname(iatom, :));
      %fprintf(fid, '%1s', ' ');
      fprintf(fid, '%4s', pdb.resname(iatom, :));
      fprintf(fid, '%1s', pdb.chainid(iatom, :));
      fprintf(fid, '%4d', mod(thisRes, 10000));
      %fprintf(fid, '%1s', pdb.icode(iatom, :));
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 1));
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 2));
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 3));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+1));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+2));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+3));
      fprintf(fid, '%6.2f', pdb.occupancy(iatom));
      if is_bfactor
        fprintf(fid, '%6.2f', bfactor(iframe, iatom));
      else
        fprintf(fid, '%6.2f', pdb.tempfactor(iatom));
      end
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      %fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', pdb.chainid(iatom, :));
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%2s', pdb.element(iatom, :));
      %fprintf(fid, '%2s', pdb.charge(iatom, :));
      fprintf(fid, '\n');
    
    else
      % PDB format with some extentions for large digits of atom_serial and resseq
      fprintf(fid, '%6s', pdb.record(iatom, :));
      if thisAtom < 100000
        fprintf(fid, '%5d', thisAtom);
        fprintf(fid, '%1s', ' ');
      else
        fprintf(fid, '%6d', thisAtom);
      end
      fprintf(fid, '%4s', pdb.name(iatom, :));
      fprintf(fid, '%1s', pdb.altloc(iatom, :));
      %fprintf(fid, '%3s', pdb.resname(iatom, :));
      %fprintf(fid, '%1s', ' ');
      fprintf(fid, '%4s', pdb.resname(iatom, :));
      fprintf(fid, '%1s', pdb.chainid(iatom, :));
      if thisRes < 10000
        fprintf(fid, '%4d', thisRes);
        %fprintf(fid, '%1s', pdb.icode(iatom, :));
        fprintf(fid, '%1s', ' ');
        fprintf(fid, '%1s', ' ');
      elseif thisRes < 100000
        fprintf(fid, '%5d', thisRes);
        fprintf(fid, '%1s', ' ');
      else
        fprintf(fid, '%6d', thisRes);
      end
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 1));
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 2));
      %fprintf(fid, '%8.3f', pdb.xyz(iatom, 3));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+1));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+2));
      fprintf(fid, '%8.3f', trj(iframe, 3*(iatom-1)+3));
      fprintf(fid, '%6.2f', pdb.occupancy(iatom));
      if is_bfactor
        fprintf(fid, '%6.2f', bfactor(iframe, iatom));
      else
        fprintf(fid, '%6.2f', pdb.tempfactor(iatom));
      end
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%1s', ' ');
      fprintf(fid, '%2s', pdb.element(iatom, :));
      fprintf(fid, '%2s', pdb.charge(iatom, :));
      fprintf(fid, '\n');
    end
  end
  
  fprintf(fid, 'TER\n');
  fprintf(fid, 'ENDMDL\n');
  
end

fprintf(fid, 'END\n');

