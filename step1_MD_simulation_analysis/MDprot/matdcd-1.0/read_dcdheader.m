
function h = read_dcdheader(filename)

% new_handle = read_dcdheader(filename)

fid = fopen(filename, 'r', 'b');
i = fread(fid, 1, 'int32');
if i ~= 84              % wrong endianism; try again
  fclose(fid);
  fid = fopen(filename, 'r', 'l');
  i = fread(fid, 1, 'int32');
  if i ~= 84
    fclose(fid);
  end
end

h.fid = fid;

%
% Figure out how big the file is
%
fseek(fid,0,1);
h.endoffile = ftell(fid);
fseek(fid,4,-1);

% Read CORD 
s = fread(fid, 4, 'char');

% Read NSET
h.NSET = fread(fid, 1, 'int32');

% read ISTART
h.ISTART = fread(fid, 1, 'int32');

% read NSAVC
h.NSAVC = fread(fid, 1, 'int32');

% Reposition to 40 from beginning; read namnf, number of free atoms
fseek(fid, 40, -1);
h.NAMNF = fread(fid, 1, 'int32');

% Figure out if we're CHARMm or not
fseek(fid, 84, -1);
i = fread(fid, 1, 'int32');
if i == 0
  h.charmm = 0;
else
  h.charmm = 1;
  h.charmm_extrablock = 0;
  h.charmm_4dims = 0;
  % check for extra block
  fseek(fid, 48, -1);
  i = fread(fid, 1, 'int32');
  if i == 1
    h.charmm_extrablock = 1;
  end
  % check for 4 dims
  fseek(fid, 52, -1);
  i = fread(fid, 1, 'int32');
  if i == 1
    h.charmm_4dims = 1;
  end
end



% Read the timestep, DELTA.  Read a float if charmm, otherwise double
fseek(fid, 44, -1);
if h.charmm == 0
  h.DELTA = fread(fid, 1, 'float64');
else
  h.DELTA = fread(fid, 1, 'float32');
end

% Get the size of the next block, and skip it
% This is the title
fseek(fid, 92, -1);
newsize = fread(fid, 1, 'int32');
numlines = fread(fid, 1, 'int32');
fseek(fid, numlines*80, 0);
newsize = fread(fid, 1, 'int32');

% Read in a 4, then the number of atoms, then another 4
i = fread(fid, 1, 'int32');
h.N = fread(fid, 1, 'int32');
i = fread(fid, 1, 'int32');


% stuff with freeindexes.  Just smile and nod.
if h.NAMNF ~= 0
  fsize = fread(fid, 1, 'int32');  % should be N-NAMNF*4
  h.FREEINDEXES = fread(fid, h.N - h.NAMNF, 'int32');
  fsize = fread(fid, 1, 'int32');  % should be N-NAMNF*4
end 

