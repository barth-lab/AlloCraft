
function xyz = readdcdmat(filename, ind)

% xyz = readdcd(filename, indices)
% reads an dcd and puts the x,y,z coordinates corresponding to indices 
% in the rows of x,y,z
% If ind is not specified, it reads all the atoms

h = read_dcdheader(filename);
nsets = h.NSET;
natoms = h.N;

% Set the default value for ind
if nargin<2
  ind = 1:natoms;
end
numind = length(ind);

x = zeros(natoms,1);
y = zeros(natoms,1);
z = zeros(natoms,1);

if nsets == 0
  xyz = zeros(1,3*numind);
  nsets = 99999;
else
  xyz = zeros(nsets, 3*numind);
end

for i=1:nsets
  pos = ftell(h.fid);
  if pos == h.endoffile 
    break;
  end
  [x,y,z] = read_dcdstep(h);
  xyz(i,1:3:3*numind) = x(ind)';
  xyz(i,2:3:3*numind) = y(ind)';
  xyz(i,3:3:3*numind) = z(ind)';
end

close_dcd(h);

