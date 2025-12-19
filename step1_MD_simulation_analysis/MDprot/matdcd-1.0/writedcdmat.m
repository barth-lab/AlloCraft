
function rc = writedcd(filename, x, y, z)

%
% rc = writedcd(filename, x, y, z)
% x, y, z are of size [natoms nsets]
%

[natoms, nsets] = size(x);

h = write_dcdheader(filename, nsets, natoms);
for i=1:nsets
  rc = write_dcdstep(h, x(:,i), y(:,i), z(:,i));
end

close_dcd(h);
 
