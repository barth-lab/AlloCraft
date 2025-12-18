
function [x,y,z] = read_dcdstep(h)

%
% [x,y,z] = read_dcdstep(handle)
%

% If this is a CHARMm file and contains an extra data block, we must skip it
if h.charmm & h.charmm_extrablock
  blocksize = fread(h.fid, 1, 'int32');
  fseek(h.fid, blocksize, 0);
  blocksize = fread(h.fid, 1, 'int32');
end

if h.NAMNF == 0

  % Get x coordinates 
  blocksize = fread(h.fid, 1, 'int32');
  x = fread(h.fid, blocksize/4, 'float32');
  blocksize = fread(h.fid, 1, 'int32');

  % Get y coordinates 
  blocksize = fread(h.fid, 1, 'int32');
  y = fread(h.fid, blocksize/4, 'float32');
  blocksize = fread(h.fid, 1, 'int32');

  % Get z coordinates 
  blocksize = fread(h.fid, 1, 'int32');
  z = fread(h.fid, blocksize/4, 'float32');
  blocksize = fread(h.fid, 1, 'int32');

else    
 
  % this is not implemented in the VMD code I copied from  

end  
 
% Skip the 4th dimension, if present
if h.charmm & h.charmm_4dims
  blocksize = fread(h.fid, 1, 'int32');
  fseek(h.fid, blocksize, 0);
  blocksize = fread(h.fid, 1, 'int32');
end

