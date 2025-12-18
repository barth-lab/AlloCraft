
function rc = write_dcdstep(h, x, y, z)

%
% rc = write_dcdstep, handle, x, y, z);
%

fwrite(h.fid, 4*length(x), 'int32');
fwrite(h.fid, x, 'float32');
fwrite(h.fid, 4*length(x), 'int32');

fwrite(h.fid, 4*length(x), 'int32');
fwrite(h.fid, y, 'float32');
fwrite(h.fid, 4*length(x), 'int32');

fwrite(h.fid, 4*length(x), 'int32');
fwrite(h.fid, z, 'float32');
fwrite(h.fid, 4*length(x), 'int32');

rc = 0;

