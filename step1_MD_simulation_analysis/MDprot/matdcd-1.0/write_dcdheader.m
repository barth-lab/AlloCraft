
function h = write_dcdheader(filename, nsets, natoms )

%
% handle = write_dcdheader(filename, nsets, natoms ) 
%

fid = fopen(filename, 'w');
h.fid = fid;

%
% write an 84 byte block header
%

fwrite(fid, 84, 'int32');             
fwrite(fid, 'CORD', 'uchar');         %  4
fwrite(fid, nsets, 'int32');          %  4
fwrite(fid, 0, 'int32');              %  4 - ISTART
fwrite(fid, 1, 'int32');              %  4 - NSAVC
fwrite(fid, zeros(6,1), 'int32');     % 24
fwrite(fid, 1, 'float64');            %  8 - DELTA
fwrite(fid, zeros(9,1), 'int32');   % 36
fwrite(fid, 84, 'int32');

%
% Now write a 164 byte block for the title
%

fwrite(fid, 164, 'int32'); 
fwrite(fid, 2, 'int32');     % number of 80-character lines in the title
count = fprintf(fid, 'REMARKS FILENAME=%s CREATED BY MATLAB',filename);
for i = count+1:80
  fwrite(fid, ' ', 'uchar');
end
count = fprintf(fid, 'REMARKS DATE: %s ', datestr(now));
for i = count+1:80
  fwrite(fid, ' ', 'uchar');
end
fwrite(fid, 164, 'int32');

%
% write the number of atoms
%

fwrite(fid, 4, 'int32');
fwrite(fid, natoms, 'int32');
fwrite(fid, 4, 'int32');

rc = 0;


