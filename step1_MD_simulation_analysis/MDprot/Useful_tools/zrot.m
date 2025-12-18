function [R] = zrot(ang,if_rad)
%zrot Rotation matrix around x
%   R = zrot(ang) creates a 3-by-3 matrix for rotating a 3-by-1 vector or 
%   3-by-N matrix of vectors around the x-axis by ang degrees. When acting
%   on a matrix, each column of the matrix represents a different vector. 
%   For the rotation matrix R and vector v, the rotated vector is given by R*v.


if exist('if_rad', 'var')
   R = [cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1]; 
else 
   R = [cosd(ang) -sind(ang) 0; sind(ang) cosd(ang) 0; 0 0 1]; 
end

end
