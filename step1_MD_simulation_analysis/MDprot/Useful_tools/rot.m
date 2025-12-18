function [R] = rot(u, ang, if_rad)
%rot Rotation matrix around a 3-dimensional vector u by an angle "ang"
%   R = rot(u, ang) creates a 3-by-3 matrix for rotating a 3-by-1 vector or 
%   3-by-N matrix of vectors around the vector u by ang degrees. When acting
%   on a matrix, each column of the matrix represents a different vector. 
%   For the rotation matrix R and vector v, the rotated vector is given by R*v.
%  Usage:
%  [R] = rot(u, ang)
%  [R] = rot(u, ang, if_rad)


if exist('if_rad', 'var')
      R = [cos(ang) + u(1)^2*(1 - cos(ang)) u(1)*u(2)*(1-cos(ang))-u(3)*sin(ang) u(1)*u(3)*(1-cos(ang))+u(2)*sin(ang) ;...
   u(1)*u(2)*(1-cos(ang))+u(3)*sin(ang) cos(ang) + u(2)^2*(1 - cos(ang)) u(2)*u(3)*(1-cos(ang))-u(1)*sin(ang); ...
   u(1)*u(3)*(1-cos(ang))-u(2)*sin(ang) u(2)*u(3)*(1-cos(ang))+u(1)*sin(ang) cos(ang) + u(3)^2*(1 - cos(ang))]; 
else 
   R = [cosd(ang) + u(1)^2*(1 - cosd(ang)) u(1)*u(2)*(1-cosd(ang))-u(3)*sind(ang) u(1)*u(3)*(1-cosd(ang))+u(2)*sind(ang) ;...
   u(1)*u(2)*(1-cosd(ang))+u(3)*sind(ang) cosd(ang) + u(2)^2*(1 - cosd(ang)) u(2)*u(3)*(1-cosd(ang))-u(1)*sind(ang); ...
   u(1)*u(3)*(1-cosd(ang))-u(2)*sind(ang) u(2)*u(3)*(1-cosd(ang))+u(1)*sind(ang) cosd(ang) + u(3)^2*(1 - cosd(ang))]; 
end

end
