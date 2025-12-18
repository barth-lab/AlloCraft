function [distMat] = calcDistMatPBC(X,L)
%calcDistMatPBC Calculate distance matrix of 3D input data taking into
% consideration periodic boundary conditions with box size L
% Having specific case for 2D/3D makes the code ~30 times faster than
% nesting a 3D loop that accounts for variable dimensions
%% Usage:
%   [distMat] = calcDistMatPBC(X,L)
% 
%% Description
% * X is the input data with dimensions: npts x 3 
% * L are the dimensions of the periodic box. Scalar input defines a square
% or vector input with [LX LY LZ]
%
%% See also
% calcDistMatPBC2D

[npts,~] = size(X);
if length(L) == 1
    LX = L;
    LY = L;
    LZ = L;
else
    assert(length(L) == 3, 'L dimensions must either be 1 or 3');
end
distMat = zeros(npts);
Xdat = X(:,1);
Ydat = X(:,2);
Zdat = X(:,3);

for pt1 = 1:npts-1
    for pt2 = pt1+1:npts

        dx = abs(Xdat(pt2) - Xdat(pt1));
        if dx > LX/2
            dx = LX - dx;
        end
        dy = abs(Ydat(pt2) - Ydat(pt1));
        if dy > LY/2
            dy = LY - dy;
        end
        dz= abs(Zdat(pt2) - Zdat(pt1));
        if dz > LZ/2
            dz = LZ - dz;
        end

        distMat(pt1,pt2) = sqrt(dx^2 + dy^2 + dz^2);
        distMat(pt2,pt1) = distMat(pt1,pt2);
    end 
end
end