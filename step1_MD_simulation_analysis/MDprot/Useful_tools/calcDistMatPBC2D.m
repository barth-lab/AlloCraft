function [distMat] = calcDistMatPBC2D(Xdat,Ydat,LX,LY)
%calcDistMatPBC Calculate distance matrix of input data taking into
% consideration periodic boundary conditions with box size [LX LY]
% Having specific case for 2D/3D makes the code ~30 times faster than
% nesting a 3D loop that accounts for variable dimensions
%
%% Usage:
%   [distMat] = calcDistMatPBC2D(Xdat,Ydat,LX,LY)
% 
%% Description
% * Xdat and Ydat are the input data
% * LX and LY are the dimensions of the periodic box
%% See also
% calcDistMatPBC

assert(length(Xdat) == length(Ydat),'Input data vectors should be the same length!');

distMat = zeros(length(Xdat));
for pt1 = 1:length(Xdat)-1

    for pt2 = pt1+1:length(Xdat)
        dx = abs(Xdat(pt2) - Xdat(pt1));
        if dx > LX/2
            dx = LX - dx;
        end
        dy = abs(Ydat(pt2) - Ydat(pt1));
        if dy > LY/2
            dy = LY - dy;
        end
        distMat(pt1,pt2) = sqrt(dx^2 + dy^2);
        distMat(pt2,pt1) = distMat(pt1,pt2) ;
    end 
end
end