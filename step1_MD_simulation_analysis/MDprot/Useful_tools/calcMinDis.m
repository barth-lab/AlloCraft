function [minDis, minDisVec] = calcMinDis(trj, set1, set2)
%% calcMinDis
% calculate minimum distance between set1 and set 2 from their Cartesian 
% coordinates
%
%% Syntax
%# minDis = calcMinDis(trj, set1, set2)
%# [minDis, minDisVec] = calcMinDis(trj, set1, set2)
%
%% Description
% Calculate distances from the input trajectory of Cartesian coordinates.
% Pairs, whose distances are calculated, can be specified via the
% variable (pair).
%
% * trj    - coordinates of atoms [nframe x natom3]
% * set1, set2 - indices of groups whose distances are calculated [nset]
% * minDis   - minimum distances between the pairs [nframe x 1]
% * minDisVec -  minimum distance vector from set1 to set 2 [nframe x 3]
%% Example
%# trj = readnetcdf('ala.nc');
%# pair = [5 6; 9 10];
%# bond = calcbond(trj, pair);
%
%% See alo
% calcbond, calcangle, calcdihedral
% 
% Mahdi Hijazi, 2022

%% initialization


nframe = size(trj, 1);

%% calculation

temptrj = trj(:, to3(set1));
tempMat  = reshape(trj(:, to3(set2)),nframe,3,[]);

tempDis = zeros(nframe,sum(set1),sum(set2));
for natom = 1:(size(temptrj,2)/3)
   tempDis(:,natom,:) = sqrt(sum(( temptrj(:,(3*(natom-1)+1):3*natom) -tempMat).^2, 2));
end
% Take the minimum along dimensions 2 (set1) and 3 (set2)
[temp1, minNdx1] = min(tempDis,[],2);
[minDis, minNdx2]  = min(temp1,[],3);

if nargout > 1 % User asked for vectors
    minNdx3 = zeros(nframe,1);
    for i=1:nframe
       minNdx3(i) =  minNdx1(i,:,minNdx2(i));
    end
    
    % XYZ for set2
    xyz2 = zeros(nframe,3);
    for i=1:nframe
       xyz2(i,:) =  tempMat(i,:,minNdx2(i));
    end
    
    % XYZ for set1
    xyz1 = zeros(nframe,3);
    for i=1:nframe
       xyz1(i,:) =  temptrj(i,to3(minNdx3(i)));
    end
    minDisVec = xyz2 - xyz1;
end

end


