function [rmsd, trj, vel, Ucell] = superimpose_double_index(ref, trj, index, index_ref, mass, mass_ref, vel, isdecentered)
%% superimpose
% least-squares fitting of structures
%
%% Syntax
%# rmsd = superimpose_double_index(ref, trj);
%# rmsd = superimpose_double_index(ref, trj, index_atom, index_ref);
%# rmsd = superimpose_double_index(ref, trj, index_atom, index_ref, mass, mass_ref);
%# rmsd = superimpose_double_index(ref, trj, [], mass);
%# [rmsd, trj] = superimpose_double_index(ref, trj, index_atom, index_ref, mass, mass_ref);
%# [rmsd, trj, vel] = superimpose_double_index(ref, trj, index_atom, mass, vel);
%# [rmsd, trj, vel] = superimpose_double_index(ref, trj, index_atom, [], vel);
%# [rmsd, trj, vel] = superimpose_double_index(ref, trj, [], [], vel);
%
%% Description
%
% * ref        - reference structure 
%                [double natom3]
% * trj        - trajectory fitted to the reference structure
%                [double nframe x natom3]
% * index_atom - index of atoms used in the calculation of fitting
%                [integer n]
% * index_ref  - index of atoms of the reference structure used in the calculation of fitting
%                useful when ref and trj have a different number of atoms,
%                and you want to use separate indices for each one
%                [integer n]
% * mass       - mass
%                [double natom]
% * mass_ref   - mass of the ref structure in case it has a different
%                number of atoms
%                [double natom]
% * vel        - velocity
%                [double nframe x natom3]
% * rmsd       - root mean square deviations after fitting
%                [double nframe]
% 
%% Example
%# trj = readnetcdf('ak.nc');
%# ref = trj(1,:);
%# [rmsd, trj] = superimpose(ref, trj);
%# plot(rmsd)
%
%% See also
% superimpose superimpose2d meanstructure
%
%% References
% m-file version (superimpose.m) using Kabsch's method:
%  W. Kabsch, "A solution for the best rotation to relate two sets of vectors." Acta. Cryst. A32, 922-923 (1976)
%  W. Kabsch, "A discussion of the solution for the best rotation to relate two sets of vectors." Acta. Cryst. A34, 827-828 (1978)
% mex-file version (suserimpose.c) using the QCP method
%  P. Liu, D. K. Agrafiotis, and D. L. Theobald, J. Comput. Chem. 31, 1561-1563 (2010).
% 

%% preparation
natom3 = size(trj, 2);
natom  = natom3/3;

natom3_ref = size(ref, 2);
natom_ref  = natom3_ref/3;

nframe = size(trj, 1);

if ~exist('index', 'var') || isempty(index)
  index = 1:natom;
else
  if islogical(index)
    index = find(index);
  end
  if iscolumn(index)
    index = index';
  end
end
index3 = to3(index);

if ~exist('index_ref', 'var') || isempty(index_ref)
  index_ref = natom_ref;
else
  if islogical(index_ref)
    index_ref = find(index_ref);
  end
  if iscolumn(index_ref)
    index_ref = index_ref';
  end
end
index_ref3 = to3(index_ref);

if ~exist('mass', 'var') || isempty(mass)
  mass = ones(1, natom);
else
  if iscolumn(mass)
    mass = mass';
  end
end

if ~exist('mass_ref', 'var') || isempty(mass_ref)
  mass_ref = ones(1, natom_ref);
else
  if iscolumn(mass_ref)
    mass_ref = mass_ref';
  end
end

if ~exist('vel', 'var')
  vel = [];
end

if ~exist('isdecentered', 'var')
  isdecentered = false;
end

if nargout >= 4
  Ucell = cell(nframe, 1);
end

%% remove the center of mass
if ~isdecentered
  trj = decenter(trj, index, mass);
  [ref, comy] = decenter(ref, index_ref, mass_ref);
  if numel(vel) ~= 0
    vel = decenter(vel, index, mass);
  end
else
  comy = [0 0 0];
end

mass = mass(index);
massxyz = repmat(mass, 3, 1);
y = reshape(ref(1, index_ref3), 3, numel(index_ref));
rmsd = zeros(nframe, 1);

%% superimpose
for iframe = 1:nframe
  % calculate R matrix
  x = reshape(trj(iframe, index3), 3, numel(index));
  rmsd(iframe) = 0.5 * sum(mass.*sum(x.^2 + y.^2));
  R = (massxyz.*y) * x';
  [V, D, W] = svd(R);
  D = diag(D);

  % check reflection
  is_reflection = det(V)*det(W');
  if(is_reflection < 0) 
    D(3) = -D(3);
    V(:, 3) = -V(:, 3);
  end
  rmsd(iframe) = rmsd(iframe) - sum(D);

  if nargout >= 2
    % calculate rotation matrix
    U = V*W';
    if nargout >= 4
      Ucell{iframe} = U;
    end

    % rotate molecule
    x = reshape(trj(iframe, :), 3, natom);
    x = U*x;

    % restore the original center of mass
    x(1, :) = x(1, :) + comy(1);
    x(2, :) = x(2, :) + comy(2);
    x(3, :) = x(3, :) + comy(3);
    trj(iframe, :) = reshape(x, 1, natom3);
  
    if numel(vel) ~= 0
      % rotate velocity
      v = reshape(vel(iframe, :), 3, natom);
      v = U*v;
      vel(iframe, :) = reshape(v, 1, natom3);
    end
  end
end

rmsd = sqrt(2.0*abs(rmsd)./sum(mass));

