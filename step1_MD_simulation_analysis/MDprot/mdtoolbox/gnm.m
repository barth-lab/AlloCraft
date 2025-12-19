function [emode, frequency, covar, covar_norm] = gnm(crd, cutoff)
%% gnm
% calculate normal modes and fluctuations by using Gaussian (Elastic) Network Model.
%
%% Syntax
%# emode = gnm(crd);
%# emode = gnm(crd, cutoff);
%# [emode, frequency] = gnm(crd, cutoff);
%# [emode, frequency, covar] = gnm(crd, cutoff);
%# [emode, frequency, covar, covar_norm] = gnm(crd, cutoff);
%% Description
% This routine performs a normal mode analysis of the Gaussian 
% (Elastic) Network Model from the input coordinates, and returns
% the normal modes (eigenvectors) and the corresponding frequencies
% (sqrt(eigenvalues)).
%
% Also, the covariance matrix of the Cartesian xyz coordinates is
% calculated from the pseudo inverse of Kirchoff matrix.
% In the calculation, the spring constant of 1.0 kcal/mol/A^2, 
% and kBT = 1.0 kcal/mol are assumed.
%
% * crd         - coordinates of atoms. 
%                 [double natom3]
% * cutoff      - cutoff distance of the model. default is 7.0 angstrom. 
%                 [double scalar]
% * emode       - normal modes, each column vector corresponds to a mode.
%                 1st column vector emode(:,1) is the lowest frequency mode.
%                 GNM produces natom - 1 modes.
%                 [double natom x natom]
% * frequency   - frequencies of the normal modes, given in ascending order.
%                 frequency(end) is zero, since GNM produces N-1 modes.
%                 [double natom]
% * covar       - covariance matrix of gaussian fluctuations. 
%                 covar(i,j) = <q_i q_j>
%                 [double natom x natom]
% * covar_norm  - normalized covariance matrix.
%                 covar_norm(i,j) = <q_i q_j>/(<q_i q_i><q_j q_j>)^0.5
%% Example
%# % normal mode analysis of Ca-based model of Lysozyme
%# [pdb, crd] = readpdb('lys.pdb');
%# index_ca = selectname(pdb.name, 'CA');
%# pdb = substruct(pdb, index_ca);
%# crd = crd(to3(index_ca));
%# crd = decenter(crd);
%# [emode, frequency, covar, covar_atom] = gnm(crd, 8.0);
%#
%# % rmsf and covariances
%# plot(diag(covar_atom)); xlabel('residue','fontsize',40); ylabel('variances [a.u.]','fontsize',40); formatplot
%# imagesc(covar_atom); axis xy; xlabel('residue','fontsize',40); ylabel('residue','fontsize',40); colorbar; formatplot2
%#
%# % mode structure
%# pdb.xyz = reshape(crd, 3, [])';
%# writepdb('lys_ca.pdb', pdb);
%# crd1 = crd + emode(:, 1)'*50;
%# pdb.xyz = reshape(crd1, 3, [])';
%# writepdb('lys_ca1.pdb', pdb);
%# % visualize the mode structure with PyMOL and modevectors.py
% 
%% See also
% anm, anm_sym, transformframe, writevmdvector
%
%% References
% Flory,P. (1976) Statistical thermodynamics of random networks. Proc. 
% R. Soc. Lond. A, 351, 351-380 
% 
% Bahar,I., Atilgan,A.R. and Erman,B. (1997) Direct evaluation of thermal 
% fluctuations in proteins using a single-parameter harmonic potential. 
% Fold. Des., 2, 173-181.

%% initial setup
KBT = 1.0;
gam = 1.0; % spring constant
x = crd(1,1:3:end);
y = crd(1,2:3:end);
z = crd(1,3:3:end);
natom = length(x);

if ~exist('cutoff', 'var')
  cutoff = 7.0;
end

%% make distance matrix S
S = zeros(natom,natom);
for i = 1:natom
  x_diff = x - x(i);
  y_diff = y - y(i);
  z_diff = z - z(i);
  S(i,:) = sqrt(x_diff.^2 + y_diff.^2 + z_diff.^2);
  S(i,i) = realmax; % So it does not count as a hit in the Kirchoff matrix
end

%% make kirchhoff matrix G
G = zeros(natom,natom);
for i = 1:natom
  G(i,:) = - (S(i,:) < cutoff);
  G(i,i) = - sum(G(i,:));
end

%% diagonalize kirchhoff matrix G by singular value decomposition
[~, s, v] = svd(G);
emode = v(:,(end-1):-1:1);
s = diag(s);
s = sqrt(gam*s);
% frequency = (spring cst * eigenvalue)^0.5
frequency = s((end-1):-1:1);

%% diagonalize kirchhoff matrix G by  eigenvalue decomposition (and check 
% for any difference perhaps

% [Evector,DEvalue] = eig(G);
% emode = Evector(:,1:(end-1));
% Evalue = diag(DEvalue);
% Evalue = sqrt(Evalue);
% frequency = Evalue(1:(end-1));
%% calculate covariance (= pseudo inverse of kirchhoff matrix)
if nargout >= 3 % Output covariance
  covar = 3 * (KBT/gam) * emode * diag(1./(frequency.^2)) * (emode');
end

if nargout >= 4 % Output normalized covariance
  covar_norm = covar;
  for i=1: size(covar,1)
      for j=1:size(covar,2)
          covar_norm(i,j) = covar(i,j)/sqrt(covar(i,i)*covar(j,j));
      end
  end
end
%% add the 1 external degree of freeedon in the end of output
emode = [emode v(:,end)];
frequency = [frequency; s(end)];

