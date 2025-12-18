function [cij] = corrMap(pdb, traj, allAtom, tol)
%corrMap Calculates Pearson correlations for a given trajectory similar 
% to what has been done in Schoeler et al. (Mapping Mechanical Force
% Propagation through Biomolecular complexes)
% It works now! I verified the correlations
% TO DO: ADD CUTOFF BASED ON DISTANCE
%
%% Usage:
% [cij] = corrMap(pdb, traj)
% [cij] = corrMap(pdb, traj, allAtom, tol)
%% Description:
% * cij: Pearson correlations of the trajectory, where:
% cij = mean(dot(dri,drj))/sqrt(mean(dot(dri,dri))*mean(dot(drj,drj)))
% N x N matrix
% 
% * pdb is the pdb structure obtained by pdb = readpdb('pdb.pdb'). Note
% that other structure files (.gro for example) will not work, as the way
% mdtoolbox names the structure elements is different.
%
% * traj is the trajectory, as obtained by traj = readdcdmat(traj.dcd) for
%  example. [Nframes x 3*Natoms]
%
% * allAtom = 1 specifies to calculate correlations for every atom of the
% system. Defaults to calulate C alphas only
%
% * tol is the tolerance for visualizing correlations. Any cij < tol will
% be considered 0. Defaults to 0

%% Preliminaries:
if ~exist('allAtom','var')
    allAtom = 0; % Defauls to using C alphas only
end

if ~exist('tol','var')
    tol = 0; % Defaults to using C alphas only
end

% Decenter the coordinates and align them:
traj_dec = decenter(traj); 

% Choose which frames to consider:
% this step is critical for alignment and for calculating mean values of dr
% which will affect the corralations.
frameEnd = size(traj,1);
traj = traj_dec(1:frameEnd,:);

index_CA = selectname(pdb.name, 'CA'); % find index of CAs
[~, traj] = superimpose(traj(1,:), traj, index_CA); % Use CAs to align

if allAtom == 0
    traj = traj(:, to3(index_CA)); % consider only CAs?
end

%% Calculate Pearson correlation coefficients:
% Each row in the traj file is a frame. We will work frame by frame
% dr = r - rmean
% Step 1: construct dr and drdrmean:
[Nframes,Natoms] = size(traj);
Natoms = Natoms/3;
rmean = mean(traj); 
dr = traj - rmean; % First piece of the puzzle
drmean2 = mean(dr.^2); 
cijt = zeros(Natoms,Natoms,Nframes); % Second piece, time-dependent Pearson correlations:
% N x N x Nframes
% Pearson correlation is: 
% cij = mean(dot(dri,drj))/sqrt(mean(dot(dri,dri))*mean(dot(drj,drj)))

% Calculate the denominator first (normalization coeff):
drdr = zeros(Nframes, Natoms); % Contains dot(dri,dri) for every atom
for frame=1:Nframes
    for atomi = 1:Natoms
        dri = dr(frame,3*atomi-2:3*atomi);
        drdr(frame,atomi) = dot(dri,dri);
    end
end
drdrmean = mean(drdr); % the normalization terms used in the denom

% Step 2: Calculate cijt and then cij:
for frame=1:Nframes
    for atomi = 1:Natoms
        dri = dr(frame,3*atomi-2:3*atomi);
        drdrmeani = drdrmean(atomi);
        for atomj = atomi:Natoms
            drj = dr(frame,3*atomj-2:3*atomj);
            drdrmeanj = drdrmean(atomj);
            cijt(atomi,atomj,frame) = dot(dri,drj)/sqrt(drdrmeani*drdrmeanj);
        end
    end
end
cij = mean(cijt,3);  % Get the mean across frames, that's your cij!

%% Visualize the correlations:
% Make my own sexy colormap:
max_o = 0.8;
Nsteps = 50;
small_step = max_o/Nsteps;
pure_blue_map = zeros(Nsteps,3);
pure_red_map = zeros(Nsteps,3);

for i=1:Nsteps
    pure_blue_map(i,:) = [0 0 1] + i.*[small_step small_step 0];
    pure_red_map(i,:) =  [1 max_o max_o] - i.*[0 small_step small_step];
end
my_pure_colormap = [pure_blue_map;1 1 1; pure_red_map];

% Draw the lovely correlations:

frame1 = 1;
frame2 = Nframes;

figure
% mat1 = cijt{frame1};
% mat1(abs(mat1)<tol) = 0;
% mat2 = cijt{frame2};
% mat2(abs(mat2)<tol) = 0;
% imagesc(mat1 + mat2')

mat1 = cij;
mat1(abs(mat1)<tol) = 0;
imagesc(mat1)
axis xy
caxis([-1 ;1])
colorbar
colormap(my_pure_colormap)
% title(['frames ' num2str(frame2) '(top) and ' num2str(frame1) '(bot)'])

end

