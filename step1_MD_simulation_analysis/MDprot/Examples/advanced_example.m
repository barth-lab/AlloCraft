%% Load the data:
% NOTE: before running this example, make sure to download MDtoolbox
% from https://mdtoolbox.readthedocs.io/en/latest/
% In this example we will analyze a steered molecular dynamics (SMD)
% trajectory

% Read the PDB:
[pdb, crd] = readpdb('prot.pdb');

% Read the trajectory:
traj = readdcd('traj.dcd');

% Read the forces:
forces = importdata('pullf.xvg');
%% hydrogenBondAnalysis advanced usage:

% 1- Specify the minimum and maximum distance cutoffs for the hydrogen bonds,
% as well as an angle cutoff:
cut_off_max = 2.3; % A
cut_off_min = 1.7; % A
cut_off_angle = 25; % Degrees
nHBonds = hydrogenBondAnalysis(pdb, crd, cut_off_max, cut_off_min, cut_off_angle);

%% 
% 2- If you input a trajectory rather than a coordinate file, it will
% calculate and plot the number of H-bonds per frame:

[nHBonds_traj,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] ...
    = hydrogenBondAnalysis(pdb, traj, cut_off_max, cut_off_min, cut_off_angle);


%% Visualize the shape of the energy function used to model the H-bonds:

% Let us use the first frame of traj as an example:
% Calculate the distances between atom pairs using the calcbond function:
dists = calcbond(traj(1,:),pair_cell{1});
angles = angles_cell{1}; % Just take the angles from the hydrogenBondAnalysis output

% Specify default minimum energy and equilibrium bond distance
E_0 = -1;
R_0 = 2;
[E_hb] = hydrogenBondEnergy(dists, angles, E_0 ,R_0,1);

%% Learn about the different modes of hydrogenBondPeaks
% The width factor specifies which frames the function will compare

% comp_mode specifies whether the function will compare number of H-bonds
% given a certain cut-off ('occ') or compare H-bond energies ('energy').

% Now we will test the difference ( and importance) of choosing which
% frames to compare:
% peak width mode:
width_factor = 0.5;
hydrogenBondPeaks(pdb,traj,forces,'Trajectory',width_factor,'occ')

%% if you find that the function is finding too many peaks, raise the MPP
% (minimum peak prominence) 
MPP = 200;
smoothing = 150;

hydrogenBondPeaks(pdb,traj,forces,'Trajectory',width_factor,'occ', ...
    cut_off_max,cut_off_min,cut_off_angle, 1, smoothing, MPP)

%% Peak to peak mode: width_factor = -1

width_factor = -1;
hydrogenBondPeaks(pdb,traj,forces,'Trajectory',width_factor,'occ', ...
    cut_off_max,cut_off_min,cut_off_angle, 1, smoothing, MPP)

%% Peak - peak + lag: width_factor = 0

width_factor = 0;
hydrogenBondPeaks(pdb,traj,forces,'Trajectory',width_factor,'occ', ...
    cut_off_max,cut_off_min,cut_off_angle, 1, smoothing, MPP)

%% Energy mode:
% in this mode we compare H-bond energies rather than just number of
% H-bonds, sometimes the energies tell us a different story than the
% numbers! Look at the 3rd plot for example, the N. of bonds formed is
% highe than those dissociated, but the dissociated bonds' energy is larger

width_factor = 0.5;
hydrogenBondPeaks(pdb,traj,forces,'Trajectory',width_factor,'energy', ...
    cut_off_max,cut_off_min,cut_off_angle, 1, smoothing, MPP)

%% hydrogenBondManipulate also supports 'occ' and 'energy' modes:
hydrogenBondManipulate(pdb,traj,forces,'Trajectory','occ', ...
    cut_off_max,cut_off_min)

%% In the energy mode, you can see H-bond energies fluctuating between frames
hydrogenBondManipulate(pdb,traj,forces,'Trajectory','energy', ...
    cut_off_max,cut_off_min)