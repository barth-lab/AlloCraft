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
%% Try the functions!
% Analyze the static structure of the protein
hydrogenBondAnalysis(pdb,crd);

%% Analyze the trajectory, find force peaks, and compare bond gain or loss
% of energy around force peaks
hydrogenBondPeaks(pdb,traj,forces);

%% Make an interactive H-bond map where you can visualize every frame of the trajectory:
cut_off_max = 2.3;
cut_off_min = 1.7;
hydrogenBondManipulate(pdb,traj,forces,'Example', 'energy' , cut_off_max, cut_off_min) 

% Move the slider in the bottom left of the figure to change the frame!
