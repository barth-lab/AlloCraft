function [nHBonds,peak_loc,delta_res,delta_val] = hydrogenBondPeaks(pdb, traj, forces, description, width_factor, comp_mode, cut_off_max, cut_off_min, cut_off_angle, tol_res, smoothing, MPP)
%hydrogenBondPeaks Performs Hydrogen bond analysis of the PDB for frames 
% corresponding to peaks in the force trace of a pulling simulation, such
% as steered molecular dynamics (SMD)
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
%
%% Usage:
% nHBonds = hydrogenBondPeaks(pdb, traj, forces)
% nHBonds = hydrogenBondPeaks(pdb, traj, forces, description, width_factor, comp_mode, cut_off_max, cut_off_min, cut_off_angle, tol_res, smoothing, MPP)
% [nHBonds,peak_loc,delta_res,delta_occ] = hydrogenBondPeaks(pdb, traj, forces)
% [nHBonds,peak_loc,delta_res,delta_occ] = hydrogenBondPeaks(pdb, traj, forces, description, width_factor, comp_mode, cut_off_max, cut_off_min, cut_off_angle, tol_res, smoothing, MPP)
%
%% Description:
% * nHBonds are the number of hydrogen bonds detected for every frame of
% the trajectory. [nFrames x 1] array
%
% * pdb is the pdb structure obtained by pdb = readpdb('pdb.pdb'). Note
% that other structure files (.gro for example) will not work, as the way
% mdtoolbox names the structure elements is different.
%
% * traj is the trajectory, as obtained by traj = readdcdmat(traj.dcd) for
% example. [Nframes x 3*Natoms]
%
% * forces is the force trace of the pulling (SMD) simulation corresponding
% to traj, first column contains x or time coordinate, second column
% contains the forces. [N x 2]
%
% * description is what will show as the title of the plots. [string]
%
% * width_factor specifies the method of calculating differences between
% Hbond maps at force peaks. The code supports 3 modes: 
% 1- width_factor > 0:  will calculate the Hbond maps for
%  (peak_loc +- peak_width*width_factor), and then compare the maps
% 2- width_factor = -1: will perform peak-to-peak comparison of Hbond maps  
% 3- width factor = 0 or not inserted as input: the code compares Hbond 
%  maps for peak_loc and peak_loc + 20 frames by default.
%  I recommend trying different values of width_factor (-1, 0.5, 1, and
%  1.5 for example) to find the most meaningful Hbond map.
%
% * comp_mode determines the criteria for comparing H-bond maps between
% frames, 'occ' compares H-bond occupancy between a pair of residues, where
% an H-bond is either formed or broken in a binary fashion.  'energy' 
% compares the energy of the H-bond depending on the distance and angle.
%
% * cut_off_max is the maximum cut-off for the pair-list calculations.
% Defaults to 2.5 Ang.
%
% * cut_off_min is the minimum considered distance for H-bonds, any
% distance less that this cut-off is not considered. Defaults to 1.5 Ang.
%
% * cut_off_angle is the cut-off value for the angle formed between the
% Donor-H--Acceptor, criterion for acceptance is calculated as follows: 
% abs(angle - 180) <= cut_off_angle. Defaults to 30 degrees.
% 
% * tol_res, pairs in the same residue OR within "tol_res" residues are not
% considered as H-bond forming pairs. Defaults to 1, meaning that atoms on
% the following, same, and previous residues are not considered for bonding 
% by default.
%
% * peak_loc, delta_res, and delta_occ contain the location of the force
% peaks, residue numbers (as numbered in the PDB) involved in in 
% forming or breaking Hbonds at each peak, and the difference in occupancy
% of Hbonds in these residues respectively. nPeaks x 1 cell structure, 
% where each cell is nHBonds(peak) long.
%
% * smoothing and MPP are paramters for tweaking finding the peaks in
% forces, where smoothing is performed using a moving average for the force
% trace, and MPP is the minimum peak prominence. Default to smoothing: 
% length(forces)/622 (placeholder at the moment), and MPP: 150.
%
%  See also hydrogenBondAnalysis, hydrogenBondManipulate, hydrogenBondEnergy

%% Set the default values:
if ~exist('cut_off_max','var')
    cut_off_max = 2.5; %  Angstrom
end

if ~exist('cut_off_min','var')
    cut_off_min = 1.5; % Angstrom
end

if ~exist('cut_off_angle','var')
    cut_off_angle = 30; % Degrees
end

if ~exist('tol_res','var')
   tol_res = 1;
end

if ~exist('comp_mode','var')
   comp_mode = 'energy';
end
%% This section process the force trace and identifies peaks

% Smooth the data so the peaks are clearer:
% Smoothing before finding peaks helps clean peaks that could be found in
% the noise

% Smoothing and MPP need to be balanced out manually to optimally find the peaks.
% I need to find a method that automatically finds good values for smoothing and MPP
if ~exist('smoothing','var')
    smoothing = length(forces)/622; % This is a current placeholder
end
if ~exist('MPP','var')
    MPP = 150; % Minimum peak prominence
end

forces_smooth= smoothdata(forces,'movmean', smoothing); 
ratio_stride= ceil(forces_smooth(end,1)/size(traj,1)); % factor that makes the x-axis of 
% Find the peaks
[force_peaks,peak_loc,peak_width] = findpeaks(forces_smooth(:,2),forces_smooth(:,1)/ratio_stride,'MinPeakProminence',MPP);

if ~exist('width_factor','var') ||  width_factor == 0
    width_factor = 0; % Activate peak + lag
    mode_name = 'Mode: peak + lag';
end
if width_factor > 0 % Activate peak +- width
    peak_width=peak_width*width_factor; % Scale the width using the width factor
    mode_name = ['Mode: peak +- ' num2str(width_factor) '*width'];
elseif width_factor == -1 % Activate peak-to-peak
    peak_width(1) = peak_loc(1);
    for peak=2:length(peak_loc)
       peak_width(peak) =  peak_loc(peak)-peak_loc(peak-1);
    end
    mode_name = 'Peak-to-peak';
end

figure
findpeaks(forces_smooth(:,2),forces_smooth(:,1)/ratio_stride,'MinPeakProminence',MPP, 'Annotate','extents')
xlabel('Frame Nb.', 'FontSize', 16, 'FontName', 'Helvetica')
ylabel('Pulling force [KJ/mol/nm]', 'FontSize', 16, 'FontName', 'Helvetica')
if exist('description','var') % Add title if a description was given
    title(description, 'FontSize', 20)
end
peak_loc = round(peak_loc); % round the locations to get integer indices

nSubPlots = ceil(sqrt(length(peak_loc)));
    
%% Now calculate Hbond maps for every peak and some frames after the peak,
% and then calculate the difference between them

% Initialize all the variables (lots of them I know)
res = pdb.resseq(1):pdb.resseq(end); %list of residue numbers
lag = 20; % Difference between peak and peak + lag Hbond maps
nHBonds = zeros(length(peak_loc),1);
nHBonds_lag= zeros(length(peak_loc),1);
pair = cell(length(peak_loc),1);
pair_res = cell(length(peak_loc),1);
pair_names = cell(length(peak_loc),1);
angles = cell(length(peak_loc),1);
dists = cell(length(peak_loc),1);
occ = cell(length(peak_loc),1);
occ_map = cell(length(peak_loc),1);
energy = cell(length(peak_loc),1);
energy_map = cell(length(peak_loc),1);
pair_lag = cell(length(peak_loc),1);
pair_res_lag = cell(length(peak_loc),1);
pair_names_lag = cell(length(peak_loc),1);
angles_lag = cell(length(peak_loc),1);
dists_lag = cell(length(peak_loc),1);
occ_lag = cell(length(peak_loc),1);
occ_map_lag = cell(length(peak_loc),1);
energy_lag = cell(length(peak_loc),1);
energy_map_lag = cell(length(peak_loc),1);

% Variables related to energy calculation:
E_0 = -1; % For now this is a dummy placeholder value
R_0 = 2; % The position of the minimum of the energy well 

% Plot the peak, peak - width, or peak1 (depending on the mode)
figure
lead_peak = zeros(length(peak_loc),1);
for peak=1:length(peak_loc)
    if width_factor ~= 0 % This applies in p2p and p+-w modes
        lead_peak(peak) = round(peak_loc(peak) - (max(peak_loc(peak) - peak_width(peak),1)));
    end
    subplot(nSubPlots,nSubPlots,peak)
    [nHBonds(peak),B,C,D, E] = ...
    hydrogenBondAnalysis(pdb, traj(peak_loc(peak) - lead_peak(peak),:), cut_off_max, cut_off_min, cut_off_angle, tol_res);
    title(['Frame ' num2str(peak_loc(peak)- lead_peak(peak))])
    pair{peak} = B{1};
    pair_res{peak} = C{1};
    pair_names{peak} = D{1};
    angles{peak} = E{1};
    dists{peak} = calcbond(traj(peak_loc(peak)- lead_peak(peak),:),pair{peak}); % Calculate distances for references
    
    % Calculate the energies of H-bonds
    energy{peak} = hydrogenBondEnergy(dists{peak},angles{peak},E_0,R_0);
    energy_map{peak} = zeros(length(res));
    
    % Calculate occurance of each bond between 2 residues:
    occ{peak} =  zeros(size(pair,1),1); % occurances
    occ_map{peak} = zeros(length(res));
    for row = 1:size(pair{peak},1)
        A = find(pair_res{peak}(:,1)==pair_res{peak}(row,1));
        B = find(pair_res{peak}(:,2)==pair_res{peak}(row,2));
        occ{peak}(row) = length(intersect(A,B)); % how many times does this appear?
        xRes = pair_res{peak}(row,1)-res(1) + 1;
        yRes = pair_res{peak}(row,2)-res(1) + 1;
        % Occupancy map:
        occ_map{peak}(xRes,yRes) = occ{peak}(row);
        % Energy map:
        energy_map{peak}(xRes,yRes) = energy_map{peak}(xRes,yRes) + energy{peak}(row);
    end
end

if exist('description','var') % Add title if a description was given
    sgtitle(description, 'FontSize', 20)
end
close %%%%%%%%%%%%%%%%%%%%%

% Plot the peak + lag, peaks + width, or peak2 (Depending on the mode)
figure
lag_peak = zeros(length(peak_loc),1);
for peak=1:length(peak_loc)
     if width_factor > 0 % p+-w mode
        lag_peak(peak) = round((min(peak_loc(peak) + peak_width(peak),size(traj,1))) - peak_loc(peak));
     elseif width_factor == -1 % p2p mode
        lag_peak(peak) = 0;
     elseif width_factor == 0 % p+lag mode
        if peak ~= length(peak_loc)
            if peak_loc(peak+1)-peak_loc(peak) < lag % Peaks are too close
                lag_peak(peak) = floor((peak_loc(peak+1) + peak_loc(peak))/2 - peak_loc(peak));
            else
                lag_peak(peak) = lag;
            end
        else
            if size(traj,1)-peak_loc(peak) < lag % last peak too close to end of traj
                lag_peak(peak) = floor((size(traj,1) + peak_loc(peak))/2 - peak_loc(peak));
            else
                lag_peak(peak) = lag;
            end
        end
     end
    subplot(nSubPlots,nSubPlots,peak)
    [nHBonds_lag(peak),B,C,D, E] = ...
    hydrogenBondAnalysis(pdb, traj(peak_loc(peak)+lag_peak(peak),:), cut_off_max, cut_off_min, cut_off_angle, tol_res);
    title(['Frame ' num2str(peak_loc(peak)+lag_peak(peak))])
    pair_lag{peak} = B{1};
    pair_res_lag{peak} = C{1};
    pair_names_lag{peak} = D{1};
    angles_lag{peak} = E{1};
    dists_lag{peak} = calcbond(traj(peak_loc(peak)+lag_peak(peak),:),pair_lag{peak}); % Calculate distances for references
    
    % Calculate the energies of H-bonds
    energy_lag{peak} = hydrogenBondEnergy(dists_lag{peak},angles_lag{peak},E_0,R_0);
    energy_map_lag{peak} = zeros(length(res));
    
    % Calculate occurance of each bond between 2 residues:
    occ_lag{peak} =  zeros(size(pair,1),1); % occurances
    occ_map_lag{peak} = zeros(length(res));
    for row = 1:size(pair_lag{peak},1)
        A = find(pair_res_lag{peak}(:,1)==pair_res_lag{peak}(row,1));
        B = find(pair_res_lag{peak}(:,2)==pair_res_lag{peak}(row,2));
        occ_lag{peak}(row) = length(intersect(A,B)); % how many times does this appear?
        xRes = pair_res_lag{peak}(row,1)-res(1) + 1;
        yRes = pair_res_lag{peak}(row,2)-res(1) + 1;
        occ_map_lag{peak}(xRes,yRes) = occ_lag{peak}(row);
        % Energy map:
        energy_map_lag{peak}(xRes,yRes) = energy_map_lag{peak}(xRes,yRes) + energy_lag{peak}(row);
    end
end

if exist('description','var') % Add title if a description was given
    sgtitle(description, 'FontSize', 20)
end
close %%%%%%%%%%%%%%%%%%%%%
%% Make maps of difference between Hbond map at peak + X frames and map at
% peak
figure
sz = 50; % Marker size
[X,Y] = meshgrid(pdb.resseq(1):1:pdb.resseq(end));
XA = reshape(X,[],1); % Turn it to a vector
YA = reshape(Y,[],1); % Turn it to a vector
delta_occ = cell(length(peak_loc),1);
delta_res_occ = cell(length(peak_loc),1);
delta_energy = cell(length(peak_loc),1);
delta_res_energy = cell(length(peak_loc),1);
energy_threshold = abs(E_0)/10;

% This section deals with coloring the maps:
% For energy map: make a red to blue colormap
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

% For the occupancy map
Cmat = cell(length(peak_loc),1); % Color matrix
% Reference colors to be used
Cmat_ref=[0, 0.2470, 0.5410;0, 0.4470, 0.7410;0.3010, 0.7450, 0.9330; ...
    0 0 0; 0.8500, 0.3250, 0.0980;0.6350, 0.0780, 0.1840;0.4350, 0.000, 0.0840];

for peak=1:length(peak_loc)
    % Calculate difference between occupancies of peak and peak + lag, peak
    % +- width, or peak1 and peak2 (depending on the mode)
    
    % Map differences in energy:
    map_energy_diff = energy_map_lag{peak}-energy_map{peak};
    map_energy_diff_lin = reshape(map_energy_diff,[],1); % Turn it to a vector
    % find elements above the specified energy threshold:
    ind_en_nonzero = find(abs(map_energy_diff_lin) > energy_threshold); %indices of nonzero elements
    delta_res_energy{peak} = [XA(ind_en_nonzero)  YA(ind_en_nonzero)]; 
    delta_energy{peak} = map_energy_diff_lin(ind_en_nonzero);
    
    % Map differences in occupancy:
    map_diff = occ_map{peak}-occ_map_lag{peak};
    map_diff_lin = reshape(map_diff,[],1); % Turn it to a vector
    % find nonzero elements:
    ind_nonzero = find(map_diff_lin); %indices of nonzero elements
    delta_res_occ{peak} = [XA(ind_nonzero)  YA(ind_nonzero)]; 
    delta_occ{peak} = map_diff_lin(ind_nonzero);
    
    % Define the colors:
    % If there is any delta_occ less than -3 and more than +3, put a
    % floor and ceiling  on the color range
    if  strcmp(comp_mode,'occ')
        mid = ceil(length(Cmat_ref)/2);
        Cmat_ind = max(delta_occ{peak}+mid,1);
        Cmat_ind = min(delta_occ{peak}+mid,length(Cmat_ref));
        Cmat{peak} = Cmat_ref(Cmat_ind,:);
    end
    % Plot the difference
    subplot(nSubPlots,nSubPlots,peak)
    % Draw a rectangle over the stachel residues
    rectangle('Position', [pdb.resseq(1)-1 pdb.resseq(end)-8 pdb.resseq(end)+1 9] ...
        ,'EdgeColor',[.75 .75 .75], 'FaceColor', [.75 .75 .75]);
    hold on
    rectangle('Position', [pdb.resseq(end)-8 pdb.resseq(1)-1  9 pdb.resseq(end)+1] ...
        ,'EdgeColor',[.75 .75 .75], 'FaceColor', [.75 .75 .75]);
    
    if  strcmp(comp_mode,'energy') % Plot the energies! 
        s = scatter(delta_res_energy{peak}(:,2),delta_res_energy{peak}(:,1),sz,delta_energy{peak},'s','filled', ...
                 'LineWidth',1.5);
        row = dataTipTextRow('\Delta E_{hb}',delta_energy{peak});
        s.DataTipTemplate.DataTipRows(end+1) = row;
        colormap(jet)
        hold on 
        h = zeros(2, 1); % Add a custom legend in a sneaky way
        h(1) = scatter(NaN,NaN,sz,[ 0 0 1],'s','filled');
        h(2) = scatter(NaN,NaN,sz,[1 0 0],'s','filled');
        hbonds_formed = abs(sum(delta_energy{peak}(delta_energy{peak}<0)));
        hbonds_dis = sum(delta_energy{peak}(delta_energy{peak}>0));
        legend(h, ['\Delta E_{hb}' num2str(round(hbonds_formed,1)) ' formed'],...
            ['\Delta E_{hb}' num2str(round(hbonds_dis))  ' dissociated'],'location','SouthEast')
        legend boxoff
        colormap(my_pure_colormap)
    else % Plot the occupancies!
        s = scatter(delta_res_occ{peak}(:,2),delta_res_occ{peak}(:,1),sz,Cmat{peak},'s','filled', ...
                 'LineWidth',1.5);
        row = dataTipTextRow('Occurrences',delta_occ{peak});
        s.DataTipTemplate.DataTipRows(end+1) = row;
        hold on 
        h = zeros(2, 1); % Add a custom legend in a sneaky way
        h(1) = scatter(NaN,NaN,sz,[0.3010, 0.7450, 0.9330],'s','filled');
        h(2) = scatter(NaN,NaN,sz,[0.8500, 0.3250, 0.0980],'s','filled');
        hbonds_formed = abs(sum(delta_occ{peak}(delta_occ{peak}<0)));
        hbonds_dis = sum(delta_occ{peak}(delta_occ{peak}>0));
        legend(h, [num2str(hbonds_formed) ' formed'], [num2str(hbonds_dis) ...
            ' dissociated'],'location','SouthEast')
        legend boxoff
    end
    axis([pdb.resseq(1) pdb.resseq(end) pdb.resseq(1) pdb.resseq(end)])
    xlabel('Residue number')
    ylabel('Residue number') 
    title(['Frame Diff: ' num2str(peak_loc(peak)- lead_peak(peak)) ' - ' num2str(peak_loc(peak)+lag_peak(peak))])

end

% Assign output depending on comparison mode:
if  strcmp(comp_mode,'energy') % Plot the energies! 
    delta_val = delta_energy; 
    delta_res = delta_res_energy;
else
    delta_val = delta_occ;
    delta_res = delta_res_occ;
end
if exist('description','var') % Add title if a description was given
    sgtitle([description ', ' mode_name], 'FontSize', 20)
end    
        
end

