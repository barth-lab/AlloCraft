function [nHBonds,energies,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] = hydrogenBondManipulate(pdb, traj, forces, description, comp_mode, cut_off_max, cut_off_min, cut_off_angle, tol_res)
%hydrogenBondManipulate Performs Hydrogen bond analysis of the PDB for every
%frame in traj and visualizes the results in a dynamic map
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
%
%% Usage:
% nHBonds = hydrogenBondManipulate(pdb, traj, forces)
% nHBonds = hydrogenBondManipulate(pdb, traj, forces, description, comp_mode, cut_off_max, cut_off_min, cut_off_angle, tol_res)
% [nHBonds,energies,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] = hydrogenBondManipulate(pdb, traj, forces)
% [nHBonds,energies,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] = hydrogenBondManipulate(pdb, traj, forces, description, comp_mode, cut_off_max, cut_off_min, cut_off_angle, tol_res))
%
%% Description:
% * nHBonds are the number of hydrogen bonds detected for every frame of
% the trajectory. [nFrames x 1] array
%
% *energies are the total energies of all H-bonds at each frame. 
% [nFrames x 1] array
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
% * comp_mode determines the criteria for comparing H-bond maps between
% frames, 'occ' compares H-bond occupancy between a pair of residues, where
% an H-bond is either formed or broken in a binary fashion.  'energy' 
% compares the energy of the H-bond depending on the distance and angle.
%
% * cut_off_max is the maximum cut-off for the pair-list calculations. This
% is the distance between H and the acceptor atom. Defaults to 2.5 Ang.
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
% * pair_cell, pair_res_cell, pair_names_cell, angles_cell, and is_bb_cell 
% contain the atom pairs, residue pairs, atom names, D-H--A angles, and 1/0 
% for backbone Hbond respectively, as numbered and named in PDB. 
% nFrame x 1 cell structure, where each cell is nHBonds(frame) long.
%
% See also hydrogenBondAnalysis, hydrogenBondPeaks, hydrogenBondEnergy

% Set the default values:
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

if ~exist('description','var')
   description = 'Trajectory';
end

if ~exist('comp_mode','var')
   comp_mode = 'energy';
end

nFrames =  size(traj,1);
nHBonds = zeros(nFrames, 1); 
energies = zeros(nFrames, 1);

pair_cell = cell(nFrames, 1);
pair_res_cell = cell(nFrames, 1);
pair_names_cell = cell(nFrames, 1);
angles_cell = cell(nFrames, 1);
is_bb_cell = cell(nFrames, 1);
occ_cell = cell(nFrames, 1);
occ_bb_cell = cell(nFrames, 1);
pair_res_bb_cell = cell(nFrames, 1);
energy_cell = cell(length(nFrames),1);
energy_bb_cell = cell(length(nFrames),1);
energy_sum_res_cell = cell(length(nFrames),1);

% Variables related to energy calculation:
E_0 = -1; % For now this is a dummy placeholder value
R_0 = 2; % The position of the minimum of the energy well 

% Backbone - backbone interactions: Choose C, CA, N, O, HN, and HA only
index_ca = selectname(pdb.name, 'CA');
index_c = selectname(pdb.name, 'C');
index_n = selectname(pdb.name, 'N');
index_o = selectname(pdb.name, 'O');
index_hn = selectname(pdb.name, 'HN');
index_ha = selectname(pdb.name, 'HA');
% For general calculations, choose ALL H, N, and O's
index_Hall = selectname(pdb.name, 'H*');
index_Nall = selectname(pdb.name, 'N*');
index_Oall = selectname(pdb.name, 'O*');

index_bb = index_ca | index_c | index_n | index_o | index_hn | index_ha; % 1 for backbone, 0 for not

for frame = 1:nFrames
    % Calculate a pair list:
    [pair, dists] = calcpairlist(traj(frame,:), cut_off_max); % Calculate all atom pairs within a cutoff distance 

    % Remove any bonds that are shorter than cut_off_short (1.5 A), those are not Hbonds:
    pair1 = pair(:,1);
    pair2 = pair(:,2);
    pair1(dists<cut_off_min)=[];
    pair2(dists<cut_off_min)=[];
    pair = [pair1 pair2];

    % Remove pairs within the same or within tol_res (def = 1) residues
    same_res = zeros(size(pair,1),1);
    for row = 1:size(pair,1)
        if abs(pdb.resseq(pair(row,1)) - pdb.resseq(pair(row,2))) <= tol_res
            same_res(row) = 1;
        end
    end
    pair = pair.*not(same_res); % Zeroes out pairs on same or within "tol_res" residues
    pair1 = pair(:,1);
    pair2 = pair(:,2);
    pair1 = nonzeros(pair1);
    pair2 = nonzeros(pair2);
    pair = [pair1 pair2];

    % Set rules for accepted donor - acceptor pairs: 
    % Backbone: HN -- O, HA -- O 
    % All: H* -- O*, H* -- N*
    acc_pair = zeros(size(pair,1),1); % Will contain the pairs to accept

    for row = 1:size(pair,1)
        acc_pair(row) = index_Hall(pair(row,1))*index_Oall(pair(row,2)) + ...
            index_Hall(pair(row,2))*index_Oall(pair(row,1)) + ...
            index_Hall(pair(row,1))*index_Nall(pair(row,2)) + ...
            index_Hall(pair(row,2))*index_Nall(pair(row,1));
    end

    pair = pair.*acc_pair; % Zeroes out pairs that do not form H-bonds
    pair1 = pair(:,1);
    pair2 = pair(:,2);
    pair1 = nonzeros(pair1);
    pair2 = nonzeros(pair2);
    pair = [pair1 pair2];

    pair_res = zeros(size(pair,1),2);
    pair_names = repmat('12345678',size(pair,1),1);

    for row = 1:size(pair,1)
        pair_res(row,:) = [pdb.resseq(pair(row,1)) pdb.resseq(pair(row,2))];
        pair_names(row,:) = [pdb.name(pair(row,1),:) pdb.name(pair(row,2),:)];
    end
    
    % Last criterion: angle
    triplet = zeros(size(pair,1),1); % Triplet D-H -- A
    % Hpos is either 2 or 6
    %     ' O   HN '
    %     ' HA  O  '
    for row = 1:size(pair,1) % Rearrange "pair" variable so H is in the second column
        atom_name = pdb.name(pair(row,1),:);
        atom_name = atom_name(~isspace(pdb.name(pair(row,1),:))); % Remove white spaces
        if(atom_name(1)=='H') % H is in the 1st column 
            triplet(row,1) = pair(row,2); % Switch so H is in the 2nd column
            triplet(row,2) = pair(row,1);
        else % H is in the 2nd column
            triplet(row,1) = pair(row,1);
            triplet(row,2) = pair(row,2); 
        end
        % The atom H is attached to is always before it in
        % the PDB, either -1, -2, or -3: so
        test_temp = 1;
        donor = 1;
        atom_name_test = pdb.name(triplet(row,2) - donor,:);
        atom_name_test = atom_name_test(~isspace(atom_name_test)); % Remove white spaces
        while test_temp==1
            if atom_name_test(1) ~= 'H' % it's not a Hydrogen, we found the donor!!!
                break
            end
            donor = donor + 1;
            atom_name_test = pdb.name(triplet(row,2) - donor,:);
            atom_name_test = atom_name_test(~isspace(atom_name_test)); % Remove white spaces
        end
         triplet(row,3) = triplet(row,2) - donor; % Which will give us the donor atom!
    end

    % FINALLY CALCULATE THE FREAKING ANGLES
    angles = calcangle(traj(frame,:),triplet).*180./pi; % In degrees
    angles_acc = angles;
    % Now remove bonds beyond the angle cutoff: cut_off_angle
    accept_angle = zeros(size(pair,1),1);
    for row = 1:size(pair,1)
        if abs(angles(row)-180) <= cut_off_angle
            accept_angle(row) = 1;
        end
    end
    pair = pair.*accept_angle; % Zeroes out pairs beyond the angle cutoff
    angles_acc = angles_acc'.*accept_angle;
    angles_acc = nonzeros(angles_acc);
    pair1 = pair(:,1);
    pair2 = pair(:,2);
    pair1 = nonzeros(pair1);
    pair2 = nonzeros(pair2);
    pair = [pair1 pair2];
    dists_acc = calcbond(traj(frame,:),pair); % Calculate distances for references
    
    % Calculate the energies in this frame:
    energy_cell{frame} = hydrogenBondEnergy(dists_acc,angles_acc,E_0,R_0);
    energy_sum_res_cell{frame} = energy_cell{frame}; % Contains the sum of 
    % energies over repeated residues
    energies(frame) = sum(energy_cell{frame});
    % Calculate residue numbers again
    pair_res = zeros(size(pair,1),2);
    pair_names = repmat('12345678',size(pair,1),1);

    for row = 1:size(pair,1)
        pair_res(row,:) = [pdb.resseq(pair(row,1)) pdb.resseq(pair(row,2))];
        pair_names(row,:) = [pdb.name(pair(row,1),:) pdb.name(pair(row,2),:)];
    end
    
    nHBonds(frame) = size(pair,1); % Number of H-bonds per frame
    % Now calculate occupancies and backbone pairs:
    is_bb = zeros(size(pair,1),1); % 1 for backbone, 0 for not   
    occ =  zeros(size(pair,1),1);
    for row = 1:size(pair,1)
        A = find(pair_res(:,1)==pair_res(row,1));
        B = find(pair_res(:,2)==pair_res(row,2));
        occ(row) = length(intersect(A,B)); % how many times does this appear?
        
        if occ(row) > 1 % If it appears more than once, sum the energy over 
            %repeated residues
            energy_sum_res_cell{frame}(row) = sum(energy_cell{frame}(intersect(A,B)));
            % Note that this  variable should be used for visualization and
            % NOT for calculations
        end
        if index_bb(pair(row,1))*index_bb(pair(row,2))==1
            is_bb(row) = 1; % Backbone pair!
        end
    end
    % Make arrays with BB pairs so they can be plotted separately
    pair_res_bb = [nonzeros(pair_res(:,1).*is_bb) nonzeros(pair_res(:,2).*is_bb)];
    occ_bb = nonzeros(occ.*is_bb);
    if isrow(energy_sum_res_cell{frame}) % Make sure it's a column vector
       energy_sum_res_cell{frame} = energy_sum_res_cell{frame}'; 
    end
    energy_bb_cell{frame} = nonzeros(energy_sum_res_cell{frame}.*is_bb);
    dists_acc_bb = nonzeros(dists_acc.*is_bb);
    angles_acc_bb = nonzeros(angles_acc.*is_bb);
        
    pair_cell{frame} = pair;
    pair_res_cell{frame} = pair_res;
    pair_names_cell{frame} = pair_names;
    angles_cell{frame} = angles_acc;
    is_bb_cell{frame} = is_bb;
    occ_cell{frame} = occ;
    occ_bb_cell{frame} = occ_bb;
    pair_res_bb_cell{frame} = pair_res_bb;

end

%% Visualization section:

% Plot the manipulat-able H-bond map:

% Are we doing occupancies or energies?
if  strcmp(comp_mode,'energy') % Plot the energies!
    plot_var1 = energy_sum_res_cell;
    plot_var_bb = energy_bb_cell;
    plot_des1 = '\Delta E_{hb}';
    plot_var2 = energies;
    plot_des2 = 'Energy of H-bonds';
else % we're doing occupancies!
    plot_var1 = occ_cell;
    plot_var_bb = occ_bb_cell;
    plot_des1 = 'Occurances';
    plot_var2 = nHBonds;
    plot_des2 = 'Nb. of H-bonds';
end

f = figure;
 sz = 50; % Marker size
% Divide the space into 3 tiles, 2 for the map, and 1 for the force trace
subplot(3,1,[1 2])
% Draw a rectangle over the stachel residues
rectangle('Position', [pdb.resseq(1)-1 pdb.resseq(end)-8 pdb.resseq(end)+1 9] ...
    ,'EdgeColor',[.75 .75 .75], 'FaceColor', [.75 .75 .75]);
hold on
rectangle('Position', [pdb.resseq(end)-8 pdb.resseq(1)-1  9 pdb.resseq(end)+1] ...
    ,'EdgeColor',[.75 .75 .75], 'FaceColor', [.75 .75 .75]);

% Plot energy_cell or occ_cell
s = scatter(pair_res_cell{1}(:,1),pair_res_cell{1}(:,2),sz, plot_var1{1},'s','filled', ...
    'MarkerEdgeColor',[0 0.5 0.5], 'LineWidth',0.5);
    % Add a datatip for the energy or occupancy at a given pair
new_row = dataTipTextRow(plot_des1,plot_var1{1});
s.DataTipTemplate.DataTipRows(end+1) = new_row;
        
hold on
% Plot energy_bb_cell or occ_bb_cell
sBB = scatter(pair_res_bb_cell{1}(:,1),pair_res_bb_cell{1}(:,2),sz,'.', ...
    'MarkerEdgeColor',[0.8500, 0.3250, 0.0980]		, 'LineWidth',1.5);
    % Add a datatip for the energy or occupancy at a given pair
new_row_BB = dataTipTextRow(plot_des1,plot_var_bb{1});
sBB.DataTipTemplate.DataTipRows(end+1) = new_row_BB;
        
axis([ pdb.resseq(1)-1   pdb.resseq(end)+1   pdb.resseq(1)-1  pdb.resseq(end)+1]);
t = title([description ', Frame 1']);
h = legend([num2str(size(pair,1)) ' total H-bonds'],[num2str(size(pair_res_bb,1)) ' BB-BB H-bonds'],'location','SouthEast');
legend boxoff
xlabel('Residue number')
ylabel('Residue number')
colorbar
if strcmp(comp_mode, 'energy')
    caxis([E_0 0])
end
set(gca,'FontSize',15)

% Plot the force:
subplot(3,1,3)
smoothing = length(forces)/622; % This is a current placeholder
forces_smooth= smoothdata(forces,'movmean', smoothing); 
% Correspond 1st column of forces_smooth to the length of the trajectory
ratio_stride= ceil(forces_smooth(end,1)/size(traj,1));
yyaxis left
plot(forces_smooth(:,1)/ratio_stride,forces_smooth(:,2))
hold on
% Add the marker to follow:
% Correspond the size of forces smooth to the trajectory
ratio_size= ceil(size(forces_smooth,1)/size(traj,1));
sMarker = scatter(1,forces_smooth(ratio_size,2),sz,'MarkerEdgeColor',[0.1010    0.5450    0.7330],...
              'MarkerFaceColor',[0.3010, 0.7450, 0.9330]	,...
              'LineWidth',1.5);
xlabel('Frame Nb.')
ylabel('Pulling force [KJ/mol/nm]')

% Plot the total H-bond energy OR Nb of H-bonds per frames
yyaxis right
plot(plot_var2,'LineWidth',1)
ylabel(plot_des2)
set(gca,'FontSize',15)
sMarkerH = scatter(1,plot_var2(1),sz,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],...
              'MarkerFaceColor',[0.8500, 0.3250, 0.0980],...
              'LineWidth',1.5);

% Add the slider:
c = uicontrol(f,'Style', 'slider','units','normalized','position',[0,0,0.3,0.05],...
    'SliderStep',  [1/(nFrames-1) , 5/(nFrames-1) ], ...
    'Min', 1, ...
    'Max', nFrames, ...
    'Value',1);
c.Callback = @selection; % Callback is what you change in the plot when the slider changes

function selection(src,event) 
        val = round(c.Value); % Make sure you grab an integer
       set(src, 'Value', val);
       % Update the data with as the frame number changes
        set(s, 'XData', pair_res_cell{val}(:,1))
        set(s, 'YData', pair_res_cell{val}(:,2))
        set(s, 'CData', plot_var1{val})
        new_row = dataTipTextRow(plot_des1,plot_var1{val});
        s.DataTipTemplate.DataTipRows(end) = new_row;
        set(sBB, 'XData', pair_res_bb_cell{val}(:,1))
        set(sBB, 'YData', pair_res_bb_cell{val}(:,2))
%         set(sBB, 'CData', plot_var_bb{val})
        new_row_BB = dataTipTextRow(plot_des1,plot_var_bb{val});
        sBB.DataTipTemplate.DataTipRows(end) = new_row_BB;
        set(t, 'String', [ description ', Frame: ' num2str(val)]);
        set(h, 'String', [{[num2str(size(pair_res_cell{val},1)) ' total H-bonds']} {[num2str(size(pair_res_bb_cell{val},1)) ' BB-BB H-bonds']}]);
        set(sMarker, 'XData', val)
        set(sMarker, 'YData', forces_smooth(min(val*ratio_size, length(forces_smooth)),2))
        set(sMarkerH, 'XData', val)
        set(sMarkerH, 'YData', plot_var2(val))

end
end

