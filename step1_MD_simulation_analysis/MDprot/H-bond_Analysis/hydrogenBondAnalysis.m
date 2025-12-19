function [nHBonds,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] = hydrogenBondAnalysis(pdb, traj, cut_off_max, cut_off_min, cut_off_angle, tol_res)
%hydrogenBondAnalysis Performs Hydrogen bond analysis of the PDB for every
%frame in traj.
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
%
%% Usage:
% nHBonds = hydrogenBondAnalysis(pdb, traj)
% nHBonds = hydrogenBondAnalysis(pdb, traj, cut_off_max, cut_off_min, cut_off_angle, tol_res)
% [nHBonds,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] = hydrogenBondAnalysis(pdb, traj)
% [nHBonds,pair_cell,pair_res_cell,pair_names_cell, angles_cell, is_bb_cell] = hydrogenBondAnalysis(pdb, traj, cut_off_max, cut_off_min, cut_off_angle, tol_res))
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
% See also hydrogenBondPeaks, hydrogenBondManipulate

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

nFrames =  size(traj,1);
nHBonds = zeros(nFrames, 1);

if nargout>1 % Did the user ask for the detailed pair info?
    pair_cell = cell(nFrames, 1);
    pair_res_cell = cell(nFrames, 1);
    pair_names_cell = cell(nFrames, 1);
    angles_cell = cell(nFrames, 1);
    is_bb_cell = cell(nFrames, 1);
end

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

    % Calculate residue numbers again
    pair_res = zeros(size(pair,1),2);
    pair_names = repmat('12345678',size(pair,1),1);

    for row = 1:size(pair,1)
        pair_res(row,:) = [pdb.resseq(pair(row,1)) pdb.resseq(pair(row,2))];
        pair_names(row,:) = [pdb.name(pair(row,1),:) pdb.name(pair(row,2),:)];
    end
    
    nHBonds(frame) = size(pair,1); % Number of H-bonds per frame
    
    if frame == 1 % Draw an H bond map  
        is_bb = zeros(size(pair,1),1); % 1 for backbone, 0 for not
        % Frame 1 H bond map:
        sz = 50; % Marker size
        occ =  zeros(size(pair,1),1);
        for row = 1:size(pair,1)
            A = find(pair_res(:,1)==pair_res(row,1));
            B = find(pair_res(:,2)==pair_res(row,2));
            occ(row) = length(intersect(A,B)); % how many times does this appear?
            if index_bb(pair(row,1))*index_bb(pair(row,2))==1
                is_bb(row) = 1; % Backbone pair!
            end
        end
        % Make arrays with BB pairs so they can be plotted separately
        pair_res_bb = [nonzeros(pair_res(:,1).*is_bb) nonzeros(pair_res(:,2).*is_bb)];
        occ_bb = nonzeros(occ.*is_bb);
        dists_acc_bb = nonzeros(dists_acc.*is_bb);
        angles_acc_bb = nonzeros(angles_acc.*is_bb);
   
        % Plot ALL the Hbonds
%         figure
        % Draw a rectangle over the stachel residues
        rectangle('Position', [pdb.resseq(1)-1 pdb.resseq(end)-8 pdb.resseq(end)+1 9] ...
            ,'EdgeColor',[.75 .75 .75], 'FaceColor', [.75 .75 .75]);
        hold on
        rectangle('Position', [pdb.resseq(end)-8 pdb.resseq(1)-1  9 pdb.resseq(end)+1] ...
            ,'EdgeColor',[.75 .75 .75], 'FaceColor', [.75 .75 .75]);


        s = scatter(pair_res(:,1),pair_res(:,2),sz,occ,'s','filled', ...
            'MarkerEdgeColor',[0 0.5 0.5], 'LineWidth',0.5);
        hold on
%         set(gca,'FontSize',20)
        % Add tooltips for distance, angle, and occurance
        row = dataTipTextRow('Distance',dists_acc);
        s.DataTipTemplate.DataTipRows(end+1) = row;

        row = dataTipTextRow('Angle',angles_acc);
        s.DataTipTemplate.DataTipRows(end+1) = row;

        row = dataTipTextRow('Occurrences',occ);
        s.DataTipTemplate.DataTipRows(end+1) = row;
        % Plot the BB-BB Hbond pairs
        s = scatter(pair_res_bb(:,1),pair_res_bb(:,2),sz-10,'x','filled', ...
            'MarkerEdgeColor',[0.9290, 0.6940, 0.1250]		, 'LineWidth',1.5);
        % Add tooltips for distance, angle, and occurance
        row = dataTipTextRow('Distance',dists_acc_bb);
        s.DataTipTemplate.DataTipRows(end+1) = row;

        row = dataTipTextRow('Angle',angles_acc_bb);
        s.DataTipTemplate.DataTipRows(end+1) = row;

        row = dataTipTextRow('Occurrences',occ_bb);
        s.DataTipTemplate.DataTipRows(end+1) = row;
        xlabel('Residue number')
        ylabel('Residue number')
        axis([pdb.resseq(1)-1 pdb.resseq(end)+20 pdb.resseq(1)-1 pdb.resseq(end)+20])
        
        legend([num2str(size(pair,1)) ' total H-bonds'],[num2str(size(pair_res_bb,1)) ' BB-BB H-bonds'],'location','SouthEast')
        legend boxoff
    end
    
    if nargout>1 % Did the user ask for the detailed pair info?
    pair_cell{frame} = pair;
    pair_res_cell{frame} = pair_res;
    pair_names_cell{frame} = pair_names;
    angles_cell{frame} = angles_acc;
    is_bb_cell{frame} = is_bb;
    end
end
% Visualization section:
% Plot the Nb of H-bonds vs frames
if nFrames > 1
    figure
    plot(nHBonds,'LineWidth',1)
    xlabel('Frame Nb.')
    ylabel('Nb. of H-bonds')
    set(gca,'FontSize',20)     
end
end

