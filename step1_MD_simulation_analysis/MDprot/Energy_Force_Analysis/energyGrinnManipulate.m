function [energies, Xres, Yres, f] = energyGrinnManipulate(pdb, data_imp, forces, description, tol)
%energyManipulate makes an interactive visualziation of energies obtained
% from gRINN 
% This function uses the mdtoolbox package from https://mdtoolbox.readthedocs.io/en/latest/
% This function takes the output from gRINN: https://grinn.readthedocs.io/en/latest/
%% Usage:
% [energies, Xres, Yres, f] = energyGrinnManipulate(pdb, data_imp, forces, description, tol)
%% Description:
% *energies are the non-bonded CHARMM energies of residue pairs at each frame. 
%  All energies from gRINN are in kcal/mol
%  Energies as time-series data are smoothed with a window of around 1% of
%  the length of the data
%
% * Xres and Yres are the residue pair for each column of energies. For
% example: energies(:,1) will give you a time series for pairwise energy 
% between Xres(1) and Yres(1)
%
% * f is the handle for the figure
%
% * pdb is the pdb structure obtained by pdb = readpdb('pdb.pdb'). Note
% that other structure files (.gro for example) will not work, as the way
% mdtoolbox names the structure elements is different.
%
% * data_imp is the imported data from gRINN, simply run:
% data_imp = importdata('energies_intEnTotal.csv') and give the structure
% AS IS as input to the function. (The text headers are necessary!!)
%
% * forces is the force trace of the pulling (SMD) simulation corresponding
% to traj, first column contains x or time coordinate, second column
% contains the forces. [N x 2]
%
% * description is what will show as the title of the plots. [string]
%
% * tol is the energy tolerance under which energies will not be plotted
%
% See also hydrogenBondAnalysis, hydrogenBondPeaks, hydrogenBondEnergy

% Set the default values:
if ~exist('tol','var')
   tol = 5; % kcal/mol
end

if ~exist('description','var')
   description = 'Trajectory';
end

% gRINN uses kcal/mol as energy units
%Shape of data:
frames = data_imp.data(:,1); % First column contains frame numbers
nFrames = length(frames);
energies = data_imp.data(:,2:end); % 2nd to last column contains energies
res_pairs_text = data_imp.textdata(2:end); % Names of residue pairs

% Initialize energy map to zero: Nres x Nres
Nres = length(unique(pdb.resseq));
energy_map = cell(length(frames),1);
Xres = zeros(size(energies,2),1);
Yres = zeros(size(energies,2),1);

% Extract the residue pair numbers from the column header:
for pair =1:size(energies,2)
    str = res_pairs_text(pair); % String containing residues
    pair_double = str2double(regexp(str{1},'\d+','match')); % [Xres Yres]
    Xres(pair) = pair_double(1);
    Yres(pair) = pair_double(2);
end

% Smooth the energy before filling the energy map:
% Smoothing window is 1% of the data length
smoothing_window = ceil(nFrames/100); 
energies_smooth = movmean(energies,smoothing_window);
energies = energies_smooth;

% Fill the energy map for each frame:
for frame = 1:length(frames)
    energy_map{frame} = zeros(Nres); % Initialize this cell to a map of zeros
    for pair =1:size(energies,2) % N of interacting pairs
        % Fill the energy map
        energy_map{frame}(Xres(pair),Yres(pair)) = energies(frame,pair); 
    end
end

%% Visualization section:

% Plot the manipulat-able energy map:
% Turn energies below tol to NaN so they don't show in the plot
energies_plot = energies;
energies_plot(abs(energies_plot)<=tol) = NaN;
%temp(abs(temp(:,3))<=tol,3) = NaN;

% Initialize plotting variables
plot_var1 = energies_plot;
plot_des1 = 'IE:';
% plot_var2 = energies;
plot_des2 = 'Non-bonded Energies';

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


% Plot energy
s = scatter(Xres,Yres,sz,plot_var1(1,:),'s','filled', ...
                  'LineWidth',1.5);
    % Add a datatip for the energy or occupancy at a given pair
new_row = dataTipTextRow('IE',energies(1,:));
s.DataTipTemplate.DataTipRows(end+1) = new_row;
   
hold on
        
axis([ pdb.resseq(1)-1   pdb.resseq(end)+1   pdb.resseq(1)-1  pdb.resseq(end)+1]);
t = title([description ', Frame 1, E_{cutoff} = ' num2str(tol)]);
% h = legend([num2str(size(pair,1)) ' total H-bonds'],[num2str(size(pair_res_bb,1)) ' BB-BB H-bonds'],'location','SouthEast');
% legend boxoff
xlabel('Residue number')
ylabel('Residue number')
hcb = colorbar;
hcb.Label.String = 'Interaction Energy [kcal/mol]';
% if strcmp(comp_mode, 'energy')
%     caxis([E_0 0])
% end
set(gca,'FontSize',15)

if exist('forces','var')
    % Plot the force:
    subplot(3,1,3)
    smoothing = length(forces)/622; % This is a current placeholder
    forces_smooth= smoothdata(forces,'movmean', smoothing); 
    % Correspond 1st column of forces_smooth to the length of the trajectory
    ratio_stride= ceil(forces_smooth(end,1)/nFrames);
    % yyaxis left
    plot(forces_smooth(:,1)/ratio_stride,forces_smooth(:,2))
    hold on
    % Add the marker to follow:
    % Correspond the size of forces smooth to the trajectory
    ratio_size= ceil(size(forces_smooth,1)/nFrames);
    sMarker = scatter(1,forces_smooth(ratio_size,2),sz,'MarkerEdgeColor',[0.1010    0.5450    0.7330],...
                  'MarkerFaceColor',[0.3010, 0.7450, 0.9330]	,...
                  'LineWidth',1.5);
    xlabel('Frame Nb.')
    ylabel('Pulling force [KJ/mol/nm]')
end
% Plot the total H-bond energy OR Nb of H-bonds per frames
% yyaxis right
% plot(plot_var2,'LineWidth',1)
% ylabel(plot_des2)
% set(gca,'FontSize',15)
% sMarkerH = scatter(1,plot_var2(1),sz,'MarkerEdgeColor',[0.6350, 0.0780, 0.1840],...
%               'MarkerFaceColor',[0.8500, 0.3250, 0.0980],...
%               'LineWidth',1.5);

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
%         set(s, 'XData', pair_res_cell{val}(:,1))
%         set(s, 'YData', pair_res_cell{val}(:,2))
        set(s, 'CData', plot_var1(val,:))
        new_row = dataTipTextRow(plot_des1,energies(val,:));
        s.DataTipTemplate.DataTipRows(end) = new_row;
        set(t, 'String', [ description ', Frame: ' num2str(val) ', E_{cutoff} = ' num2str(tol)]);
%         set(h, 'String', [{[num2str(size(pair_res_cell{val},1)) ' total H-bonds']} {[num2str(size(pair_res_bb_cell{val},1)) ' BB-BB H-bonds']}]);
        if exist('forces','var')
        set(sMarker, 'XData', val)
        set(sMarker, 'YData', forces_smooth(min(val*ratio_size, length(forces_smooth)),2))
        end
%         set(sMarkerH, 'XData', val)
%         set(sMarkerH, 'YData', plot_var2(val))

end
end

