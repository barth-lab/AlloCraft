function [E_hb] = hydrogenBondEnergy(dists,angles,E_0,R_0, to_plot, pairs)
%hydrogenBondEnergy Calculates hydrogen bond energies based on a simple
%function depending on bond distances and angles
%% Usage:
% [E_hb] = hydrogenBondEnergy(dists, angles)
% [E_hb] = hydrogenBondEnergy(dists, angles, E_0, R_0, to_plot, pairs)
%% Description:
% * E_hb is a vector containing the energies of the bonds 
%
% * dists are the input distances between the hydrogen and the acceptor.
%
% * angles are the angles formed between D-H--A
%
% * E_0 is the minimum energy value of the bond, defaults to -1.
%
% * R_0 is the equilibrium distance between the hydrogen and the acceptor.
% Defaults to 2 A.
%
% * to_plot,when present, plots the shape of the well as a 2D contour plot,
% and the shows the positions of the H-bonds on the energy well.
%
% * pair are the residue pairs between which the H-bonds are formed, this
% variable is used to add a tooltip to the plot.
%
%  See also hydrogenBondAnalysis, hydrogenBondManipulate, hydrogenBondPeaks


% Set default values
if ~exist('E_0','var')
    E_0 = -1; % For now this is a dummy placeholder value
end
if ~exist('R_0','var')
    R_0 = 2; % The position of the minimum of the energy well
end
if ~exist('to_plot','var')
    to_plot = 0;
end
if ~exist('pairs','var')
    pairs = 0;
end
% Make sure both dists and angles are row vectors:
if isrow(dists) == 0
   dists = dists'; % Transpose!
end
if isrow(angles) == 0
   angles = angles'; % Transpose! 
end
% Model hydrogen bond energy from PM7: 
% Stewart, J. J. (2013). Optimization of parameters for semiempirical 
% methods VI: more modifications to the NDDO approximations and 
% re-optimization of parameters. Journal of molecular modeling, 19(1), 1-32.
E_hb =E_0*((cosd(angles)).^4).*exp(-80.*(dists - R_0).^2);

if to_plot == 1
    figure
    % Grid for plotting
    r = linspace(R_0-0.4,R_0+0.4); 
    th = linspace(90,270); 
    [R,TH] = meshgrid(r,th);
    E_hb_grid = E_0*((cosd(TH)).^4).*exp(-80.*(R - R_0).^2);
    % Plot the grid:
    contourf(R,TH,E_hb_grid,50,'LineColor','none')
     colormap(flipud(jet))
     colorbar
     hold on
     s = scatter(dists,angles,[],E_hb,'filled','MarkerEdgeColor',[1 1 1]);
     xlabel('d_{H-A}')
     ylabel('\theta_{DH-A}')
     if pairs ~= 0
         pair_res_cellpairs = num2cell(pairs,2);
         row = dataTipTextRow('Residue pairs',pair_res_cellpairs);
         s.DataTipTemplate.DataTipRows(end+1) = row;
     end
     s.DataTipTemplate.DataTipRows(1).Label= 'd_{H-A}';
     s.DataTipTemplate.DataTipRows(2).Label= '\theta_{DH-A}';
end
end

