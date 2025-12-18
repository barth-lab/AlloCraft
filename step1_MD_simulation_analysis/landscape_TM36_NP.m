%% Calculate the potential of mean force:


xi = 7.5:0.05:21; % grids in x-axis
yi = -0.3:0.02:5.5; % grids in y-axis

% Calculate the pmf
gaussian_width = [0.1 0.1]; %the bandwidth of the Gaussian kernel function
pmf_ar = calcpmf2d([vertcat(x{:}) vertcat(np{:})], xi, yi, gaussian_width);
% pmf_b4 = calcpmf2d([x_b4_mat np_b4_mat], xi, yi, gaussian_width);

s   = getconstants(); % get Boltzmann constant in kcal/mol/K
T   = 310.0;          % set temperature
pmf_ar = s.KB*T*pmf_ar;     % convert unit from KBT to kcal/mol
% pmf_b4 = s.KB*T*pmf_b4;
% Plot the landscapes

contour_levels = 0:0.25:3.5;

%figure
%subplot(1,2,1)
landscape(xi, yi, pmf_ar, contour_levels); colorbar;
h = colorbar;
% ylabel(h, 'kcal/mol')
axis([ xi(1) xi(end) yi(1) yi(end)]);
Ang = char(197);
xlabel(['TM3-6 distance (' Ang ')'], 'FontSize', 20, 'FontName', 'Helvetica');
ylabel(['NPxxY rmsd to inactive (' Ang ')'], 'FontSize', 20, 'FontName', 'Helvetica');
%title('5HT-1B')
hold on

% Calculate the average starting point
x_ar_beg = 0;
np_ar_beg = 0;

for runi=1:sim.runCount
    x_ar_beg = x_ar_beg + x{runi}(1);
    np_ar_beg = np_ar_beg + y{runi}(1);
end
x_ar_beg = x_ar_beg/sim.runCount;
np_ar_beg = np_ar_beg/sim.runCount;

% Here we plot the inactive state reference point, taken from PDB 6CM4:
scatterRed = [0.8500, 0.3250, 0.0980];
mrkredgecolor=[0.8500, 0.3250, 0.0980]+[0.15 0.3 0.3];
mrkredgecolor_green=[0.4660, 0.6740, 0.1880]-[0.4 0.4 0.188];
scatter(x_ar_beg,np_ar_beg,150,'MarkerEdgeColor',mrkredgecolor,...
              'MarkerFaceColor',scatterRed, ...
              'LineWidth',1.5)
hold on
% mrkr_end='s';
% scatter(8.74,0,150,mrkr_end,'MarkerEdgeColor',mrkredgecolor_green,...
%               'MarkerFaceColor',[0.4660, 0.6740, 0.1880]-[0.2 0.2 0.188],...
%               'LineWidth',1.5)
% hold on
 % Add the legend in a sneaky way
h = zeros(1, 1);
h(1) = scatter(NaN,NaN,150,'MarkerEdgeColor',mrkredgecolor,...
              'MarkerFaceColor',scatterRed, ...
              'LineWidth',1.5);
% h(2) = scatter(NaN,NaN,150,mrkr_end,'MarkerEdgeColor',mrkredgecolor_green,...
%               'MarkerFaceColor',[0.4660, 0.6740, 0.1880]-[0.2 0.2 0.188],...
%               'LineWidth',1.5);
% h(4) = scatter(NaN,NaN,100,'MarkerEdgeColor',mrkredgecolor_purple,...
%               'MarkerFaceColor',[0.4940    0.1840    0.5560],...
%               'LineWidth',1.5);
% legend(h,'4IAR', 'FontSize', 15, 'FontName', 'Helvetica');
% legend('boxoff')
