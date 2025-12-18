% This module assesses the convergence of the system's configuational
% entropies by calculating the first and second order entropies at
% different points of the simulation.
function assessEntropyConvergence(simulation, options)
  arguments
    simulation
    options.Name = simulation.entry.name
    options.PointCount = 50
    options.SavePath
    options.T = 310; % Kelvin
    options.dihCountS2 = 100; % Number of dihedrals (sorted by S2) to consider for S2
  % convergence test
    options.convCriterion = 0.05 % Fraction of acceptable S2 change in the last options.framesConv frames
    options.framesConv = 500
  end

  r =  1.987e-3; % kcal⋅K−1⋅mol−1 -> Gas constant

  [~, dihOrder] = sort(sum(simulation.mi), 'descend');

  % Randomly shuffle rows
  P = randperm(size(simulation.dihedralsMat,1));
  data = zeroStretchtotwopi(simulation.dihedralsMat(P, :));

  frameCount = size(data, 1);

  endFrameIndices = round(linspace(1, frameCount, options.PointCount + 1));
  endFrameIndices = endFrameIndices(2:end);

  plotData = zeros(length(endFrameIndices), 2);

  for endFrameIndicesIndex = 1:length(endFrameIndices)
    endFrameIndex = endFrameIndices(endFrameIndicesIndex);
    s1 = computeEntropies(data(1:endFrameIndex, :));
    [~, s2, ~] = computeEntropies(data(1:endFrameIndex, dihOrder(1:options.dihCountS2)));
    plotData(endFrameIndicesIndex, 1) = sum(s1);
    plotData(endFrameIndicesIndex, 2) = sum(s2, 'all');
  end

  figure;
  hold on;
  plot(endFrameIndices, -r * options.T * plotData(:, 1), '-o', 'LineWidth', 1.5,'MarkerSize',3);
  ylabel("-TS^{(1)} [kcal mol^{-1}]");

  yyaxis right;
  plot(endFrameIndices, -r * options.T * plotData(:, 2), '-o', 'LineWidth', 1.5,'MarkerSize',3);

  xlabel("Number of concatenated frames");
  ylabel("-TS^{(2)} [kcal mol^{-1}]");

  title(sprintf("Convergence of %s", options.Name));
    
  % plot two horizantal lines
  s2End = -r * options.T *plotData(end, 2);
  axisLimits = axis;
  convLimit =  s2End - options.convCriterion*( s2End);
  line(axisLimits(1:2),[convLimit convLimit],'Color','red','LineStyle','--')
  line(axisLimits(1:2),[s2End s2End],'Color','red','LineStyle','-')


    % Plot vertical line 
  line([frameCount-options.framesConv frameCount-options.framesConv],axisLimits(3:4),'Color','magenta','LineStyle',':')

  legend("S_1", ['S_2: top ' num2str(options.dihCountS2) ' dihedrals'], ['\DeltaS2 = ' num2str(options.convCriterion*100) '%' ], ...
      '',['\DeltaFrames = ' num2str(options.framesConv)], 'Location', 'best');
  legend boxoff;

  A = axis;
  axis([A(1) frameCount A(3) A(4)])
  set(gca, 'FontSize', 16);

  if isfield(options, 'SavePath')
    savefig(options.SavePath);
    print2pdf(options.SavePath);
  end
end
