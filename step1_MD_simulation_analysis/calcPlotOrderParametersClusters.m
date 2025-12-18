function [params, ndxXY, clusterNdx_filtered] = calcPlotOrderParametersClusters(database, refEntry, clusterNdx, options)
  arguments
    database
    refEntry
    clusterNdx % For example, PCA cluster index: indexOfCluster_pca
    options.SaveName = "%s"
    options.SavePath
    options.frames2skip = 1 % 1 does not skip any frames, ironically
    options.nClusters =  min(max(clusterNdx),3) % Defaults to plotting 3 clusters, or less if there are less than 3
    options.cullFrames = false  % Remove frames from simulation according to user input 
  end

  res3 = database.findResidue("3.50");
  res6 = database.findResidue("6.30");
  res7 = database.findResidue("7.53");


  % Parameters

  params.npxxy.data = database.calcRmsd(database.findFeature("7.53", 4, 0), refEntry);
  params.npxxy.label = ["RMSD of NPxxY with respect", "to inactive (" + refEntry.name + ") (Å)"];

  params.tm36.data = database.calcDistance(res3, res6);
  params.tm36.label = "TM3-6 distance (Å)";

  params.tm37.data = database.calcDistance(res3, res7);
  params.tm37.label = "TM3-7 distance (Å)";


  % Plots
  plots = cell(2, 1);
  plots{1}.x = params.tm36;
  plots{1}.y = params.npxxy;
  plots{1}.title = "TM3-6 vs NPxxY";
  plots{1}.name = "param_tm36_npxxy_clusters";

  plots{2}.x = params.tm36;
  plots{2}.y = params.tm37;
  plots{2}.title = "TM3-6 vs TM3-7";
  plots{2}.name = "param_tm36_tm37_clusters";

  ndxXY = cell(2,1); % Index for state culling
%   clusterNdx_filtered= cell(2,1); % Cluster index for state culling
  for currentPlotIndex = 1:size(plots)
    currentPlot = plots{currentPlotIndex};

    tempX = cell(size(currentPlot.x.data{1},1),1);
    tempY = cell(size(currentPlot.x.data{1},1),1);
    for thisRun=1:size(currentPlot.x.data{1},1)
        tempX{thisRun} = currentPlot.x.data{1}{thisRun}((options.frames2skip):end);
        tempY{thisRun} = currentPlot.y.data{1}{thisRun}((options.frames2skip):end);
    end
    mainX = cell2mat(tempX(:));
    mainY = cell2mat(tempY(:));
%     mainX = [];
%     mainY = [];
%     for thisRun = 1:size(tempX,1) % Make mainX into a linear vector
%         mainX = [ mainX tempX{thisRun}];
%         mainY = [ mainY tempY{thisRun}];
%     end

    xl = minmaxAll([mainX; cell2mat(currentPlot.x.data(2:end))']);
    yl = minmaxAll([mainY; cell2mat(currentPlot.y.data(2:end))']);

    figure('Name', currentPlot.title, 'Position', [10 10 1500 600]);
    tiledlayout(1, 3); % Plot the top 3 clusters

    xi = xl(1):0.05:xl(2);
    yi = yl(1):0.02:yl(2);

    gaussianWidth = [0.1 0.1];
    
    for thisCluster = 1:options.nClusters
    nexttile;
%      scatter(mainX, mainY, 10, clusterNdx, 'o','filled');
    pmf = calcpmf2d([mainX(clusterNdx==thisCluster), mainY(clusterNdx==thisCluster)], xi, yi, gaussianWidth);
    % contourLevels = linspace(min(pmf, [], 'all'), max(pmf, [], 'all'), 10);
    contourLevels = 0:0.25:3.5;
    
  
    

    landscape(xi, yi, pmf, contourLevels);
    hold on;
    
    if thisCluster ==3  
        colorbar;    
    end
    
    handles = zeros(1, length(database.entries) - 1);

    for entryIndex = 2:length(database.entries)
      h = scatter(currentPlot.x.data{entryIndex}, currentPlot.y.data{entryIndex}, 100, 'o', 'DisplayName', database.entries{entryIndex}.name);

      set(h, 'MarkerFaceColor', get(h, 'MarkerEdgeColor'));
      set(h, 'MarkerEdgeColor', 'w');

      handles(entryIndex - 1) = h;
    end

    xlabel(currentPlot.x.label);
    if thisCluster==1
        ylabel(currentPlot.y.label);
    
        legend(handles, 'location','best');
        legend boxoff;
    end
    title(['Cluster C' num2str(thisCluster)])
    xlim(xl);
    ylim(yl);
    formatplot2;

    end


    sgtitle(currentPlot.title, 'FontName','Helvetica', 'FontSize', 25);

    if isfield(options, 'SavePath')
      figPath = fullfile(options.SavePath, sprintf(options.SaveName, currentPlot.name));
      savefig(figPath);
      print2pdf(figPath);
    end

     if options.cullFrames % Cull frames according the user input:
        activeLimits = repmat([0 20 0 20], options.nClusters, 1); % Hardcoded limits as requested by user
        %activeLimits = input("Input X and Y limits for " + plots{1}.title +" for " + num2str(options.nClusters) ... 
        %     +" clusters; Input format: [XminC1 XmaxC1 YminC1 YmaxC1; XminC2 XmaxC2 YminC2 YmaxC2] ...");

        figure( 'Position', [10 10 1500 600]); % May need some tinkering
        tiledlayout(1,options.nClusters);
        clusterNdx_filtered(:,currentPlotIndex) = zeros(size(clusterNdx));

        ndxXY{currentPlotIndex} = zeros(length(mainX), options.nClusters); % Index for frames within limits

        for thisCluster = 1:options.nClusters
        
            ndx1 = mainX < activeLimits(thisCluster,2) & mainX > activeLimits(thisCluster,1);
            ndx2 = mainY < activeLimits(thisCluster,4) & mainY > activeLimits(thisCluster,3);
            
            ndxXY{currentPlotIndex}(:,thisCluster) = ndx2 & ndx1; 

            % Plot excluded frames
            nexttile
            temp = ndxXY{currentPlotIndex}(:,thisCluster) & clusterNdx==thisCluster;
            clusterNdx_filtered(temp,currentPlotIndex) = thisCluster;

            scatter(mainX(clusterNdx==thisCluster),mainY(clusterNdx==thisCluster),5,'filled')
            
            hold on
            scatter(mainX(temp),mainY(temp),5,'filled')

            xlim(xl);
            ylim(yl);
            legend(['All frames: ' num2str(sum(clusterNdx==thisCluster))], ['Included frames: ' num2str(sum(temp))])
            legend boxoff
            title(['Cluster C' num2str(thisCluster)])
            formatplot2;
            if isfield(options, 'SavePath')
              figPath = fullfile(options.SavePath, sprintf(options.SaveName, currentPlot.name + "_excluded_frames"));
              savefig(figPath);
              print2pdf(figPath);
            end
        end
     else % All frames are taken

        clusterNdx_filtered(:,currentPlotIndex) = clusterNdx;
     end
  end
end


%% Functions

function output = minmaxAll(arr)
  output = [min(arr, [], 'all'), max(arr, [], 'all')];
end
