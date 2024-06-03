function [excessMi, binMiArray, iiArray, significanceArray] = filterCorrectMI(dihedralsMat, mi, options)
    arguments
        dihedralsMat
        mi
        options.BinCount = 50
        options.DihTypes = []
        options.InterpBinCount = 100
        options.InterpPointCountPerBin = 10
        options.PermutationCount = 100
        options.SaveDir
        options.SignificanceThreshold = 0.01
    end

    data = zeroStretchtotwopi(dihedralsMat);
    dihCount = size(data, 2);

    % Store all S1 values for quick access later
    s1 = arrayfun(@(dih) computeS1(data(:, dih), options.BinCount), 1:dihCount);

    % Scaling function for normalization, here log10()
    fn = @(x) log10(x);
    fnInv = @(x) 10 .^ x;

    % Create groups

    if isempty(options.DihTypes)
        groupLabels = "All pairs";
        pairsByGroup = {nchoosek(1:dihCount, 2)};
    else
        dihTypesUnique = sort(unique(options.DihTypes));

        % if ~all(dihTypesUnique == [0; 1])
        %     error("Invalid dihedral types")
        % end

        dihByType = cell(length(dihTypesUnique), 1);

        for typeIndex = 1:length(dihByType)
            dihByType{typeIndex} = find(options.DihTypes == dihTypesUnique(typeIndex));
        end

        typeIndicesByGroup = nchoosek(1:(length(dihTypesUnique) + 1), 2) - [0 1];
        groupLabels = strings(size(typeIndicesByGroup, 1), 1);
        pairsByGroup = cell(length(dihTypesUnique), 1);

        for groupIndex = 1:length(typeIndicesByGroup)
            groupTypeIndices = typeIndicesByGroup(groupIndex, :);
            groupTypes = dihTypesUnique(groupTypeIndices);
            groupLabels(groupIndex) = sprintf("Pairs %d-%d", groupTypes(:));

            if groupTypeIndices(1) == groupTypeIndices(2)
                pairsByGroup{groupIndex} = nchoosek(dihByType{groupTypeIndices(1)}, 2);
            else
                [ga, gb] = ndgrid(dihByType{groupTypeIndices(1)}, dihByType{groupTypeIndices(2)});
                pairsByGroup{groupIndex} = [ga(:) gb(:)];
            end
        end
    end

    excessMi = zeros(size(mi));

    for groupIndex = 1:length(pairsByGroup)
        groupPairs = pairsByGroup{groupIndex};
        miFlat = mi(sub2ind(size(mi), groupPairs(:, 1), groupPairs(:, 2)));

        % Calculate the histogram using the scaling function
        [bins, edges] = discretize(fn(miFlat), options.InterpBinCount);
        edges = fnInv(edges);

        figure;
        hold on;

        plot([1e-3 1e-1], [1e-3 1e-1]);

        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');

        title(groupLabels(groupIndex));
        xlabel("Mutual information");
        ylabel("Independent information");

        significanceArray = zeros(options.InterpBinCount,1);
        binMiArray = zeros(options.InterpBinCount,1);
        iiArray = zeros(options.InterpBinCount,1); 
        for binIndex = 1:options.InterpBinCount
            % Values and indices (with respect to miFlat) of all pairs contained in this bin
            binIndices = find(bins == binIndex);
            binMi = miFlat(binIndices);

            % Skip this bin if it is empty
            if length(binMi) < 1
                continue
            end

            % Find the values and indices (with respect to binMi and then miFlat) all the pairs at the center of the bin
            midBinIndices = randperm(length(binMi));
            midBinIndices = midBinIndices(1:min(length(midBinIndices), options.InterpPointCountPerBin));
            midIndices = binIndices(midBinIndices);
            midMi = binMi(midBinIndices);

            % Run the permutations
            permMi = zeros(length(midMi), options.PermutationCount);

            for midIndex = 1:length(midMi)
                pairx = groupPairs(midIndices(midIndex), 1);
                pairy = groupPairs(midIndices(midIndex), 2);

                x = data(:, pairx);
                y = data(:, pairy);

                s2 = zeros(options.PermutationCount, 1);

                for permutationIndex = 1:options.PermutationCount
                    s2(permutationIndex) = computeS2(x(randperm(length(x))), y, [options.BinCount, options.BinCount]);
                end

                permMi(midIndex, :) = s1(pairx) + s1(pairy) - s2;
            end

            % Compute the independent information for this bin
            ii = mean(permMi, 'all');

            % Test every pair using these permutations
            plotColors = zeros(length(binMi), 3);

            for pairBinIndex = 1:length(binMi)
                pairx = groupPairs(binIndices(pairBinIndex), 1);
                pairy = groupPairs(binIndices(pairBinIndex), 2);

                pvalue = sum(permMi(:) > binMi(pairBinIndex)) / numel(permMi);

                pairExcessMi = max(mi(pairx, pairy) - ii, 0);
                significant = pvalue <= options.SignificanceThreshold;

                excessMi(pairx, pairy) = pairExcessMi * significant;
                plotColors(pairBinIndex, 1) = significant;
            end

            plot(edges(binIndex:(binIndex + 1)), [ii ii], 'k');
            scatter(binMi, repmat(ii, length(binMi), 1), 3, plotColors);
            
            significanceArray(binIndex) = significant;
            binMiArray(binIndex) = mean(binMi);
            iiArray(binIndex) = ii;
        end

        if isfield(options, 'SaveDir')
            figPath = fullfile(options.SaveDir, sprintf("highpass_%d", groupIndex));
            % savefig(figPath); % Takes too much time
            print2pdf(figPath);
            save(fullfile(options.SaveDir, 'MI_filtered'),'excessMi');
        end
        % Add a legend in a sneaky way
        hold on
        h(1) = scatter(NaN,NaN,100,'black','filled','DisplayName',['p>' num2str(options.SignificanceThreshold)]);
        h(2) = scatter(NaN,NaN,100,'red','filled','DisplayName',['p<' num2str(options.SignificanceThreshold)]);
        
        legend(h);
        legend('boxoff')
    end

    excessMi = excessMi + excessMi';
end


% Utility functions

function output = binLimits(x)
    output = [max(min(x), 0), min(max(x), 2*pi)];
end

function output = computeS1(x, binCount)
    n = length(x);
    [hist, edges] = histcounts(x, binCount, 'BinLimits', binLimits(x), 'Normalization', 'pdf');
    width = edges(2) - edges(1);
    hist = hist(hist > 0);
    error = (length(hist) - 1) / (2 * n);
    output = error - width * sum(hist .* log(hist));
end

function output = computeS2(x, y, binCount)
    n = length(x);
    [hist, edgesx, edgesy] = histcounts2(x, y, binCount, 'XBinLimits', binLimits(x), 'YBinLimits', binLimits(y), 'Normalization', 'pdf');

    widthx = edgesx(2) - edgesx(1);
    widthy = edgesy(2) - edgesy(1);
    hist = hist(hist > 0);
    error = (length(hist) - 1) / (2 * n);
    output = error - widthx * widthy * sum(hist .* log(hist));
end
