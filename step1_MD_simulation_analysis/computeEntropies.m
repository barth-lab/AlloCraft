function [s1, s2, i2] = computeEntropies(subset, options)
    arguments
        subset
        options.AdaptivePartitioning = false
        options.BinCount = 50
    end

    dihCount = size(subset, 2);

    if options.AdaptivePartitioning
        [~, order] = sort(subset);

        for dih = 1:dihCount
            subset(order(:, dih), dih) = 1:length(subset);
        end
    end

    s1 = arrayfun(@(dih) computeS1(subset(:, dih), options.BinCount), 1:dihCount);

    if nargout > 1 % Don't calculate what's not asked
        s2 = zeros(dihCount, dihCount);
        i2 = zeros(dihCount, dihCount);
        tic
            parfor dihx = 1:dihCount
                x = subset(:, dihx);
    
                for dihy = 1:dihCount
                    if dihy >= dihx
                        break;
                    end
    
                    y = subset(:, dihy);
    
                    s2(dihx, dihy) = computeS2(x, y, [options.BinCount options.BinCount]);
                    i2(dihx, dihy) = s1(dihx) + s1(dihy) - s2(dihx, dihy);
                end
            end
        toc
    end
end


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
