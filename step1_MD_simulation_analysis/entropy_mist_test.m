%% MIE/MIST computatiom

f = 100; % test frames
b = 50; % frame per test frame

m = 10; % degrees of freedom
n = f * b; % time frames
binCount = 100;

cov = unifrnd(1, 2, m, m);
cov = cov * cov';
% cov = [4 0; 0 4];

data = mvnrnd(zeros(1, m), cov, n);

s1r = computeS1Normal(sqrt(cov));
snr = computeSnNormal(cov);

s1 = arrayfun(@(x) computeS1(data(:, x), binCount), 1:m);

% display(s1r);
% display(s1);
display(sum(abs(s1r - s1)));


%


output = zeros(f, 4);

for frame = 1:f
    disp(frame);

    subset = data(1:(frame * b), :);
    s1 = arrayfun(@(x) computeS1(subset(:, x), binCount), 1:m);

    s2 = zeros(m, m);
    i2 = zeros(m, m);

    for i = 1:m
        for j = 1:m
            s2(i, j) = computeS2(subset(:, i), subset(:, j), [binCount binCount]);
            i2(i, j) = s1(i) + s1(j) - s2(i, j);
        end
    end

    % Fix imprecision
    i2 = 0.5 * (i2 + i2');
    g = graph(-i2, 'omitselfloops'); % We need max spanning tree so we take negative edge weights
    [T, pred] = minspantree(g);

    output(frame, 1) = frame * b;
    output(frame, 2) = sum(s1) - sum(tril(i2, -1), 'all');
    output(frame, 3) = sum(s1) - sum(arrayfun(@(x, i) i2(x, i), pred(2:end), 1:(m - 1)));
    output(frame, 4) = sum(s1);
end


%%

figure;
hold on;

yline(snr);
plot(output(:, 1), output(:, 2), '-o');
plot(output(:, 1), output(:, 3), '-o');
plot(output(:, 1), output(:, 4), '-o');

legend("Exact S_n", "S_n MIE", "S_n MIST", "\Sigma S_1");

% display(snr);
% display(sna);
% display(snb);

display(snr);
display(output(end, 2));
display(output(end, 3));


%%

figure;
plot(data(:, 1), data(:, 2), '.');


%%

figure;

p = plot(g, 'EdgeLabel', g.Edges.Weight);


highlight(p, T);



%% Functions

function rank = sortRank(data)
    [~, order] = sort(data);
    rank = 1:length(data);
    rank(order) = rank;
    rank = rank';
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
    % figure;
    % hist3([x y]);

    n = length(x);
    [hist, edgesx, edgesy] = histcounts2(x, y, binCount, 'XBinLimits', binLimits(x), 'YBinLimits', binLimits(y), 'Normalization', 'pdf');

    widthx = edgesx(2) - edgesx(1);
    widthy = edgesy(2) - edgesy(1);
    hist = hist(hist > 0);
    error = (length(hist) - 1) / (2 * n);
    output = error - widthx * widthy * sum(hist .* log(hist));
end

function output = computeI2(x, y, binCount)
    output = computeS1(x, binCount(1)) + computeS1(y, binCount(2)) - computeS2(x, y, binCount);
end


function output = computeS1p(x, binCount)
    n = length(x);
    [hist, edges] = histcounts(x, binCount);
    hist = hist(hist > 0);
    output = sum(hist / n .* (log(n) - psi(hist) - (-1).^hist ./ (hist + 1)));
end

function output = computeS2p(x, y, binCount)
    n = length(x);
    [hist, edgesx, edgesy] = histcounts2(x, y, binCount);
    hist = hist(hist > 0);
    output = sum(hist / n .* (log(n) - psi(hist) - (-1).^hist ./ (hist + 1)));
end

function output = computeI2p(x, y, binCount)
    output = computeS1p(x, binCount(1)) + computeS1p(y, binCount(2)) - computeS2p(x, y, binCount);
end



function [s1, s2, i2, pv] = computeEntropies(subset, options)
    arguments
        subset
        options.AdaptivePartitioning = false
        options.BinCount = 10
        options.ComputeIndependentInformation = false
    end

    dihCount = size(subset, 2);
    s2 = zeros(dihCount, dihCount);
    i2 = zeros(dihCount, dihCount);
    pv = zeros(dihCount, dihCount);

    if options.AdaptivePartitioning
        [~, order] = sort(subset);

        for dih = 1:dihCount
            subset(order(:, dih), dih) = 1:length(subset);
        end
    end

    s1 = arrayfun(@(dih) computeS1(subset(:, dih), options.BinCount), 1:dihCount);


    % Compute independent information

    if options.ComputeIndependentInformation
        % indepDih = randi([1 dihCount], 1000, 2);
        % indepInfos = zeros(1, length(indepDih));

        % tic
        %     parfor indepIndex = 1:length(indepDih)
        %         dihx = indepDih(indepIndex, 1);
        %         dihy = indepDih(indepIndex, 2);

        %         x = subset(:, dihx);
        %         y = subset(:, dihy);

        %         indepInfos(indepIndex) = s1(dihx) + s1(dihy) - computeS2(x(randperm(length(x))), y, [options.BinCount options.BinCount]);
        %     end
        % toc

        % ii = mean(indepInfos);

        % permCount = 10;
        % ii = mean(arrayfun( ...
        %     @(dihx, dihy) s1(dihx) + s1(dihy) - computeS2(subset(randperm(length(subset)), dihx), subset(:, dihy), [options.BinCount options.BinCount]), ...
        %     repmat((1:dihCount)', 1, permCount), ...
        %     randi([1 dihCount], dihCount, permCount) ...
        % ), 2);
    end


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
                pv(dihx, dihy) = testS2(x, y, s2(dihx, dihy), options.BinCount);
            end
        end
    toc
end

function output = testI2s(x, y, binCount)
    permCount = 100;
    s2 = arrayfun(@(i) computeS2(x(randperm(length(x))), y, [binCount, binCount]), 1:permCount);
    output = computeS1(x, binCount) + computeS1(y, binCount) - s2;
end

function output = testI2(x, y, binCount)
    permCount = 100;
    s2 = zeros(permCount, 1);

    for j = 1:permCount
        s2(j) = computeS2(x(randperm(length(x))), y, [binCount, binCount]);
    end

    output = computeS1(x, binCount) + computeS1(y, binCount) - s2;
end

function output = testI2p(x, y, binCount)
    permCount = 100;
    s2 = arrayfun(@(i) computeS2p(x(randperm(length(x))), y, [binCount, binCount]), 1:permCount);
    output = computeS1p(x, binCount) + computeS1p(y, binCount) - s2;
end


function output = computeS1Normal(cov)
    output = 0.5 * log(2 * pi * diag(cov)'.^2) + 0.5;
end

function output = computeSnNormal(cov)
    output = 0.5 * length(cov) * (1 + log(2 * pi)) + 0.5 * log(abs(det(cov)));
end


% y: array of size [1 x m]
%   - sorted with respect to the x axis
% outputIndices: array of size [1 x n] where n <= options.BinCount * options.SampleCountPerBin and n <= m
function outputIndices = sampleUniformly(y, options)
    arguments
        y
        options.BinCount = 50
        options.SampleCountPerBin = 20
    end

    [binIndices, edges] = discretize(y, options.BinCount);
    outputIndices = zeros(1, 0);

    for i = 1:(length(edges) - 1)
        indices = find(binIndices == i);

        if ~isempty(indices)
            % outputIndices = [outputIndices round(linspace(min(indices), max(indices), options.SampleCountPerBin))];

            addIndices = linspace(min(indices), max(indices) + 1, options.SampleCountPerBin + 1);
            outputIndices = [outputIndices round(addIndices)];
            % outputIndices = [outputIndices round(addIndices(1:(end - 1)))];
        end
    end

    % outputIndices = unique(outputIndices)';
    outputIndices = unique(min(outputIndices, length(y)))';
end