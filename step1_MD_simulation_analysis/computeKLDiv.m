function [kl1, kl2, m2, pVal, pVal2, kl2H0 ] = computeKLDiv(dihedralsRef, dihedralsTest, nBlocks, options)
%computeKLDiv: Computes marginal and joint Kullback-Leibler (KL) divergences  
% as well as mutual divergence of a test sample and reference data. Assumes
% input are dihedrals
%% Usage:
% [kl1] = computeKLDiv(dihedralsRef,dihedralsTest, nBlocks)
% [kl1, kl2, m2] = computeKLDiv(dihedralsRef,dihedralsTest, nBlocks)
% [kl1, kl2, m2] = computeKLDiv(dihedralsRef,dihedralsTest, nBlocks, options)
%% Description:
% * kl1: marginal Kullback-Leibler divergence for one degree of freedom, if
% only this output is given, then  kl2 and m2 will not be calculated and
% the function will run much faster
%
% * kl2: Kullback-Leibler divergence for a joint distribution of two
% degrees of freedom
%
% * m2: mutual divergence, the mutual information analog for KL div. This
% term could be understood as: if you know how divergent d.o.f. i is in
%  test ensemble from reference ensemble, what does this tell us about
%  divergence of d.o.f. j between test and reference ensembles
%
% * dihedralsRef,dihedralsTest: input dihedral time series from molecular
% simulations. [N x Ndih], where N is the number of frames/observations and
% Ndih is the number of dihedrals. Ndih must be the same in both inputs
%
% * nBlocks: number of blocks to split the reference dihedral time series
% into for significance testing using the bootstrapping method in McClendon  
% et al
%
% * [kl1, kl2, m2] = computeKLDiv(..., 'BinCount', nBins)
%   Sets the number of bins for histogramming. 2D histograms will use the
%   same number of bins in both dimensions
%
% * [kl1, kl2, m2] = computeKLDiv(..., 'Grassberger', true)
%    Correct for finite sample size in kl1 with the Grassberger 
%    (1988 Phys Lett A) estimate, defaults to FALSE
%
% * [kl1, kl2, m2] = computeKLDiv(..., 'SignificanceThreshold', 0.1)
%   Specifies threshold for considering KL div of a d.o.f. or
%   a pair of d.o.f.s significant, defaults to 0.05
%
% * [kl1, kl2, m2] = computeKLDiv(..., 'quick2ndOrderCalc', false)
%   Calculates the full matrix of 2nd order KLDiv terms. This option
%   defaults to true, which ignores 2nd order terms when both 1st order
%   d.o.f.s are insignificant and thus makes the calculation faster
%
% The code is based on the derivation of McClendon et al., J. Chem. Theory
% Comput. (2012)
    arguments
        dihedralsRef
        dihedralsTest
        nBlocks
        options.AdaptivePartitioning = false % Correct this for KL
        options.BinCount = 50
        options.Grassberger = false % Histogram finite size effect correction
        options.SignificanceThreshold = 0.05
        options.quick2ndOrderCalc = true
    end


%% Main function body:

    dihCount = size(dihedralsRef, 2);
    kl2 = zeros(dihCount, dihCount);
    m2 = zeros(dihCount, dihCount);

    if options.AdaptivePartitioning % Correct this for KL
        [~, order] = sort(dihedralsRef);

        for dih = 1:dihCount
            dihedralsRef(order(:, dih), dih) = 1:length(dihedralsRef);
        end
    end
    
     % If nBlocks is not even, give an error!
     assert( mod(nBlocks,2) == 0,'Error: nBlocks must be an even number!!!')   
     
    kl1 = arrayfun(@(dih) computeKL1(dihedralsRef(:, dih), dihedralsTest(:, dih), ... 
         options.BinCount,'Grassberger', options.Grassberger), 1:dihCount);

    % Statistical significance:
    % Calculate KL div in the reference ensemble to quantify intra-ensemble
    % variability 
    interface = floor(linspace(0, (size(dihedralsRef,1) - mod(size(dihedralsRef,1),nBlocks)), nBlocks+1));
    
    % Make temp dihedrals variable divided by blocks
    dihedralsTemp = zeros((size(dihedralsRef,1) - mod(size(dihedralsRef,1),nBlocks))/nBlocks, dihCount, nBlocks);
    for iblock = 1:nBlocks
        istart = interface(iblock)+1;
        iend = interface(iblock+1);
        dihedralsTemp(:,:,iblock) = dihedralsRef(istart:iend,:);
    end
    
    combinations =  nchoosek(1:nBlocks,nBlocks/2); 
    % Do we repeat 1 2 with 3 4 and then 3 4 with 1 2? YES! Since KLDiv is not
    % symmetric
    combinations = [combinations flipud(combinations)];
    nCombs = nchoosek(nBlocks,nBlocks/2);
    kl1H0 = zeros(nCombs,dihCount);
    for iblock = 1:nCombs % Combinations of blocks
    
        kl1H0(iblock,:) = arrayfun(@(dih) computeKL1(reshape(dihedralsTemp(:,dih,combinations(iblock,1:nBlocks/2)),[],1,1) ...
            , reshape(dihedralsTemp(:,dih,combinations(iblock,(nBlocks/2+1):end)),[],1,1), ... 
                 options.BinCount,'Grassberger', options.Grassberger), 1:dihCount);
    end
    kl1H0Mean = sum(kl1H0)/nCombs; % Average intra-ensemble variability
    
    % Now get the p-value for every d.o.f., zero out non-significant KL divs,
    % and subtract "excess" divergence from significant ones
    pVal = sum((kl1 - kl1H0)<0) / nCombs;
    
    kl1 = kl1 - kl1H0Mean;
    kl1(pVal >= options.SignificanceThreshold) = 0;

    % 2nd order domain: do I need to calcualte 2nd order if both 1st order
    % dihedrals show insignificant KLDiv?
    
    if nargout > 1 % Don't calculate what's not asked
        tic
        kl2H0 = zeros(dihCount, dihCount,nCombs);
        pVal2 = ones(dihCount, dihCount);

        parfor dihx = 1:dihCount
            xRef = dihedralsRef(:, dihx);
            x = dihedralsTest(:,dihx);

            for dihy = 1:dihCount
                if dihy >= dihx
                    break;
                end

                if options.quick2ndOrderCalc % Skip if both dofs are insignificant
                   if pVal(dihx) >= options.SignificanceThreshold && pVal(dihy) >= options.SignificanceThreshold
                        continue;
                   end
                end

                yRef = dihedralsRef(:, dihy);
                y = dihedralsTest(:,dihy);

                kl2(dihx, dihy) = computeKL2(xRef, x, yRef, y, [options.BinCount options.BinCount],'Grassberger', options.Grassberger);

                % Statistical significance testing:
                for iblock = 1:nCombs% Combinations of blocks
    
                kl2H0(dihx, dihy,iblock) = computeKL2(reshape(dihedralsTemp(:,dihx,combinations(iblock,1:nBlocks/2)),[],1,1), ...
                    reshape(dihedralsTemp(:,dihx,combinations(iblock,(nBlocks/2+1):end)),[],1,1), ...
                     reshape(dihedralsTemp(:,dihy,combinations(iblock,1:nBlocks/2)),[],1,1), ...
                     reshape(dihedralsTemp(:,dihy,combinations(iblock,(nBlocks/2+1):end)),[],1,1), ...
                     [options.BinCount options.BinCount],'Grassberger', options.Grassberger);
                end
            end
        end
        kl2H0Mean = sum(kl2H0,3)/nCombs; % Average intra-ensemble variability

        % Calculate p-value for 2nd order shit 
        parfor dihx = 1:dihCount
            for dihy = 1:dihCount
                if dihy >= dihx
                    break;
                end
                if options.quick2ndOrderCalc % Skip if both dofs are insignificant
                   if pVal(dihx) >= options.SignificanceThreshold && pVal(dihy) >= options.SignificanceThreshold
                        continue;
                   end
                end
                pVal2(dihx, dihy) = sum((kl2(dihx, dihy) - kl2H0(dihx, dihy,:))<0) / nCombs;
            end
        end
%         pVal2 = sum((kl2 - kl2H0)<0,3)/ nCombs;
    
        kl2 = kl2 - kl2H0Mean;
        kl2(pVal2 >= options.SignificanceThreshold) = 0;

        % Mutual divergence:
        m2 = tril(repmat(kl1,dihCount,1),-1) + tril(repmat(kl1',1,dihCount),-1) - kl2;
        toc
    end
end

%% Support functions:

function output = binLimits(xref,x) 
    % Creates common bin limits for two dihderal series while keeping them
    % in the [ 0 2*pi] limit
    minCommon = min(min(x), min(xref));
    maxCommon = max(max(x), max(xref));
    output = [max(minCommon, 0), min(maxCommon, 2*pi)];

end

function output = computeKL1(xRef, x, binCount, options) 
    arguments
        xRef
        x
        binCount
        options.Grassberger = false
    end

    if options.Grassberger == false
        normalization = 'pdf';
    else
        normalization = 'count'; 
    end
    [hist, edges] = histcounts(x, binCount, 'BinLimits', binLimits(xRef,x), 'Normalization', normalization);
    [histRef, ~] = histcounts(xRef, binCount, 'BinLimits', binLimits(xRef,x), 'Normalization', normalization);
    edgeWidth = edges(2) - edges(1);

    % Remove bins with zero probability from ref. state
    hist = hist(histRef > 0);
    histRef = histRef(histRef > 0);

    % Remove bins with zero probability from test state
    histRef = histRef(hist > 0);
    hist = hist(hist>0);

    if options.Grassberger == false
        output = edgeWidth * sum(hist .* log(hist./histRef));
    else % Calculation of KL1 with Grassberger estimate
        output = 1/length(xRef)*sum(hist.*(psi(hist) - psi(histRef)-1./histRef));
    end
end


function output = computeKL2(xRef, x, yRef, y, binCount, options) 
    arguments
        xRef
        x
        yRef
        y
        binCount
        options.Grassberger = false
    end

    if options.Grassberger == false
        normalization = 'pdf';
    else
        normalization = 'count'; 
    end
    [hist, edgesx, edgesy] = histcounts2(x, y, binCount, 'XBinLimits', binLimits(xRef,x), 'YBinLimits', binLimits(yRef,y), 'Normalization', normalization);
    [histRef] = histcounts2(xRef, yRef, binCount, 'XBinLimits', binLimits(xRef,x), 'YBinLimits', binLimits(yRef,y), 'Normalization', normalization);

    widthx = edgesx(2) - edgesx(1);
    widthy = edgesy(2) - edgesy(1);

    % Remove bins with zero probability from ref. state
    hist = hist(histRef > 0);
    histRef = histRef(histRef > 0);

    % Remove bins with zero probability from test state
    histRef = histRef(hist > 0);
    hist = hist(hist>0);

    if options.Grassberger == false
        output =   widthx * widthy * sum(hist .* log(hist./histRef));
    else % Calculation of KL2 with Grassberger estimate
        output = 1/length(xRef)*sum(hist.*(psi(hist) - psi(histRef)-1./histRef));        
    end
end
