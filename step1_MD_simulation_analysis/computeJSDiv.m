function [JS1, pVal ] = computeJSDiv(dihedralsRef, dihedralsTest, nBlocks, options)
%computeKLDiv: Computes marginal Jensen−Shannon Divergence of a test sample 
% and reference data. JS divergence is a symmetric form of KL divergence
% where the reference distribution is a superposition between test and
% reference. 
% Assumes input are dihedrals. (For now)
%
%% Usage:
% [JS1] = computeKLDiv(dihedralsRef,dihedralsTest, nBlocks)
% [JS1, pVal ] = computeJSDiv(dihedralsRef, dihedralsTest, nBlocks, options)
%% Description:
% * JS1: marginal Jensen−Shannon divergence for one degree of freedom
%
% * dihedralsRef,dihedralsTest: input dihedral time series from molecular
% simulations. [N x Ndih], where N is the number of frames/observations and
% Ndih is the number of dihedrals. Ndih must be the same in both inputs
%
% * nBlocks: number of blocks to split the reference dihedral time series
% into for significance testing using the bootstrapping method in McClendon  
% et al
%
% * [kl1, pVal] = computeKLDiv(..., 'BinCount', nBins)
%   Sets the number of bins for histogramming. 2D histograms will use the
%   same number of bins in both dimensions
%
% * [kl1, pVal] = computeKLDiv(..., 'Grassberger', true)
%    Correct for finite sample size in kl1 with the Grassberger 
%    (1988 Phys Lett A) estimate, defaults to FALSE
%
% * [kl1, pVal] = computeKLDiv(..., 'SignificanceThreshold', 0.1)
%   Specifies threshold for considering KL div of a d.o.f. or
%   a pair of d.o.f.s significant, defaults to 0.05
%
% The code is based on the derivation of McClendon et al., J. Chem. Theory
% Comput. (2012)
    arguments
        dihedralsRef
        dihedralsTest
        nBlocks
        options.BinCount = 50
        options.Grassberger = false % Histogram finite size effect correction
        options.SignificanceThreshold = 0.05
    end


%% Main function body:

    dihCount = size(dihedralsRef, 2);

     % If nBlocks is not even, give an error!
     assert( mod(nBlocks,2) == 0,'Error: nBlocks must be an even number!!!')   
     
    JSp = arrayfun(@(dih) computeKL1([dihedralsRef(:, dih) ;dihedralsTest(:, dih)], ...
        dihedralsTest(:, dih), options.BinCount,'Grassberger',... 
        options.Grassberger), 1:dihCount);
    JSq = arrayfun(@(dih) computeKL1([dihedralsRef(:, dih) ;dihedralsTest(:, dih)], ...
        dihedralsRef(:, dih), options.BinCount,'Grassberger',... 
        options.Grassberger), 1:dihCount);
    JS1 = 0.5*(JSp + JSq);

    % Statistical significance:
    % Calculate KL div in the reference and target ensembles to have the
    % null hypothesis JS
    interfaceRef = floor(linspace(0, (size(dihedralsRef,1) - mod(size(dihedralsRef,1),nBlocks)), nBlocks+1));
    interfaceTest = floor(linspace(0, (size(dihedralsTest,1) - mod(size(dihedralsTest,1),nBlocks)), nBlocks+1));

    % Make temp dihedrals variable divided by blocks
    dihedralsTemp = zeros((size(dihedralsRef,1) - mod(size(dihedralsRef,1),nBlocks))/nBlocks, dihCount, nBlocks);
    dihedralsTempTest = zeros((size(dihedralsTest,1) - mod(size(dihedralsTest,1),nBlocks))/nBlocks, dihCount, nBlocks);
    for iblock = 1:nBlocks
        istart = interfaceRef(iblock)+1;
        iend = interfaceRef(iblock+1);
        dihedralsTemp(:,:,iblock) = dihedralsRef(istart:iend,:);

        istart = interfaceTest(iblock)+1;
        iend = interfaceTest(iblock+1);
        dihedralsTempTest(:,:,iblock) = dihedralsTest(istart:iend,:);
    end
    
    combinations =  nchoosek(1:nBlocks,nBlocks/2); 
    % Do we repeat 1 2 with 3 4 and then 3 4 with 1 2? YES! Since KLDiv is not
    % symmetric
    combinations = [combinations flipud(combinations)];
    nCombs = nchoosek(nBlocks,nBlocks/2);
%     kl1qH0 = zeros(nCombs,dihCount);
%     kl1pH0 = zeros(nCombs,dihCount);
    JS1H0 = zeros(nCombs,dihCount);
    for iblock = 1:nCombs % Combinations of blocks
   
        kl1qH0 = arrayfun(@(dih) computeKL1(reshape(dihedralsTemp(:,dih,combinations(iblock,1:nBlocks/2)),[],1,1) ...
            , reshape(dihedralsTemp(:,dih,combinations(iblock,(nBlocks/2+1):end)),[],1,1), ... 
                 options.BinCount,'Grassberger', options.Grassberger), 1:dihCount);

        kl1pH0 = arrayfun(@(dih) computeKL1(reshape(dihedralsTempTest(:,dih,combinations(iblock,1:nBlocks/2)),[],1,1) ...
            , reshape(dihedralsTempTest(:,dih,combinations(iblock,(nBlocks/2+1):end)),[],1,1), ... 
                 options.BinCount,'Grassberger', options.Grassberger), 1:dihCount);
        JS1H0(iblock,:) = 0.5*(kl1qH0 + kl1pH0);
    end
    JS1H0Mean = sum(JS1H0)/nCombs; % Average intra-ensemble variability
    
    % Now get the p-value for every d.o.f., zero out non-significant KL divs,
    % and subtract "excess" divergence from significant ones
    pVal = sum((JS1 - JS1H0)<0) / nCombs;
    
    JS1 = JS1 - JS1H0Mean;
    JS1(pVal >= options.SignificanceThreshold) = 0;

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