% Entropy: Returns entropy (in bits) of each column of 'X', assumes X is a
% continuous set of data that will be binned in a histogram
% by Mahdi Hijazi
%
%% Usage:
% H = EntropyBonds(X)
% [H, histX] = EntropyBonds(X, nbins, binrange)
%
%% Description:
% * H = row vector of calculated entropies (in bits)
%
% * X = data to be analyzed
%
% * nbins (optional) = number of bins for the histogram
%
% * binrange (optional) = range of data to be divided into bins, useful when
% dealing with angles or torsions, which we know will be limited into a
% range of [0, 2*pi] for example.
%
% * histX (optional) = histogram used to calculate the entropies
%
% Note: A correction for understampling has been added to calculate
% entropy; otherwise, estimated entropy values are slightly less
% than true, due to finite sample size.
%
% Correction based on : Niesen, M. J. et al (2013).
% The Journal of Physical Chemistry B, 117(24), 7283-7291.
%
% Last modified: Oct-29-2021

function [H, histX] = EntropyBonds(X, nbins, binrange)

if ~exist('nbins','var') || isempty(nbins)
    usenbins = 0; % Let Matlab figure out n of bins
else
    usenbins = 1; % User generated nbins
end

if ~exist('binrange','var') 
    usenbinrange = 0; % Let Matlab figure out the range
else
    usenbinrange = 1; % User generated range
end

% Size of data
[N,M] = size(X);
H = zeros(1,M);

for Column = 1:M
    
    if usenbinrange % User generated range
        if usenbins % user specified nbins
            [histX,edges]  = histcounts(X(:,Column), nbins, 'BinLimits',binrange); 
        else % Let Matlab figure it out
            [histX,edges]  = histcounts(X(:,Column),'BinLimits',binrange); 
            nbins = length(histX); 
        end
    else % Let Matlab figure out the range
        % OR add our own little trick:
        maxX = max(X(:,Column));
        minX = min(X(:,Column));
        assert( minX >= 0, ' Negative bonds exist, check data!!!')
        % Expand boundaries slightly
        maxX = maxX + 0.000000005;
        minX = minX - 0.000000005;
        minX = max(minX,0);
        
        % If distribution is too narrow, widen it:
        if (maxX - minX) <1e-4
            maxX =  maxX + 0.05;
        end
        if usenbins % user specified nbins
            [histX,edges]  = histcounts(X(:,Column), nbins, 'BinLimits', ...
                [minX maxX]); 
        else % Let Matlab figure it out
            [histX,edges]  = histcounts(X(:,Column), 'BinLimits', ...
                [minX maxX]);  
            nbins = length(histX); 
        end
    end
    
 	edgew = edges(2)-edges(1); % edge width, used in entropy calculation
    
    % Calculate sample class probabilities
    P = histX / N;
	
    % Calculate entropy in bits
    % Eq 2 for entropy: using this equation directly produces NaNs for
    % empty bins in the histogram
%    H(Column) = -sum(P .* log(P/edgew)) + (nbins-1)/(2*n);

    E = (sum(histX>0)-1)/(2*N); % error for undersampling (Only nonzero
    %  bins contribute to the error)
    % (R.Steuer et al. (2002 Bioinf.), B.Killian et al. (2012 JCTC))
    
    % When calculating H, don't forget to account for the Jacobian term for
    % the transformation from cartesian to BAT coordinates
    bi = minX + edgew/2;
    for bin = 1:nbins
        if P(bin) > 0 % Non zero prob, add it to sum
        H(Column) = H(Column) - P(bin) .* log(P(bin)/(bi*bi*edgew));
        end
        bi = bi + edgew;
    end
    % Correct for undersampling:
    H(Column) = H(Column) + E;
end
