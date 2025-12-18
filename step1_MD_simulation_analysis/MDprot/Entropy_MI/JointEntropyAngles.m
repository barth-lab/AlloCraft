% JointEntropyCont: Returns Joint entropy (in bits) of each column of 'X'
% with every column of Y. Assumes X and Y are continuous sets of data that
% will be binned in histograms
% by Mahdi Hijazi
%
%% Usage:
% H = JointEntropyAngles(X, Y)
% [H, histXY] = JointEntropyAngles(X, Y, nbins, binrangeX)
% [H, histXY] = JointEntropyAngles(X, Y, nbins, binrangeX, binrangeY)
%
%% Description:
% * H = Matrix of calculated entropies (in bits) of size ncolumns(X) x
% ncolumns(Y). Contains the joint entropies of every column of X with every
% column of Y
%
% * X, Y = data to be analyzed, in this case, X and Y are the
% angle data
%
% * nbins (optional) = [nbinsX nbinsY] number of bins for the histogram
%
% * binrangeX, binrangeY (optional) = range of data to be divided into
% bins, useful when dealing with angles or torsions, which we know will
% be limited into a range of [0, 2*pi] for example. When only binrangeX is
% given, binrangeY defaults to be equal to binrangeX
%
% * histXY (optional) = histogram used to calculate the entropies
%
% Note: A correction for understampling has been added for Eq 2 to 
% calculate entropy; otherwise, estimated entropy values are slightly less
% than true, due to finite sample size.
%
% Eq 2 based on : Niesen, M. J. et al (2013).
% The Journal of Physical Chemistry B, 117(24), 7283-7291.
%
% Last modified: Oct-29-2021

function [H, histXY] = JointEntropyAngles(X, Y, nbins, binrangeX, binrangeY)

if ~exist('nbins','var') || isempty(nbins)
    usenbins = 0; % Let Matlab figure out n of bins
else
    usenbins = 1; % User generated nbins
end

if ~exist('binrangeX','var') 
    usenbinrange = 0; % Let Matlab figure out the range
else
    usenbinrange = 1; % User generated range
end

if exist('binrangeX','var') & ~exist('binrangeY','var') 
    binrangeY = binrangeX; % set both X and Y ranges to be the same
end


    
% error('The crank''s length (%f) cannot exceed that of the slider (%f)', r2, r3)

% Size of data
[N,M1] = size(X);
[N2,M2] = size(Y);
assert( N2==N,'Data has different number of frames!!!')
H = zeros(M1,M2); % Contains joint pairwise entropy

for ColumnX = 1:M1
%     ColumnX
    if X == Y % If the data is identical, H2 matrix will be symmetrical
        jStart = ColumnX + 1;
    else
        jStart = 1;
    end
    for ColumnY = jStart:M2

        if usenbinrange % User generated range
            if usenbins % user specified nbins
                [histXY,xEdges,yEdges] = histcounts2(X(:,ColumnX),Y(:,ColumnY), ...
                    nbins, 'XBinLimits',binrangeX, 'YBinLimits',binrangeY);
                minX = binrangeX(1);
                minY = binrangeY(1);
            else % Let Matlab figure it out
                [histXY,xEdges,yEdges] = histcounts2(X(:,ColumnX),Y(:,ColumnY), ...
                    'XBinLimits',binrangeX, 'YBinLimits',binrangeY);
                nbins = size(histXY); 
                minX = binrangeX(1);
                minY = binrangeY(1);
            end
        else % Let Matlab figure out the range
            % OR add our own little trick:
            % X are the bonds
            maxX = max(X(:,ColumnX));
            minX = min(X(:,ColumnX));
            % Expand boundaries slightly
            maxX = maxX + 0.000000005;
            minX = minX - 0.000000005;
            maxX = min(maxX,pi); % Make sure the max and min are within [0 pi]
            minX = max(minX,0);
            
            % Y are the angles
            maxY = max(Y(:,ColumnY));
            minY = min(Y(:,ColumnY));
            % Expand boundaries slightly
            maxY = maxY + 0.000000005;
            minY = minY - 0.000000005;
            maxY = min(maxY,pi); % Make sure the max and min are within [0 pi]
            minY = max(minY,0);
            
            if usenbins % user specified nbins
                [histXY,xEdges,yEdges] = histcounts2(X(:,ColumnX),Y(:,ColumnY), ...
                    nbins, 'XBinLimits', [minX maxX], ...
                    'YBinLimits',[minY maxY]); 
            else % Let Matlab figure it out
                [histXY,xEdges,yEdges] = histcounts2(X(:,ColumnX),Y(:,ColumnY), ...
                    'XBinLimits', [minX maxX], ...
                    'YBinLimits',[minY maxY]);  
                nbins = size(histXY); 
            end
        end

        xEdgew = xEdges(2)-xEdges(1); % edge width, used in entropy calculation
        yEdgew = yEdges(2)-yEdges(1);
        
        % Calculate sample class probabilities
        P = histXY / N;

        % Calculate entropy:

        E = (sum(sum(histXY>0))-1)/(2*N); % error for undersampling (Only nonzero
        %  bins contribute to the error)
        % (R.Steuer et al. (2002 Bioinf.), B.Killian et al. (2012 JCTC))

        % When calculating H, don't forget to account for the Jacobian term for
        % the transformation from cartesian to BAT coordinates
        thi = minX + xEdgew/2;
        for binX = 1:nbins(1) % angle part i
            thj = minY + yEdgew/2;
            for binY = 1:nbins(2) % angle part j
                if P(binX,binY) > 0 % Non zero prob, add it to sum
                H(ColumnX,ColumnY) = H(ColumnX,ColumnY) - P(binX,binY) .* log(P(binX,binY)/(sin(thi)*sin(thj)*xEdgew*yEdgew));
                end
                thj = thj + yEdgew;
            end
            thi = thi +  xEdgew;
        end

        % Correct for undersampling:
        H(ColumnX,ColumnY) = H(ColumnX,ColumnY) + E;
    end
end


