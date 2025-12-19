function X = zeroStretchtotwopi(X, modfitbins)
%zeroStretchto2pi moves the zero stretch in dihedral data toward the 2pi
%end. This helps make peaks clearer.
%% Usage:
% [X] = zeroStretchto2pi(X)
% [X] = zeroStretchto2pi(X, modfitbins)

if ~exist('modfitbins','var') || isempty(modfitbins)
   modfitbins = 100; % Set default value for modfit
end
    
% Size of data
[N,M] = size(X);

for Column = 1:M % Loop over all the dihedrals in question
    histo = zeros(modfitbins,1); % initialize the histogram
    
    binSize =(2*pi+0.000000005)/modfitbins;
    
    for frame = 1:N % build the histogram
        thisbin = floor(X(frame,Column)/binSize) + 1;
        histo(thisbin) = histo(thisbin) + 1;
    end
    zeroesPos = find(histo==0); % find if there are ay empty bins
    if zeroesPos % if any of the bins of the histogram is empty find the 
        %longest consecutive stretch of  empty bins
        longestZeroStretch=0;
        currentZeroStretch=0;
        longestZeroStretchPos=-1;
        for k = 1:2*modfitbins
           l = mod(k,modfitbins) + 1; % we go over the data twice to consider stretches that loop over 2pi
            if (currentZeroStretch==0) && (histo(l)==0)
                currentZeroStretch=1;
                currentZeroStretchPos=k+1;
            end 
            
            if (currentZeroStretch>0) && (histo(l)==0)
                currentZeroStretch = currentZeroStretch + 1;
            end
            if (currentZeroStretch>0) && (histo(l)~=0)
                if (currentZeroStretch>longestZeroStretch)
                    longestZeroStretch = currentZeroStretch;
                    longestZeroStretchPos = currentZeroStretchPos;
                end
                currentZeroStretch = 0;
            end
            
        end
        
    else % If none of the bins is empty
        % misuse the zeroStretch variables for determining the minimum
        [longestZeroStretch, longestZeroStretchPos] = min(histo);
    end
    modFit = 2*pi - (longestZeroStretchPos+0.5)*binSize; % calc shift to throw zero stretch toward 2pi end
    for frame = 1:N % Make the shift
        X(frame,Column) = X(frame,Column) + modFit - 2*pi*floor((X(frame,Column) + modFit)/(2*pi));
    end
end

end


