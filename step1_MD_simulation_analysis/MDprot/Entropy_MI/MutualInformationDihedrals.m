function I = MutualInformationDihedrals(X,Y, nbins, binrangeX, binrangeY)
%MutualInformation Calcaultes the mutual information between sets X and Y
% by Mahdi Hijazi
%% Usage:
% I = MutualInformationDihedrals(X, Y)
% I = MutualInformationDihedrals(X, Y, nbins, binrangeX, binrangeY)
%
%% Description:
% * H = Matrix of calculated entropies (in bits) of size ncolumns(X) x
% ncolumns(Y). Contains the joint entropies of every column of X with every
% column of Y
%
% * X, Y = data to be analyzed
%
% * nbins (optional) = [nbinsX nbinsY] number of bins for the histogram
%
% * binrangeX, binrangeY (optional) = range of data to be divided into
% bins, useful when dealing with angles or torsions, which we know will
% be limited into a range of [0, 2*pi] for example. When only binrangeX is
% given, binrangeY defaults to be equal to binrangeX

% Initialize deafults for the variables
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

% Size of data
[~,M1] = size(X);
[~,M2] = size(Y);

if usenbinrange % User generated range
    if usenbins % user specified nbins
        HX = EntropyDihedrals(X,nbins(1),binrangeX);
        HY = EntropyDihedrals(Y,nbins(2),binrangeY);
        HXY = JointEntropyDihedrals(X,Y, nbins, binrangeX, binrangeY);
    else % Let Matlab figure it out
        HX = EntropyDihedrals(X,[],binrangeX);
        HY = EntropyDihedrals(Y,[],binrangeY);
        HXY = JointEntropyDihedrals(X,Y, [], binrangeX, binrangeY);
    end
else % Let Matlab figure out the range
    if usenbins % user specified nbins
        HX = EntropyDihedrals(X,nbins(1));
        HY = EntropyDihedrals(Y,nbins(2));
        HXY = JointEntropyDihedrals(X,Y, nbins); 
    else % Let Matlab figure it out
        HX = EntropyDihedrals(X);
        HY = EntropyDihedrals(Y);
        HXY = JointEntropyDihedrals(X,Y);
    end
end
        
HXmat = repmat(HX',1,M2); % X differs in rows, same in columns
HYmat = repmat(HY,M1,1);  % Y differs in columns, same in rows

% Axiom of information theory:
I = HXmat + HYmat - HXY;
end

