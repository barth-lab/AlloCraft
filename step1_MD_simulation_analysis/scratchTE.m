% Transfer entropy calculation concepts:
% Kamberaj and van der Vaart, Biophys J. 2009
% Schreiber, PRL 2000

% Time delayed embedding of time series:

% Time series could either be fluctuations r(tn) = abs(xi(tn) - mean(xi)),
% where i is the atom and tn is the time at which 
% 
% I could also construct a time series of dihedrals for my analysis

% Embedding parameters: m and tau
% tau is a multiple of 

% Input 
trj = mainSim.traj{1};
index = 1; % Atom/residue/dihedlal index
tn = 100; % Time step at which r is evaluated

rtn = norm(trj(tn,to3(index)) - mean(trj(:,to3(index))));


% Embedding parameters
m = 10; % Length of embedding vector 
tau = 5; % Time lag between samples 

mTau = 1:tau:m*tau;
Iembed = arrayfun(@(x) computeRtn(trj,x,index),mTau);



function rtn = computeRtn(trj,tn,index) % Compute positional fluctuations
    rtn = norm(trj(tn,to3(index)) - mean(trj(:,to3(index))));
end