function [sysState,index] = calcState(dihedrals,dihedrals_active,dihedrals_inactive,tol)
%calcState: Checks whether the dihedrals are in an active or inactive
%conformation
%
%% Usage:
% [sysState] = calcState(dihedrals,dihedrals_active,dihedrals_inactive)
% [sysState] = calcState(dihedrals,dihedrals_active,dihedrals_inactive, tol)
% [sysState,index] = calcState(dihedrals,dihedrals_active,dihedrals_inactive)
%
%% Description:
% * dihedrals: a nRotamers x nRuns cell structure where each cell contains
% nFrames x 1 dihedrals. MAKE SURE that the dihedrals contain no jumps due
% to periodicity, as the function will NOT account for those. For protein
% dihedrals, using [0,360] window rather than [-180,+180] is better for
% avoiding these periodicity jumps.
%
% * dihedrals_active: reference dihedral values for the active state,
% nRotamers x 1 cell structure. Can take more than 1 active state in the
% form of a vector in the cell element.
%
% * dihedrals_inactive: reference dihedral values for the inactive state,
% nRotamers x 1 cell structure.
%
% * tol:  data is considered in active/inactive state if it is within tol*std
% from the reference, where std is max([standard deviation of the data, 30]);
% default = 1.
%
% * sysState: nRotamers x nRuns cell structure where each cell contains 
% (3+k) x 1 vector: 1 for inactive state, 2 for intermediate/other state, 
% and 3 (and higher) for active states. k is the number of addiitonal
% active states.
%
% * index: nRotamers x nRuns cell structure where each cell contains
% nFrames x 1 vector indicating whether this frame is in the inactive (1),
% intermediate/other (2) state, or active state (3,4,...,k).
%%

% Set default value for tolerance
if nargin<4
    tol = 1;
end

numRuns = size(dihedrals,2);
nRotamers = size(dihedrals,1);
states = ['Active      ';'Inactive    ';'Intermediate'];
numStates = size(states,1); %Active, inactive, intermediate/other
sysState = cell(nRotamers,numRuns);
index = cell(nRotamers,numRuns); % Index to hold active/inactive data for every frame

% Rationale: if the dihedral is within active ref +- std, it is in active
% state, within inactive ref +- std, it is in inactive, else it is in
% intermediate/other

for i=1:nRotamers
    % How many states for this rotamer? 1 is for intermediates
    numStates = 1 + length(dihedrals_active{i}) +length(dihedrals_inactive{i});
    for j =1:numRuns
%         std_dihedral = std(dihedrals{i,j});
        std_dihedral = max([std(dihedrals{i,j}) 30]); % Degrees
        % 30 degress was decided after studying the average dispersion/width
        % of data in our simulations
        sysState{i,j} = zeros(numStates,1);
        
        % As many active states as there are:
        index_active = cell(length(dihedrals_active{i}),1);
        index_inactive = cell(length(dihedrals_active{i}),1);
        index_active_tot = zeros(length(dihedrals{i,j}),1); % To hold ALL active states
        index_inactive_tot = ones(length(dihedrals{i,j}),1); 
        
        if nargout>1 % Did the user ask for the detailed index?
            index{i,j} = zeros(length(dihedrals{i,j}),1);
        end
        
        for k=1:length(dihedrals_active{i})
        % Are active and inactive refs within +- stds of each other?
        if abs(dihedrals_inactive{i}-dihedrals_active{i}(k)) < 2*tol*std_dihedral
            
        % Options: 1- Take midpoint between them as boundary
        %          2- Take a moving average over a window and see if it's
        %          closer to active/inactive
        % 1- Midpoint method
        midpt = (dihedrals_inactive{i}+dihedrals_active{i}(k))/2;
        
          if dihedrals_active{i}(k)>dihedrals_inactive{i}
            index_inactive{k} = (dihedrals{i,j}>(dihedrals_inactive{i}-tol*std_dihedral) & ...
                dihedrals{i,j}<midpt);
            index_active{k} = (dihedrals{i,j}>midpt & ...
                dihedrals{i,j}<(dihedrals_active{i}(k)+tol*std_dihedral));
          else % Inactive > active
            index_inactive{k} = (dihedrals{i,j}>midpt & ...
                dihedrals{i,j}<(dihedrals_inactive{i}+tol*std_dihedral));
            index_active{k} = (dihedrals{i,j}>(dihedrals_active{i}(k)-tol*std_dihedral) & ...
                dihedrals{i,j}<midpt);
          end
          
        % 2- Moving average method:
        %
        % TO BE ADDED HERE
        %
        
        else % Active and inactive are far enough
            
            index_inactive{k} = (dihedrals{i,j}>(dihedrals_inactive{i}-tol*std_dihedral) & ...
                dihedrals{i,j}<(dihedrals_inactive{i}+tol*std_dihedral));
            index_active{k} = (dihedrals{i,j}>(dihedrals_active{i}(k)-tol*std_dihedral) & ...
                dihedrals{i,j}<(dihedrals_active{i}(k)+tol*std_dihedral));
        end
        
        if nargout>1 % Did the user ask for the detailed index?
            index{i,j} = index{i,j} + (2+k)*index_active{k}; % +3/4/5/... for active
        end
        
        % Define the active system state:
        sysState{i,j}(2+k) = sum(index_active{k});
        % Add to the total index:
        index_active_tot = (index_active_tot | index_active{k});
        % Keep the memory of the inactive state: (Useful when you have an
        % active state that overlaps with inactive and another that does
        % not)
        index_inactive_tot = (index_inactive_tot & index_inactive{k});
        end
        
        % Index other is not inactive and not ALL the active
        index_other = not(index_inactive_tot | index_active_tot);
        sysState{i,j}(1) = sum(index_inactive_tot);
        sysState{i,j}(2) = sum(index_other);
        
        if nargout>1 % Did the user ask for the detailed index?
            index{i,j} = index{i,j} + 1*index_inactive_tot; % +1 for inactive
            index{i,j} = index{i,j} + 2*index_other; % +2 for other states
        end
    end
end

end

