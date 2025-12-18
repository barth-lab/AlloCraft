% Make a matrix out of channels
mydir = 'D:\Telework_library\dopamine_phase_3\md2pathMeta\buildGraphFromPipeline\';
sysName = 'double_mut';
nodeScore = zeros(Nres,1);
edgeWeight1a = zeros(Nres,Nres);
edgeWeight1b = zeros(Nres,Nres);

edgeWeight2 = MIres;
% Get the hubscores
for thisChnl = 1:length(channelstruc)
    nodeScore(channelstruc(thisChnl).hub) = nodeScore(channelstruc(thisChnl).hub)' + channelstruc(thisChnl).hubstrength;
end

% Way 1 of calculating edge connections: use interClsTopPathway to populate
% the connections between subclusters of channels
for i = 1:length(outstruc)
    for j = 1:(outstruc(i).Npath-1)
        edgeWeight1a(outstruc(i).path(j), outstruc(i).path(j+1)) = edgeWeight1a(outstruc(i).path(j), outstruc(i).path(j+1)) + 1;
    end
end

% Way 2: use clsConnectMat to weight connections between nodemajor of
% subclusters within channels, calculating the number of pathways
% connecting any given pair

for thisChnl = 1:length(channelstruc)
    connectMatHere = channelstruc(thisChnl).clsConnectMat;
    edgeWeight1b( [channelstruc(thisChnl).cls.nodemajor], [channelstruc(thisChnl).cls.nodemajor]) = ...
        edgeWeight1b( [channelstruc(thisChnl).cls.nodemajor], [channelstruc(thisChnl).cls.nodemajor]) + connectMatHere;
end

% Write out the variables:
writematrix(nodeScore, [mydir 'nodeScore_' sysName]);
writematrix(edgeWeight1a, [mydir 'edgeWeight1a_' sysName]);
writematrix(edgeWeight1b, [mydir 'edgeWeight1b_' sysName]);
writematrix(edgeWeight2, [mydir 'MIres_' sysName]);

%% For each cluster, write major nodes and pathways among them
count = 1;
for cls=1:length(channelstruc)
    clsConnectMat = channelstruc(cls).clsConnectMat;
    interClsTopPathway =channelstruc(cls).interClsTopPathway;
    ncls = channelstruc(cls).ncls;
    edgecount = 1;
    for i=1:ncls-1
        for j=i+1:ncls
            npath = clsConnectMat(i,j);
            if npath
                % get the major residues corresponding to both clusters
                majorResi = channelstruc(cls).cls(i).nodemajor;
                majorResj = channelstruc(cls).cls(j).nodemajor;
                
%                 % find the pathway that has resi & j as terminii
%                 ind = find( prod(double(ismember(endnode,...
%                     [majorResi,majorResj])),2) );
%                 outstruc(count) = pathstruc(I(ind));
%                 outstruc(count).cls = cls;
%                 [dum1,path,dum2] = graphshortestpath(G,majorResi,majorResj,...
%                     'Directed', 'false');
                pathind = interClsTopPathway(i,j);
                path = pathstruc(I(pathind)).path;
                outstruc(count).path = path;
                outstruc(count).Npath = length(path);
                %outstruc(count).MI = MIres(majorResi,majorResj);
                outstruc(count).MI = pathstruc(I(pathind)).MI;
                outstruc(count).meanpathMI = 0;
                outstruc(count).cls = cls;
                BondStrength(count,1) = cls;
                BondStrength(count,2) = edgecount;
                BondStrength(count,3) = edgecount+length(path)-1;
                BondStrength(count,4) = npath;
                edgecount = edgecount+length(path);
                count = count+1;
            end
        end
    end
end


%% Adding more than one edge weight to matlab graph object

% Make a random digraph using an adjacency matrix
rng default
A = rand(10) < 0.5;
D = digraph(A);
% Give the edges some weights
D.Edges.W1 = randi(10, numedges(D), 1);
D.Edges.W2 = randi(10, numedges(D), 1);
D.Edges.W3 = randi(10, numedges(D), 1);
% Because D.Edges has no Weight variable, shortestpath considers it to be unweighted
noweight = shortestpath(D, 1, 3)
% When we give the edges a Weight attribute, the shortest path changes
% The additional weight on the edges that made up the previous best path (noweight)
% make another route a shorter option
D.Edges.Weight = D.Edges.W3;
withweight = shortestpath(D, 1, 3)
