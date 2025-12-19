resno = [PDB.Model.Atom(:).resSeq];
atomname = char(PDB.Model.Atom(:).AtomName);
CAind = find( strncmp(cellstr(atomname),'CA',3) );
reslist = unique(resno);
Nres = length(reslist);

% compute residue distance matrix
dismat = zeros(Nres,Nres);
for i=1:Nres-1
    xi = PDB.Model.Atom(CAind(i)).X;
    yi = PDB.Model.Atom(CAind(i)).Y;
    zi = PDB.Model.Atom(CAind(i)).Z;
    for j=i:Nres
        xj = PDB.Model.Atom(CAind(j)).X;
        yj = PDB.Model.Atom(CAind(j)).Y;
        zj = PDB.Model.Atom(CAind(j)).Z;
        dis = norm([xi yi zi] - [xj yj zj]);
        dismat(i,j) = dis;
        dismat(j,i) = dis;
    end
end
CAcoord = [PDB.Model.Atom(CAind).X; PDB.Model.Atom(CAind).Y;...
    PDB.Model.Atom(CAind).Z]';

% analyze cluster sizes
clusters = [pathstruc(I(1:Npath)).cls];
ncls = max(clusters);
clsSize = zeros(ncls,1);
for i=1:ncls
    paths = find(clusters==i);
    clsSize(i) = length(paths);
    %members = unique([pathstruc(I(paths)).path]);
end
[temp,clsSizeOrder] = sortrows(clsSize,-1);

% find MIs of each pathway
count = 1;
for i=1:ncls
    paths = find(clusters==i);
    clear temp
    for j=1:length(paths)
        resi = pathstruc(I(paths(j))).path(1);
        resj = pathstruc(I(paths(j))).path(end);
        clsMI(count,:) = [i,clsSize(i),MIres(resi,resj)];
        temp(j) = MIres(resi,resj);
        count = count+1;
    end
    [temp1,index] = max(temp);
    pathMaxMI(i) = (paths(index));
end

clsMImean = zeros(ncls,1);
clsMIstd = zeros(ncls,1);
clsMIsorted = sortrows(clsMI,[1,-2,-3]);
for i=1:ncls
    clsMImean(i) = mean(clsMI(clsMI(:,1)==i,3));
    clsMIstd(i) = std(clsMI(clsMI(:,1)==i,3));
end
if diagnosticsOn
    figure; 
    scatter(1:ncls,clsMImean,clsSize,'filled','MarkerFaceAlpha',0.5)
    hold on
    errorbar(1:ncls,clsMImean,clsMIstd./max(sqrt(clsSize-1),1),'.')
    xlabel('Cluster')
    ylabel('Mean MI')
    title('Mean cluster MI')
    set(gca,'FontSize',16)
    savefig(fullfile(pathCalcdir,"cluster_mean_MI"));
    print2pdf(fullfile(pathCalcdir,"cluster_mean_MI"));
end


% compute spread of each cluster
for cl=1:ncls
    sumdis = 0;
    count = 0;
    paths = find(clusters==cl);
    for i=1:clsSize(cl)-1
        path1 = pathstruc(I(paths(i))).path;
        for j=i+1:clsSize(cl)
            path2 = pathstruc(I(paths(j))).path;
            pathDistances = dismat(path1,path2);
            [N1,N2] = size(pathDistances);
            if(N1<N2)
                mindis = min(pathDistances');
            else
                mindis = min(pathDistances);
            end
            sumdis = sumdis+sum(mindis);
            count = count+length(mindis);
        end
    end
    clsspread(cl) = sumdis/count;
    if isnan(clsspread(cl))
        clsspread(cl) = 0;
    end
end


% write out cluster stat
% writematrix([clsSizeOrder';clsSize(clsSizeOrder)';clsMImean(clsSizeOrder); ...
%     clsspread(clsSizeOrder)]', fullfile(pathCalcdir,"cls_stat.txt"),'Delimiter'," ")
writematrix([clsSizeOrder';clsSize(clsSizeOrder)';clsMImean(clsSizeOrder)'; ...
    clsspread(clsSizeOrder)]', fullfile(pathCalcdir,"cls_stat.txt"),'Delimiter'," ")

% Find which residues are heavily involved in cross-channel communication
% By counting times in which residue shows up in every cluster
res_cross_freq = zeros(Nres,ncls+2);
for i=1:Nres
    freq = zeros(ncls,1);
    for j=1:ncls
        ind = ([pathstruc(I(1:Npath)).cls]==j);
        freq(j) = sum([pathstruc(I(ind)).path]==i);
    end
    res_cross_freq(i,1) = sum(freq);
    res_cross_freq(i,2) = std(freq);
    res_cross_freq(i,3:end) = freq(:);
end
% sort cross residues by standard deviation
[temp,res_cross_freq_order] = sortrows(res_cross_freq,-2);

writematrix( [res_cross_freq_order(1:20),...
    res_cross_freq(res_cross_freq_order(1:20),:)], ...
    fullfile(pathCalcdir,"res_cross_freq.txt"),'Delimiter'," ");

% cluster nodes in each cluster, size>=10
clear channelstruc
for cls = find(clsSize>=10)'
    
    % collect nodes
    paths = find(clusters==cls);
    for i=1:length(paths)
        nodes(2*i-1) = pathstruc(I(paths(i))).path(1);
        nodes(2*i) = pathstruc(I(paths(i))).path(end);
    end
    nodes = unique(nodes);
    
    % cluster nodes
    X = CAcoord(nodes,:);
    % determine optimum clustering; cycle through no. of clusters 2-50
    maxcls = min(10,length(X));
    %nodedis = zeros(maxcls,1);
    avgsil = zeros(maxcls,1);
    for n=2:maxcls
        T = kmeans(X,n,'Replicates',15);
        sil = silhouette(X,T);
        avgsil(n) = mean(sil);
    end
    
    % find local maxima in average silhouette for ncluster < npoints/2
    optsil = 0;
    optncluster = 2;
    for i=3:length(avgsil)-1
        [~,maxindex] = max( avgsil(i-1:i+1) );
        if maxindex == 2
            if avgsil(i) > optsil
                optsil = avgsil(i);
             optncluster = i;
            end
        end
    end
    ncluster = optncluster;
    
    % do the actual clustering
    [T,centroid] = kmeans(X,ncluster,'Replicates',15);
    
    % determine connectivity among clusters
    clsConnectMat = zeros(ncluster,ncluster);
    
    % top pathways in between different cluster 
    interClsTopPathway = zeros(ncluster,ncluster);
    
    pathnodes = zeros(length(paths),2);
    for i=1:length(paths)
        pathnodes(i,:) = [pathstruc(I(paths(i))).path(1) ...
            pathstruc(I(paths(i))).path(end)];
    end
    for i=1:ncluster-1
        nodesi = nodes(T==i);
        connpartner = 0;
        maxconn = 0;
        for j=i+1:ncluster
            if i~=j
                nodesj = nodes(T==j);
                memberpaths = ( ...
                      (ismember(pathnodes(:,1),nodesi)...
                    | ismember(pathnodes(:,2),nodesi))...
                    & (ismember(pathnodes(:,1),nodesj)...
                    | ismember(pathnodes(:,2),nodesj)) );
                clsConnectMat(i,j) = sum(memberpaths);
                clsConnectMat(j,i) = sum(memberpaths);
                
                if sum(memberpaths)
                    % find the topmost pathway in terms of MI
                    pathindex = paths(memberpaths);
                    [temp,ind] = max([pathstruc(I(pathindex)).MI]);
                    interClsTopPathway(i,j) = pathindex(ind);
                    interClsTopPathway(j,i) = pathindex(ind);
                else
                    interClsTopPathway(i,j) = 0;
                    interClsTopPathway(j,i) = 0;
                end
            end
        end
    end
    
    % store cluster data in structure
    % find strength of each node in each of the clusters
    for i=1:ncluster
        nodesi = nodes(T==i);
        channelstruc(cls).cls(i).node = nodesi;
        
        % for each residue, find how involved they are, i.e. how many
        % connections they are involved in
        strength = zeros(1,length(nodesi));
        nodeMI = zeros(1,length(nodesi));
        for j=1:length(nodesi)
            strength(j) = ...
                sum((pathnodes(:,1)==nodesi(j)) | (pathnodes(:,2)==nodesi(j)));
            nodeMI(j) = sum(MIres(nodesi(j),nodes));
        end
        channelstruc(cls).cls(i).strength = strength;
        channelstruc(cls).cls(i).MI = nodeMI;
        [temp,ind] = max(strength);
        channelstruc(cls).cls(i).nodemajor = nodesi(ind);
    end
    channelstruc(cls).ncls = ncluster;
    channelstruc(cls).clsConnectMat = clsConnectMat;
    channelstruc(cls).interClsTopPathway = interClsTopPathway;
    channelstruc(cls).size = clsSize(cls);
    channelstruc(cls).spread = clsspread(cls);
    channelstruc(cls).MI = clsMImean(cls);
    
end

% Identify hub residues
respathwaycount = zeros(Nres,1);
for i=1:length(channelstruc)
    ind = find([pathstruc(I(1:Npath)).cls] == i);
    pathresidues = unique([pathstruc(I(ind)).path]);
    channelstruc(i).hub = pathresidues;
    for j=1:length(pathresidues)
       count = 0;
       for k=ind          
            if MIWeightPaths
                count = count + ismember(pathresidues(j),pathstruc(I(k)).path)*pathstruc(I(k)).MI;
            else
                count = count + ismember(pathresidues(j),pathstruc(I(k)).path);  
            end
       end
       channelstruc(i).hubstrength(j) = count;
       respathwaycount(pathresidues(j)) = respathwaycount(pathresidues(j))+count;
    end
end

cls = 1;
clsConnectMat = channelstruc(cls).clsConnectMat;
ncluster = channelstruc(cls).ncls;

% write channel info to file
outstring = [];
for i=1:length(channelstruc)
    for j=1:length(channelstruc(i).cls)
        for k=1:length(channelstruc(i).cls(j).node)
            outstring = [outstring; i channelstruc(i).size channelstruc(i).spread channelstruc(i).MI ...
                channelstruc(i).cls(j).node(k) ...
                channelstruc(i).cls(j).strength(k) ...
                channelstruc(i).cls(j).MI(k) ];
        end
    end
end
writematrix( outstring, fullfile(pathCalcdir,"channelinfo.dat"),'Delimiter',"\t");

% write hub info
outstring = [];
for i=1:length(channelstruc)
    hub = channelstruc(i).hub;
    hubstrength = channelstruc(i).hubstrength;
    for j=1:length(hub)
        outstring = [outstring; i channelstruc(i).size channelstruc(i).spread channelstruc(i).MI ...
            hub(j) hubstrength(j)];
    end
end
writematrix( outstring, fullfile(pathCalcdir,"hubinfo.dat"),'Delimiter',"\t");
