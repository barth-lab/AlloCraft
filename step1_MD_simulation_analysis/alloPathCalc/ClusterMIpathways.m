% Clustering algorithms:
% Spectral clustering performs fine with pathoverlap
% But it does not work at all with similarity from step 2 of clustering,
% and it simply gives 1 huge cluster 


% compute inter-pathway distances and overlaps
diffFactor = 0.2; % Scales penalty to path length differences during overlap calculations
% diffFactor > 1 - overlapcutoff --> smaller paths cannot integrate with
% larger overlapping paths
% diffFactor < 1 - overlapcutoff --> smaller paths can integrate with
% larger overlapping paths

% Define overlapcutoff for jaccard:
% overlapcutoff =(min([pathstruc(I(1:Npath)).Npath])/max([pathstruc(I(1:Npath)).Npath]))+0.1;

pathdis = zeros(Npath,Npath);
pathoverlap = ones(Npath,Npath);
pathlendiff = max([pathstruc(I(1:Npath)).Npath]) - min([pathstruc(I(1:Npath)).Npath]);
for i=1:Npath-1
    path1 = pathstruc(I(i)).path;
    if mod(i,100) == 0 % Debug
        i
    end
    for j=i+1:Npath
        path2 = pathstruc(I(j)).path;
        pathDistances = dismat(path1,path2);
        [N1,N2] = size(pathDistances);
        mindis = min( mean(min(pathDistances)),mean(min(pathDistances')) ); 
        pathdis(i,j) = mindis;
        pathdis(j,i) = mindis;
        near1 = sum(min(pathDistances')<nearcutoff);
        near2 = sum(min(pathDistances )<nearcutoff);


        % pathlendiff tries to include the width of the distribution of
        % path lengths into the pathoverlap metric as a kind of
        % normalization, what other metrics could I use?

        pathoverlap(i,j) = max(1.0*near1/N1,1.0*near2/N2) - diffFactor*abs(N1-N2)/pathlendiff; % Why? 
        pathoverlap(j,i) = pathoverlap(i,j);

        % Try "raw" Jaccard similarity:
%         pathoverlap(i,j) = sum(ismember(path1,path2))/length(unique([path1 path2])); % The equations work for paths containing non-repeating elements
%         pathoverlap(j,i) = pathoverlap(i,j);

        % Try Jaccard with distance based overlap as intersection
%         pathoverlap(i,j) = max(near1,near2)/length(unique([path1 path2]));
%         pathoverlap(j,i) = pathoverlap(i,j);
    end
end

% Clustering diagnostics:

% highlight(pHere,path1,'NodeColor','r','EdgeColor','r')
% highlight(pHere,path2,'NodeColor','g','EdgeColor','g')
% 
% 
% figure; dendrogram(linkagemat,1000)

pathoverlap(pathoverlap<overlapcutoff) = 0; % default overlapcutoff would be 0.75
% generate pdist format similarity matrix (required when using similarity!)
save(fullfile(pathCalcdir,"workspace.mat"),'pathoverlap','-append');

overlapvect = squareform(1-pathoverlap,'tovector');
%overlapvect = squareform(pathdis,'tovector');
linkagemat = linkage(overlapvect,'average'); 

cophCorr = cophenet(linkagemat,overlapvect); % How well do the distances in 
% linkage reflect the distances in the data? A low number means we need to
% change the way we define overlap or linkage method


% determine inconsistent coeffs: clusters with high inconsistency naturally
% divide into separate clusters, since the height difference along the tree
% is higher. Using a specific value of inconsistency cutoff is a BAD idea,
% since we will end of either with 1 cluster (if cutoff = 1.16) or with
% >1000 clusters (if cutoff = 1.15) in the data I tested. So it's HIGHLY
% sensitive to the choice of cutoff

% inconscoeff = inconsistent(linkagemat);
% inconscoeff = inconscoeff(:,4);
% dendrogram(linkagemat)
% T = cluster(linkagemat,'cutoff',1.1);



% determine optimal 'cutoff' with another method: 

%% cutofflist = 2:min(200,Npath/10);
cutofflist = [2:1:50 52:2:100 105:5:200 210:10:round(Npath/10)];

meansepratio = zeros(length(cutofflist),1);
meanintraoverlap = zeros(length(cutofflist),1);
meaninteroverlap = zeros(length(cutofflist),1);
maxinteroverlap = zeros(length(cutofflist),1);
minintraoverlap = zeros(length(cutofflist),1);
clssize = zeros(length(cutofflist),1);
k=1;

figure
hold
set(gca,'Fontsize',16)
xlabel('No. of clusters','Fontsize',16)
ylabel('Pathway overlap','Fontsize',16)
title('Clustering step 1')

for cut=cutofflist

    % Hierarchical:
    %T = cluster(linkagemat,'cutoff',cut);
    T = cluster(linkagemat,'maxclust',cut);
    
    % Spectral:
%     T = spectralcluster(1-pathoverlap,cut,'Distance','precomputed','ClusterMethod','kmeans');


    ncls = max(T);
    % compute intra-cluster overlaps
    intraoverlap = zeros(ncls,1);
    for i=1:ncls
        ind = (T==i);
        overlapmat = tril(pathoverlap(ind,ind),-1)...
            +triu(2*ones(sum(ind),sum(ind)));
        intraoverlap(i) = mean(overlapmat(overlapmat<2));
        if isnan(intraoverlap(i))
            intraoverlap(i) = 1;
        end
    end
    minoverlap = min(intraoverlap);
    
    % compare inter & intra cluster overlaps
    %interoverlap = ones(ncls,ncls);
    avginter = 0; avgsep = 0; count = 0;
    maxoverlap = 0;
    for i=1:ncls-1
        indi = (T==i);
        for j=i+1:ncls
            indj = (T==j);
            overlapmat = pathoverlap(indi,indj);
            interoverlap = mean(mean(overlapmat));
            if isnan(interoverlap)
                interoverlap = 1;
            end
%             sepratio = (intraoverlap(i)+intraoverlap(j))...
%                 /2/interoverlap;
            sepratio = (intraoverlap(i)+intraoverlap(j))/2-1.*interoverlap;
            if isinf(sepratio)
                sepratio = 0;
            end
            %interoverlap(j,i) = interoverlap(i,j);
            if interoverlap
                avginter = avginter+interoverlap;
                avgsep = avgsep + sepratio;
                count = count+1;
            end
            if interoverlap>maxoverlap
                maxoverlap = interoverlap;
            end
%             if interoverlap<minoverlap
%                 minoverlap = interoverlap;
%             end
        end
    end
    meansepratio(k) = avgsep/count;
    meaninteroverlap(k) = avginter/count;
    meanintraoverlap(k) = mean(intraoverlap);
    maxinteroverlap(k) = maxoverlap;
    minintraoverlap(k) = minoverlap;
    %meansepratio(k) = minoverlap;
    clssize(k) = ncls;
    if ncls==1
        meansepratio(k) = 1;
        meaninteroverlap(k) = 0;
    end
    plot(cutofflist(1:k),meanintraoverlap(1:k),'b-')
    plot(cutofflist(1:k),meaninteroverlap(1:k),'r-')
    plot(cutofflist(1:k),minintraoverlap(1:k),'b--')
    plot(cutofflist(1:k),maxinteroverlap(1:k),'r--')
    plot(cutofflist(1:k),meansepratio(1:k),'k-')
    drawnow
    k = k+1;
end

axis([0 cutofflist(end) 0 1.1])
legend('Mean intra-cluster overlap','Mean inter-cluster overlap',...
    'Min intra-cluster overlap','Max inter-cluster overlap',...
    'Cluster separation efficiency')
box on
%%
% optimum clustering should be intraoverlap>=75%, interoverlap<=50%,
% max(sep ratio)
[temp,ind] = max(meansepratio(meanintraoverlap>0.75...
    & meaninteroverlap<0.5));
temp = cutofflist(meanintraoverlap>=0.75 & meaninteroverlap<=0.5);
nclsoptim = temp(ind);

% Hierarchical
T = cluster(linkagemat,'maxclust',nclsoptim);

% Spectral:
% T = spectralcluster(1-pathoverlap,nclsoptim,'Distance','precomputed','ClusterMethod','kmeans');


% assign clusters to pathstruc
for i=1:Npath
    pathstruc(I(i)).cls = 0;
end
for i=1:Npath
    pathstruc(I(i)).cls = T(i);
end
ncls = nclsoptim;

% identify hubs
for i=1:ncls
    ind = find([pathstruc(I(1:Npath)).cls] == i);
    pathresidues = unique([pathstruc(I(ind)).path]);
    channelstruc(i).hub = pathresidues; % This channelstruc is temporary 
    for j=1:length(pathresidues)
       count = 0;
       for k=ind % this hubstrength is used for channel similarity, so need
           % to weigh it by MI
            count = count + ismember(pathresidues(j),pathstruc(I(k)).path);
       end
       channelstruc(i).hubstrength(j) = count; 
       %respathwaycount(pathresidues(j)) = respathwaycount(pathresidues(j))+count;
    end
end
% merge overlapping clusters
similarity = ones(ncls,ncls);
clear hub hubstrength
% identify hubs
for i=1:ncls
    ind = find([pathstruc(I(1:Npath)).cls] == i);
    pop = length(ind);
    pathresidues = unique([pathstruc(I(ind)).path]);
    hub{i} = pathresidues;
    for j=1:length(pathresidues)
       count = 0;
       for k=ind
            count = count + ismember(pathresidues(j),pathstruc(I(k)).path);
       end
       hubstrength{i}(j) = count/pop*1.0;
    end
end

% calculate channel similarity
for i=1:ncls-1
    resi = hub{i};
    strengthi = zeros(1,Nres);
    strengthi(resi) = hubstrength{i};
    for j=i+1:ncls
        resj = hub{j};
        strengthj = zeros(1,Nres);
        strengthj(resj) = hubstrength{j};
        common = resj(ismember(resj,resi));
        all = unique([resi resj]);
        weight_common = (strengthi(common)+strengthj(common))/2;
        weight_all = (strengthi(all)+strengthj(all))/2;
        similarity(i,j) = 1.0*sum(weight_common)/sum(weight_all);
        %similarity(i,j) = 1.0*length(common)/length(all);
        similarity(j,i) = similarity(i,j);
    end
end
similarity(similarity<0.5) = 0;

% determine optimum clustering given channel similarity 
simvect = squareform(1-similarity,'tovector');
linkagemat = linkage(simvect,'average');
inconscoeff = inconsistent(linkagemat);
inconscoeff = inconscoeff(:,4);
cutofflist = 2:1:ncls;
meansepratio = zeros(length(cutofflist),1);
for cut=cutofflist
    
    
    % Hierarchical
    T = cluster(linkagemat,'maxclust',cut);
    
    % Spectral:
%     T = spectralcluster(1-similarity,cut,'Distance','precomputed','ClusterMethod','kmeans');

    ncls1 = max(T);
    intraoverlap = zeros(ncls1,1);
    for i=1:ncls1
        ind = (T==i);
        temp = similarity(ind,ind);
        intraoverlap(i) = mean(temp(temp<1));
        if(isnan(intraoverlap(i)))
            intraoverlap(i) = 0;
        end
    end
    avginter = 0; avgsep = 0; count = 0;
    for i=1:ncls1-1
        indi = (T==i);
        for j=i+1:ncls1
            indj = (T==j);
            temp = similarity(indi,indj);
            interoverlap = mean(temp(temp<1));
            if isnan(interoverlap)
                interoverlap = 1;
            end
            sepratio = (intraoverlap(i)+intraoverlap(j))/2-1.*interoverlap;
            if isinf(sepratio)
                sepratio = 0;
            end
            if interoverlap
                avginter = avginter+interoverlap;
                avgsep = avgsep + sepratio;
                count = count+1;
            end
        end
    end
    meansepratio(cut) = avgsep/count;
end
figure
plot(meansepratio,'-')
set(gca,'Fontsize',16)
xlabel('No. of clusters','Fontsize',16)
ylabel('Cluster separation efficiency','Fontsize',16)
title('Clustering step 2')

[maxsepratio,ncls1] = max(meansepratio);
hold on
scatter(ncls1,maxsepratio,50,'fill','r')
legend('meansepratio','Optimum cluster nbr','location','best')
legend boxoff

% Hierarchical
T = cluster(linkagemat,'maxclust',ncls1);

% Spectral:
% T = spectralcluster(1-similarity,ncls1,'Distance','precomputed','ClusterMethod','kmeans');

% OR try DBscan
% epsilon = 0.005; % no easy way to tell epsilon
% minpts = 20; 
% 
% [T, corepts] = dbscan(1-pathoverlap,epsilon,minpts,'Distance','precomputed');
% [max(T) sum(T==-1)]
% 
% for i=1:Npath
%     pathstruc(I(i)).cls = T(i);
% end
% ncls = max(T);


% assign new cls numbers
for i=1:Npath
    temp = pathstruc(I(i)).cls;
    pathstruc(I(i)).cls = T(temp);
end
ncls = ncls1;



% sort clusters by size
clslist = [pathstruc(I(1:Npath)).cls];
clssize = zeros(ncls,1);
for i=1:ncls
    clssize(i) = sum(clslist==i);
end
[temp,ind] = sortrows(clssize,-1);

% reassign clusters
for i=1:Npath
    clslist(i) = find(ind==clslist(i));
end
for i=1:Npath
    pathstruc(I(i)).cls = clslist(i);
end
