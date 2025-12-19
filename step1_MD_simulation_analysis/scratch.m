% Clustering
% Possible change points:

% Filtering of MI (removeAverageMI?)

% Definition of pathOverlap

% Choosing number of clusters 
% The current method for specifying optimal cluster numbers is inherently
% flawed as the cluster separation efficiency keeps increasing to
% ridiculous numbers of clusters (~1000 in some cases)

% Currently I am using spectral cluster same as hierarchal by using the
% same overlap metric. Is it possible to make spectral act on the
% graph/pathways directly?

load('D:\Telework_library\dopamine_phase_3\1-d2_dop_WT\md2pathdev\alloPathCalc_Culled_data\workspace.mat');

% Question about construction of matrix
% Edge weights?

maxInd = 3;
L = laplacian(G);
[v,D] = eigs(L,maxInd);
Lambda = diag(D);

t = 1; % parameter for scaling the eigenvalues
g = @(x)(exp(-x)); % Function handle 
% g = @(x)(x^(-2)); % Function handle 

D2 = zeros(size(G.Nodes,1));
for i = 1:size(G.Nodes,1)-1
    for j = i+1:size(G.Nodes,1)
%         D2(i,j) = sum( g(t*Lambda') .* (v(i,:) - v(j,:)).^2);

        D2(i,j) = sum( arrayfun(@(x)g(x), t*Lambda') .* (v(i,:) - v(j,:)).^2);
    end
end
% Caclulate low dimensional embedding:
Psit = [];
for k = 1:maxInd
    Psit = [Psit g(t*Lambda(k))^0.5* v(:,k)];
end
figure; 
s = scatter(Psit(:,1),Psit(:,2),10,'filled');

dtRows = dataTipTextRow("Datapt idx",1:size(Psit,1));
s.DataTipTemplate.DataTipRows(end+1) = dtRows;

% D = sqrt(D2);
% figure
% imagesc(D)
%% Intro to spectral clustering:

rng('default'); % For reproducibility
n = 20;
y = [randn(n,2)*0.5+3;
    randn(n,2)*0.5;
    randn(n,2)*0.5-3];

%

dist_temp = pdist(y);
dist = squareform(dist_temp);

S = exp(-dist.^2);
issymmetric(S)

S_eps = S;
S_eps(S_eps<0.5) = 0;

%
k=3;
idx3 = spectralcluster(S_eps,k,'Distance','precomputed');
figure

gscatter(y(:,1),y(:,2),idx3);

%% Experiment with spectal clusters directly on graph of protein:

% Matrix of 1.1*maxMI - MI(i,j)
Gmatmajor;

% Try different formulations of edge weights?
% exp(-MI(i,j))
% 1/MI(i,j)
L = laplacian(G);




k = 10;
idx10 = spectralcluster(Gmatmajor,k,'Distance','precomputed');
% gscatter(X(:,1),X(:,2),idx3);

figure
pHere = plot(G,'XData',CAcoord(:,1),'YData',CAcoord(:,2),'ZData',CAcoord(:,3),'MarkerSize',5);
colormap hsv
xlabel('X [A]');
ylabel('Y [A]');
zlabel('Z [A]');
axis('equal')
pHere.NodeCData = idx10;

title('Spectral clustering on G','FontSize',20)
set(gca,'zdir','reverse')
set(gca,'ydir','reverse')

%% Betweenness centrality?
wbc = centrality(G,'betweenness','Cost',G.Edges.Weight);
n = numnodes(G);



pHere.NodeCData = 2*wbc./((n-2)*(n-1));
colormap(flip(autumn,1));
title('Betweenness Centrality Scores - Weighted','FontSize',20)
%% Try to cluster pathoverlap:

S = 1-pathoverlap;

G_S = graph(S(1:500,1:500));
figure;
plot(G_S)



k = 10;
[~,V_temp,D_temp] = spectralcluster(S,k,'Distance','precomputed','ClusterMethod','kmeans');

idx10 = spectralcluster(S,k,'Distance','precomputed','ClusterMethod','kmeans');

% How the hell do we visualize it?


%%
wbc = centrality(G_S,'betweenness','Cost',G_S.Edges.Weight);
n = numnodes(G_S);

figure
p = plot(G_S,'MarkerSize',5);
title('My GPCR')

p.NodeCData = 2*wbc./((n-2)*(n-1));
colormap(flip(autumn,1));
title('Betweenness Centrality Scores - Weighted')


%% Path overlap:

i = a(10);
j= b(10);

path1 = pathstruc(I(i)).path;
path2 = pathstruc(I(j)).path;
pathDistances = dismat(path1,path2);
[N1,N2] = size(pathDistances);
mindis = min( mean(min(pathDistances)),mean(min(pathDistances')) ); % POSSIBLE ISSUE: definition of pathdis 
pathdis(i,j) = mindis;
pathdis(j,i) = mindis;
near1 = sum(min(pathDistances')<nearcutoff);
near2 = sum(min(pathDistances )<nearcutoff);
% pathlendiff tries to include the width of the distribution of
% path lengths into the pathoverlap metric as a kind of
% normalization, what other metrics could I use?
pathoverlap(i,j) = max(1.0*near1/N1,1.0*near2/N2) - 0.4*abs(N1-N2)/pathlendiff; % Why? 
%         pathoverlap(i,j) = max(1.0*near1/N1,1.0*near2/N2) - 0.2*abs(N1-N2)/pathlendiff; % Why? 
pathoverlap(j,i) = pathoverlap(i,j);



LWidths =  G.Edges.Weight/max(G.Edges.Weight);
figure
% pHere = plot(G,'XData',CAcoord(:,1),'YData',CAcoord(:,2),'ZData',CAcoord(:,3),'MarkerSize',5,'EdgeCData',LWidths);
% colormap hot
% colorbar
pHere = plot(G,'XData',CAcoord(:,1),'YData',CAcoord(:,2),'ZData',CAcoord(:,3),'MarkerSize',5);
xlabel('X [A]');
ylabel('Y [A]');
zlabel('Z [A]');
axis('equal')
resText = mainEntry.chains{1}.formatResidues(resMain,'BWonly',true);
row = dataTipTextRow('Residue',resText);
pHere.DataTipTemplate.DataTipRows(end+1) = row;

pHere.NodeCData = 1:Nres;

% [T,p] = minspantree(G);
% highlight(pHere,T,'EdgeColor','r','LineWidth',1.5)


highlight(pHere,path1,'NodeColor','r','EdgeColor','r','LineWidth',1.5)
highlight(pHere,path2,'EdgeColor','g','LineWidth',1.5)

title(['Path overlap: ' num2str(pathoverlap(j,i))])

%% Effect of pathlendiff term:

figure('Position',[10 100 1600 400])
tiledlayout(1,3)
nexttile
imagesc(pathoverlap)
title('0.4')
nexttile
imagesc(pathoverlap02)
title('0.2')
nexttile
diff = pathoverlap02 - pathoverlap;
imagesc(diff)
title('Diff')
colorbar
colormap turbo


%% Try max flow rather than shortest path?
res1 = 162;
res2 = 261;
% For this application, edge weights need to increase with MI:
MIresDisCut = zeros(size(MIres));

for i=1:Nres-1
    for j=i+1:Nres
        if (dismat(i,j)<=disCutOff)
            MIresDisCut(i,j) = MIres(i,j);
%             MIresDisCut(i,j) = MIresHere{2}(i,j);
            MIresDisCut(j,i) =  MIresDisCut(i,j);
        end
    end
end
GMI = graph(MIresDisCut);

figure
pHere = plot(GMI,'XData',CAcoord(:,1),'YData',CAcoord(:,2),'ZData',CAcoord(:,3),'MarkerSize',5);
xlabel('X [A]');
ylabel('Y [A]');
zlabel('Z [A]');
axis('equal')

[mf, GF] = maxflow(GMI,res1,res2);
title(['Max flow = ' num2str(mf) ', s = ' num2str(res1) ', t = ' num2str(res2)])

resText = mainEntry.chains{1}.formatResidues(1:Nres,'BWonly',true);
row = dataTipTextRow('Residue',resText);
pHere.DataTipTemplate.DataTipRows(end+1) = row;
pHere.EdgeLabel = {};
highlight(pHere,GF,'EdgeColor',[0.4660 0.6740 0.1880],'LineWidth',2);
highlight(pHere,[res1 res2],'NodeColor',[0.8500 0.3250 0.0980]	,'LineWidth',1.5)

st = GF.Edges.EndNodes;
% labeledge(pHere,st(:,1),st(:,2),GF.Edges.Weight);

%% Max flow between BS and GPI:
mflow = zeros(length(BSres), length(GPIres));

 for i = 1:length(BSres)
     res1 = BSres(i);
     for j = 1:length(GPIres)
        res2 = GPIres(j);
        mflow(i,j) = maxflow(GMI,res1,res2);
     end
 end

 figure
 h = heatmap(BSres,GPIres,mflow');
 h.FontSize=14;
 xlabel('BSres')
 ylabel('GPIres')
 title('\fontsize{25} Max flow map')
 colormapCustom('seq','inputColor1', [0 0.65 0.65])
 savefig(fullfile(md2pathdir,'max_flow_bs_gpi'));
% 
%  figure
%  h = heatmap(BSres,GPIres,MIres(BSres,GPIres)');
%  h.FontSize=14;
%  xlabel('BSres')
%  ylabel('GPIres')
%  title('\fontsize{25} Max flow map')
%  colormapCustom('seq','inputColor1', [0 0.65 0.65])

%% Max flow full map between "allosteric" residues: does this even make sense?
mflowFull = zeros(Nres);
for i=1:Nres
    for j=i+1:Nres
%         if MIres(i,j)>0 && (dismat(i,j)>disCutOff)
            mflowFull(i,j) = maxflow(GMI,i,j);
%         end
            mflowFull(j,i) = mflowFull(i,j);
    end
end

 figure
%  h = heatmap(1:Nres,1:Nres,mflowFull);
imagesc(mflowFull); colormap turbo; colorbar
%  h.FontSize=14;
 xlabel('Residue')
 ylabel('Residue')
 title('\fontsize{25} Max flow map')
 drawTMhelices(1:Nres,settings.helices,1:Nres)
 drawTMhelices(1:Nres,settings.helices,1:Nres,'Y')
  savefig(fullfile(md2pathdir,'max_flow_map'));

%  colormapCustom('seq','inputColor1', [0 0.65 0.65])

%% Diffusion maps

load('D:\Telework_library\dopamine_phase_3\1-d2_dop_WT\md2pathdev\alloPathCalc_Culled_data\workspace.mat');

% Question about construction of matrix
% Edge weights?

maxInd = 3;
L = laplacian(G);
[v,D] = eigs(L,maxInd);
Lambda = diag(D);

t = 1; % parameter for scaling the eigenvalues
g = @(x)(exp(-x)); % Function handle 
% g = @(x)(x^(-2)); % Function handle 

D2 = zeros(size(G.Nodes,1));
for i = 1:size(G.Nodes,1)-1
    for j = i+1:size(G.Nodes,1)
%         D2(i,j) = sum( g(t*Lambda') .* (v(i,:) - v(j,:)).^2);

        D2(i,j) = sum( arrayfun(@(x)g(x), t*Lambda') .* (v(i,:) - v(j,:)).^2);
    end
end
% Caclulate low dimensional embedding:
Psit = [];
for k = 1:maxInd
    Psit = [Psit g(t*Lambda(k))^0.5* v(:,k)];
end
figure; 
s = scatter(Psit(:,1),Psit(:,2),10,'filled');

dtRows = dataTipTextRow("Datapt idx",1:size(Psit,1));
s.DataTipTemplate.DataTipRows(end+1) = dtRows;

% D = sqrt(D2);
% figure
% imagesc(D)
%% Read sequence from database and pdb

[~, seqRef] = fastaread(fullfile(settings.databasePath, settings.systemName + ".fasta"));
% Or fetch sequence from online database:
% This one fetches the DD2R sequence for example
  D = getgenpept('NP_000786')

 seqTest = database.entries{1}.seq{Chains.receptor};

% align and compare
[a,b] = nwalign(seqRef,seqTest);

tempseq1 = b(1,1:end);
tempseq2 = b(3,1:end);


% Any mutation/difference in sequence 
% matches = (tempseq1==tempseq2);
% mutations = [tempseq1(matches==0); tempseq2(matches==0)];
% mutPos =  find(matches==0);
% tempCell = cell(length(mutPos),1);
% for i=1:length(mutPos)
%     tempCell{i} = [mutations(2,i) num2str(mutPos(i)) mutations(1,i)];
% end

% Only substitutions (what we usually care about):
substitutions = tempseq1~=tempseq2 & tempseq1~='-' & tempseq2~='-';
mutPos = find(substitutions);
tempCell = cell(length(mutPos),1);
for i=1:length(mutPos)
    
    tempCell{i} = [tempseq1(mutPos(i)) num2str(mutPos(i)) tempseq2(mutPos(i))];
end


%% Attempt at trajectory clustering with hierarchical clustering:

% modift simulation.computeRmsdAllFramePairs
% Choose atoms 
receptorAtoms = mainChain.getAtoms();
RMSD = mainSim.computeRmsd(receptorAtoms);

% 
% atom_indices = [a.index for a in traj.topology.atoms if a.element.symbol != 'H']
% distances = np.empty((traj.n_frames, traj.n_frames))
% for i in range(traj.n_frames):
%     distances[i] = md.rmsd(traj, traj, i, atom_indices=atom_indices)

%% Contact GPCRdb programatically:

% Examples
% Get protein data:
url = 'https://gpcrdb.org/services/protein/adrb2_human/';
myProt = webread(url);

% Get protein alignment by specifying protein and segments
url = 'https://gpcrdb.org/services/alignment/protein/adrb1_human,adrb2_human,adrb3_human/TM3,TM5,TM6/';
response = webread(url);

% Get list of all species in the DB:
url = 'https://gpcrdb.org/services/species/';
speciesList = webread(url);

% Get list of available structures in GPCRdb:
url = 'https://gpcrdb.org/services/structure/';
structList = webread(url);

% Get list of representative structures for each protein and activation
% state:
url = 'https://gpcrdb.org/services/structure/representative/';
repStructList = webread(url);

% Get a single structure instance
url = 'https://gpcrdb.org/services/structure/6VMS/';
strucHere = webread(url);

% get  residue list: (This is NOT BW, this is GPCRdb(A) )
url = 'https://gpcrdb.org/services/residues/adrb2_human/';
resList = webread(url);

% get extended residue list: (contains all residue numbering schemes)
url = 'https://gpcrdb.org/services/residues/extended/adrb2_human/';
resList = webread(url);

% Get snake plot: snake, SNAKE, SNAAAAAAAAAAAAAAAAAAAAKE
url = 'https://gpcrdb.org/services/plot/snake/adrb1_human/';
snek = webread(url);


% Get list of protein families:
url = 'https://gpcrdb.org/services/proteinfamily/';
familyList = webread(url);
% Get Protein from certain family: Aminergic recceptors:
url = ['https://gpcrdb.org//services/alignment/family/' familyList(3).slug '/'];
options = weboptions('Timeout',30);
familyAlign = webread(url,options);


% Upload PDB file and get generic numbering: (Internal error now :( )
% Read file contents
dataFile = 'D:\Telework_library\dopamine_phase_3\1-d2_dop_WT\prot.pdb';
try
    fid = fopen(dataFile, 'r');
    data = char(fread(fid)');
    fclose(fid);
catch someException
    throw(addCause(MException('unableToReadFile','Unable to read input file.'),someException));
end
% Upload with webwrite
options = weboptions('MediaType', 'application/json');
tempOutput = webwrite('https://gpcrdb.org/services/structure/assign_generic_numbers', data,options);