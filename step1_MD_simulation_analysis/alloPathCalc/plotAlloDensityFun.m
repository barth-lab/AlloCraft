function [density] = plotAlloDensityFun(PDB,pathstruc,setRes,options)
%plotAlloDensityFun Visualizes density of allosteric pathways in 3D
%containing setRes and (optionally) optRes. setRes could be binding site
%residues while optRes could be effector binding residues for example
% 
%% Usage:
% plotAlloDensityFun(PDB,pathstruc,setRes)
% [density] = plotAlloDensityFun(PDB,pathstruc,setRes,options)
%
% * density: 
%
% * PDB: pdb file obtained by: PDB = pdbread(pdbfile)
% 
% * pathstruc: structure containing pathways obtained by running: : graphanalysis.m
%
% options:
% 
% 'div'
% 'sigma'
% 'optRes'
% 'cut': pathway density cutoff for isosurface, could be a vector or a
% scalar
% 'MIFractionCutoff'
% 'title'
% 'mainEntry': To label residues by reference or generic numbering
% 'includeLigPathwayCalc'

arguments
    PDB
    pathstruc
    setRes
    options.div = 1
    options.sigma = 2
    options.optRes = []
    options.cut = 0.2
    options.MIFractionCutoff = 0.85
    options.title = []
    options.mainEntry = []
    options.includeLigPathwayCalc = false
end

%% Get Nres and CAcoords from PDB:

resno = [PDB.Model.Atom(:).resSeq];
atomname = char(PDB.Model.Atom(:).AtomName);
CAind = find( strncmp(cellstr(atomname),'CA',3) );
reslist = unique(resno);
Nres = length(reslist);

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

%% Sort pathways by MI and get Npaths corresponding to MIFractionCutoff


[~,I] = sortrows([pathstruc(:).MI]',-1);

cumMI = cumsum([pathstruc(I).MI]);
Npath = find(cumMI > cumMI(end) * options.MIFractionCutoff, 1);
%% obtain allosteric points density
div = options.div;
sigma = options.sigma;

%get allosteric points from pathways that contact the input set of residues

points = [];
if options.optRes
    for i=1:Npath
        pathpoints = pathstruc(I(i)).path;
        if sum(ismember(setRes,pathpoints)) && sum(ismember(options.optRes,pathpoints))
            points = [points pathpoints];
        end
    end
else
    for i=1:Npath
        pathpoints = pathstruc(I(i)).path;
        if sum(ismember(setRes,pathpoints))
            points = [points pathpoints];
        end
    end
end

%points = [pathstruc(I(1:Npath)).path];
mincoord = min(CAcoord);
maxcoord = max(CAcoord);
xcoord = mincoord(1)-5:div:maxcoord(1)+5;
ycoord = mincoord(2)-5:div:maxcoord(2)+5;
zcoord = mincoord(3)-5:div:maxcoord(3)+5;
density = zeros(length(xcoord),length(ycoord),length(zcoord));
[N1,N2,N3] = size(density);
for i=1:N1
    for j=1:N2
        for k=1:N3
            pi = [xcoord(i),ycoord(j),zcoord(k)];
            density(i,j,k) = sum( ...
                  pdf('norm',(CAcoord(points,1)-pi(1))/div,0,sigma)...
                .*pdf('norm',(CAcoord(points,2)-pi(2))/div,0,sigma)...
                .*pdf('norm',(CAcoord(points,3)-pi(3))/div,0,sigma) )...
                ;
        end
    end
end
% writematrix(density,'PathwayDensity.dat');

[X, Y, Z] = meshgrid(ycoord,xcoord,zcoord);

%% Make the isosurface from the density
cut = options.cut;

colormapHere = turbo;
figure
count = 1;
if length(cut)<9 % Just to make the isosurfaces look good
    alphaList = [0.1:0.1:(0.1*(length(cut)-1)) 1];
else
    alphaList = [linspace(0.1,0.8,length(cut) - 1) 1];
end

% Isovalues are the cutoffs of the maximum density
isovalues = cut* (max(density(:))-min(density(:)))+min(density(:)); 

for i=isovalues
    
    p = patch(isosurface(X,Y,Z,density,i));
     isonormals(X,Y,Z,density,p)
     if length(cut) == 1
        set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1);
     else
        colorHere = colormapHere(round(length(colormapHere)/length(cut)*(count-1))+1,:);
        set(p,'FaceColor',colorHere,'EdgeColor','none','FaceAlpha',alphaList(count));
     end
    
    daspect([1,1,1])
    view(-158,12); axis tight
    camlight
    lighting gouraud
    count = count + 1;
end
hold

% Plot the backbone
p3 = plot3(CAcoord(:,2),CAcoord(:,1),CAcoord(:,3),...
'color',[0.8 0.5 0.],'Linewidth',6);
if ~isempty(options.mainEntry)
    resMain = options.mainEntry.chains{1}.resIds; % chain 1 is always main chain
    if options.includeLigPathwayCalc % Make list from both main Chain and ligand chain
        resLig = options.mainEntry.chains{3}.resIds;  % chain 3 is always ligand
        resText = [options.mainEntry.chains{1}.formatResidues(resMain,'BWonly',true); ...
            options.mainEntry.chains{3}.formatResidues(resLig)];
    else
        resText = options.mainEntry.chains{1}.formatResidues(resMain,'BWonly',true);
    end
    row = dataTipTextRow('Residue',resText);
else
    row = dataTipTextRow('Residue',1:Nres);
end

p3.DataTipTemplate.DataTipRows(end+1) = row;

hold on

% Scatter the CAs colored by sequence
s3 = scatter3(CAcoord(:,2),CAcoord(:,1),CAcoord(:,3),...
    [],1:Nres,'o','filled','Linewidth',2);
if ~isempty(options.mainEntry)
    row = dataTipTextRow('Residue',resText);
else
    row = dataTipTextRow('Residue',1:Nres);
end
s3.DataTipTemplate.DataTipRows(end+1) = row;

colormap turbo

leg = legend(num2str(cut(:)),'Location','best','FontSize',16);
legend boxoff
title(leg,'Pathway density cutoff')

title(['Pathway density: ' options.title ', MIcut = ' num2str(options.MIFractionCutoff)],'FontSize',16)

axis xy
axis('equal')
axis off

end