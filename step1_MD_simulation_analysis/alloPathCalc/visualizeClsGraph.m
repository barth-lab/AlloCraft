function  [pHere] = visualizeClsGraph(PDB,pathstruc,Gmatmajor,cls,options)
%visualizeClsGraph Visualizes allosteric channels on a 3D scatter of the
%protein structure. CAs are plotted as spheres colored by sequence from
%Nterm to Cterm according to the default colormap
% 
%% Usage:
% visualizeClsGraph(PDB,pathstruc,Gmatmajor,cls)
% [pHere] = visualizeClsGraph(PDB,pathstruc,Gmatmajor,cls,options)
%
%
% * PDB: pdb file obtained by: PDB = pdbread(pdbfile)
% 
% * pathstruc: structure containing pathways obtained by running: : graphanalysis.m
%
% * Gmatmajor: matrix with edgeweights: 1.1*MImax-MIres(i,j) obtained by
%  running: : graphanalysis.m
%
% * cls: channels to plot, can take either a scalar or a vector of channels
% to be visualized (matlab can get laggy if too many are plotted at the
% same time!)
%
% options:
% 
% 'MIFractionCutoff'

arguments
    PDB
    pathstruc
    Gmatmajor
    cls
    options.MIFractionCutoff = 0.85
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

%% Visualize the structure and the channels


clusters = [pathstruc(I(1:Npath)).cls];

figure

p3 = plot3(CAcoord(:,1),CAcoord(:,2),CAcoord(:,3),...
    'color',[0.8 0.5 0.],'Linewidth',2);
hold on
p3 = scatter3(CAcoord(:,1),CAcoord(:,2),CAcoord(:,3),...
    [],1:Nres,'o','filled','Linewidth',2);

% hcb = colorbar;
% hcb.Label.FontSize = 16;
% hcb.Label.String = "Sequence";

row = dataTipTextRow('Residue',1:Nres);
p3.DataTipTemplate.DataTipRows(end+1) = row;
 hold on
legendText{1} = 'Backbone';
legendText{2} = '';
count = 3;
for clsHere = cls
    paths = find(clusters==clsHere);
    members = unique([pathstruc(I(paths)).path]);
    
    
    G = graph(Gmatmajor(members,members));
    
    graphHere=G;
    x = CAcoord(members,1);
    y = CAcoord(members,2);
    z = CAcoord(members,3);
    LWidths = 2*graphHere.Edges.Weight/max(graphHere.Edges.Weight);
    
    pHere = plot(graphHere,'XData',x,'YData',y,'ZData',z,'LineWidth', LWidths,'MarkerSize',5);
%     pHere.NodeColor  = [0.8500 0.3250 0.0980];
%     pHere.EdgeColor  = [0.8500 0.3250 0.0980];
    pHere.NodeLabel = members; % Format by residue number
    legendText{count} = ['Channel ' num2str(clsHere)];
    count = count + 1;
end
xlabel('X [A]');
ylabel('Y [A]');
zlabel('Z [A]');
title(['Channels: ' num2str(cls) ', MIcut = ' num2str(options.MIFractionCutoff)],'FontSize',16)
legend(legendText);
formatplot2

axis('equal')
% axis off