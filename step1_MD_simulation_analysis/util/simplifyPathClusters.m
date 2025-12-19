% after running AlloDy, you can run this script to visualize major hubs of
% channels on the structure
% Parameters
k = 4; %cluster ndx
topHubs = 10;
paperColors;

% Needed variables in workspace
options.mainEntry = mainEntry;
options.includeLigPathwayCalc = includeLigPathwayCalc;
receptorGpResIds;
receptorLigandResIds;
PDB; 
% or
CAcoord;

% Plot the backbone
figure
p3 = plot3(CAcoord(:,2),CAcoord(:,1),CAcoord(:,3),...
'color',paperColorPalette.Hex(4),'Linewidth',3);
if ~isempty(options.mainEntry)
    resMain = options.mainEntry.chains{1}.resIds; % chain 1 is always main chain
    if options.includeLigPathwayCalc % Make list from both main Chain and ligand chain
        resLig = options.mainEntry.chains{3}.resIds;  % chain 3 is always ligand
        resText = [options.mainEntry.chains{1}.formatResidues(resMain,'BWonly',true); ...
            options.mainEntry.chains{3}.formatResidues(resLig)];
    else
        resText = options.mainEntry.chains{1}.formatResidues(resMain,'BWonly',true);
    end
    row(1) = dataTipTextRow('Residue',resText);
    row(2) = dataTipTextRow('Residue',1:Nres);
    p3.DataTipTemplate.DataTipRows(end+1:end+2) = row;

else
    row = dataTipTextRow('Residue',1:Nres);
    p3.DataTipTemplate.DataTipRows(end+1) = row;

end


hold on


% set(gca,'zdir','reverse')


% Now do the channelstrucs:
hubHstrength = [ channelstruc(k).hub' channelstruc(k).hubstrength'];
sortedHubs = sortrows(hubHstrength,2,'descend');

% Are there any ligand binding resdiues here?
if find(ismember(sortedHubs(1:topHubs,1),receptorLigandResIds)) % we good
    hubsToPlot = sortedHubs(1:topHubs,:);
elseif find(ismember(sortedHubs(1:topHubs+10,1),receptorLigandResIds)) % Loosen up the constraints a bit
    ligBinding = find(ismember(sortedHubs(1:topHubs+10,1),receptorLigandResIds));
    hubsToPlot = [sortedHubs(1:topHubs,:); sortedHubs(ligBinding(1),:)];
else     
    hubsToPlot = sortedHubs(1:topHubs,:);
end


y = CAcoord(hubsToPlot(:,1),1);
x = CAcoord(hubsToPlot(:,1),2);
z = CAcoord(hubsToPlot(:,1),3);
s3 = scatter3(x,y,z,hubsToPlot(:,2),'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor',paperColorPalette.Hex(1), ...
    'MarkerFaceColor',paperColorPalette.Hex(2));

resText = [];
for i =1:size(hubsToPlot,1)
    temp = options.mainEntry.chains{1}.formatResidues(hubsToPlot(i,1),'BWonly',true);
    resText = [resText temp];
end
 row(1) = dataTipTextRow('Residue',resText);
 row(2) = dataTipTextRow('Strength',hubsToPlot(:,2));
 s3.DataTipTemplate.DataTipRows(end+1:end+2) = row;

xRange = range(CAcoord(:,2));
yRange = range(CAcoord(:,1));
zRange = range(CAcoord(:,3));
text(x + xRange/50, y + yRange/50, z + zRange/50, resText,'FontSize',16)
title(sprintf("Cluster %i, size = %.2f",k,channelstruc(k).size),'Fontsize',20)

axis xy
axis('equal')
 set(gca,'zdir','reverse')
% 
% channelstruc
% 
% receptorGpResIds
% receptorLigandResIds
% length(channelstruc(k).hub)
% 
% channelstruc
% 
% receptorGpResIds
% receptorLigandResIds



% Scatter the ligand and Gp binding residues

s3 = scatter3(CAcoord(receptorGpResIds,2),CAcoord(receptorGpResIds,1),CAcoord(receptorGpResIds,3),...
    [],'o','filled','Linewidth',2,'MarkerFaceColor',paperColorPalette.Hex(5));
if ~isempty(options.mainEntry)
    gpResText = options.mainEntry.chains{1}.formatResidues(receptorGpResIds,'BWonly',true);
    row = dataTipTextRow('Residue',gpResText);
else
    row = dataTipTextRow('Residue',receptorGpResIds);
end
s3.DataTipTemplate.DataTipRows(end+1) = row;

s3 = scatter3(CAcoord(receptorLigandResIds,2),CAcoord(receptorLigandResIds,1),CAcoord(receptorLigandResIds,3),...
    [],'o','filled','Linewidth',2,'MarkerFaceColor',paperColorPalette.Hex(9));
if ~isempty(options.mainEntry)
    lbResText = options.mainEntry.chains{1}.formatResidues(receptorLigandResIds,'BWonly',true);
    row = dataTipTextRow('Residue',lbResText);
else
    row = dataTipTextRow('Residue',receptorLigandResIds);
end
s3.DataTipTemplate.DataTipRows(end+1) = row;
