clear outstruc outstring BondStrength endnode
% get terminal ends of all pathways
endnode = zeros(length(I),2);
for i=1:length(I)
    if pathstruc(I(i)).MI
        endnode(i,1) = pathstruc(I(i)).path(1);
        endnode(i,2) = pathstruc(I(i)).path(end);
    end
end

% For each cluster, write major nodes and pathways among them
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

% Identify hub residues
respathwaycount = zeros(Nres,1);
for i=1:length(channelstruc)
    ind = find([pathstruc(I(1:Npath)).cls] == i);
    pathresidues = unique([pathstruc(I(ind)).path]);
    channelstruc(i).hub = pathresidues;
    for j=1:length(pathresidues)
       count = 0;
       for k=ind
            count = count + ismember(pathresidues(j),pathstruc(I(k)).path);
       end
       channelstruc(i).hubstrength(j) = count;
       respathwaycount(pathresidues(j)) = respathwaycount(pathresidues(j))+count;
    end
end

writepathwayClusters(PDB,pathstruc,I,Npath,fullfile(char(pathCalcdir) ,'cls'))

writepathwayClusters(PDB,outstruc,1:length(outstruc),...
    length(outstruc), fullfile(char(pathCalcdir) ,'channel'))

% Write pymol script to adjust network appearance
colorlist = {'red' 'green' 'blue' 'cyan' 'magenta' 'yellow'...
    'brown' 'orange'};

% load the structures, set background colors etc.
outstring = { 'load protein_MI.pdb'...
              'show cartoon'...
              'hide lines'...
              'color grey'...
              'set cartoon_cylindrical_helices, 1'...
              'set cartoon_transparency, 0.6'...
              'bg white'...
                      };
for i=1:length(channelstruc)
    outstring(end+1) = {sprintf('load channel%d.pdb',i)};
    if i<=length(colorlist)
        color = colorlist{i};
    else
        color = 'grey40';
    end
    outstring(end+1) = {sprintf('color %s, channel%d',color,i)};
end

for i=1:length(channelstruc)
    outstring(end+1) = {sprintf('load cls%d.pdb',i)};
end
    
outstring(end+1) = {'show sticks, channel*'};
outstring(end+1) = {'hide lines, cls*'};

% Adjust bond thickness according to inter-cluster connectivity
% bondscale = 0.05;
 bondscale = 0.225;
for i=1:length(BondStrength)
    cls = BondStrength(i,1);
    node1 = BondStrength(i,2);
    node2 = BondStrength(i,3);
    strength = BondStrength(i,4);
    for j=node1:node2-1
%         outstring(end+1) = { sprintf('set_bond stick_radius, %f, channel%d & id %d+%d',...
%             bondscale*strength,cls,j,j+1) };
        outstring(end+1) = { sprintf('set_bond stick_radius, %f, channel%d & id %d+%d',...
            bondscale*sqrt(strength),cls,j,j+1) };
    end
end

% spherescale = 0.0625/2.5;
spherescale = 0.15;
if MIWeightPaths % Also weight cut_off
    hubstrength_cutoff = mean([channelstruc.MI])*20;
else
    hubstrength_cutoff = 20;
end
for i=1:length(channelstruc)
    str = '';
    for j=1:length(channelstruc(i).hub)
        resid = num2str(channelstruc(i).hub(j));
        str = [str,resid,'+'];
        outstring(end+1) = { sprintf('set sphere_scale, %f, resi %s & cls%d',...
            spherescale*sqrt(channelstruc(i).hubstrength(j)),resid,i) };
        if channelstruc(i).hubstrength(j)>hubstrength_cutoff
            outstring(end+1) = { sprintf('show spheres, resi %s & protein_MI',...
                resid) };
        end
%     outstring(end+1) = { sprintf('set sphere_scale, %f, resi %s & cls%d',...
%             1,resid,i) };
    end
    str = str(1:end-1);
    outstring(end+1) = { sprintf('show spheres, resi %s & cls%d',str,i) };
    if i<=length(colorlist)
        color = colorlist{i};
    else
        color = 'grey40';
    end
    outstring(end+1) = { sprintf('color %s, resi %s & cls%d',color,str,i) };
end

outstring(end+1) = {'reset'};

FID = fopen(fullfile(pathCalcdir,"pymol.pml"),'w');
for i=1:length(outstring)
    fprintf(FID,'%s\n',    outstring{i});
end
fclose(FID);
