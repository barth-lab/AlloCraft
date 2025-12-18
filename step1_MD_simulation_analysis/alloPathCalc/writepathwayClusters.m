function writepathwayClusters(PDB,pathstruc,K,topcutoff,prefix)
clear conectrecord
ncls = max( [pathstruc(K(1:topcutoff)).cls] );
atomno = [PDB.Model.Atom(:).AtomSerNo];
resno = [PDB.Model.Atom(:).resSeq];
atomname = char(PDB.Model.Atom(:).AtomName);
CAind = find( strncmp(cellstr(atomname),'CA',3) );
reslist = unique(resno);
Nres = length(reslist);
%charset = repmat(['A':'Z','a':'z','0':'9'],1,Nres*Nres);
charset = [['A':'Z','0':'9'],repmat(' ',1,Nres*Nres)];
CAcoord = [PDB.Model.Atom(CAind).X; PDB.Model.Atom(CAind).Y;...
    PDB.Model.Atom(CAind).Z]';
for i=1:ncls
    members = find([pathstruc(K(1:topcutoff)).cls] == i);
    count = 1;
    count_conect = 1;
    MIconnectFile = [prefix,int2str(i),'.pdb'];
    FID = fopen(MIconnectFile,'w');
    clear conectrecord
    for k=1:length(members)
        index = K(members(k));
        resi = pathstruc(index).path(1);
        resj = pathstruc(index).path(end);
        clsno = pathstruc(index).cls;
        xi = PDB.Model.Atom(CAind(resi)).X;
        yi = PDB.Model.Atom(CAind(resi)).Y;
        zi = PDB.Model.Atom(CAind(resi)).Z;
        xj = PDB.Model.Atom(CAind(resj)).X;
        yj = PDB.Model.Atom(CAind(resj)).Y;
        zj = PDB.Model.Atom(CAind(resj)).Z;
        
        fprintf(FID,'ATOM  %5d O    NOD %s %3d    %8.3f%8.3f%8.3f  1.00%6.2f\n',...
            count,charset(pathstruc(index).cls),resi,xi,yi,zi...
            ,pathstruc(index).MI);
        count = count+1;
        %[dum1,path,dum2] = graphshortestpath(G,resi,resj);
        %[i/2 dum1]
        for j=2:length(pathstruc(index).path)-1
            xp = PDB.Model.Atom(CAind(pathstruc(index).path(j))).X;
            yp = PDB.Model.Atom(CAind(pathstruc(index).path(j))).Y;
            zp = PDB.Model.Atom(CAind(pathstruc(index).path(j))).Z;
            fprintf(FID,'ATOM  %5d O    PTH %s %3d    %8.3f%8.3f%8.3f  1.00%6.2f\n',...
                count,charset(pathstruc(index).cls),pathstruc(index).path(j),xp,yp,zp,...
                pathstruc(index).MI);
            %if j<length(path)
                conectrecord(count_conect) = cellstr( sprintf('CONECT%5d%5d\n',...
                    count-1,count) );
                count = count+1;
                count_conect = count_conect+1;
            %end
        end
        conectrecord(count_conect) = cellstr( sprintf('CONECT%5d%5d\n',...
                count-1,count) );
        fprintf(FID,'ATOM  %5d O    NOD %s %3d    %8.3f%8.3f%8.3f  1.00%6.2f\n',...
            count,charset(pathstruc(index).cls),resj,xj,yj,zj,pathstruc(index).MI);
        count = count+1;
        count_conect = count_conect+1;
    end
    for j=1:length(conectrecord)
        fprintf(FID,'%s\n',char(conectrecord(j)));
    end
    fprintf(FID,'END\n');
    fclose(FID);
end