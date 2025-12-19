function writepathways(PDB,pathstruc,K,topcutoff)
clear conectrecord
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
for i=1:topcutoff
    %resi = rescontact(I(i),1);
    %resj = rescontact(I(i),2);
    index = K(i);
    resi = pathstruc(index).path(1);
    resj = pathstruc(index).path(end);
    clsno = pathstruc(index).cls;
    xi = PDB.Model.Atom(CAind(resi)).X;
    yi = PDB.Model.Atom(CAind(resi)).Y;
    zi = PDB.Model.Atom(CAind(resi)).Z;
    xj = PDB.Model.Atom(CAind(resj)).X;
    yj = PDB.Model.Atom(CAind(resj)).Y;
    zj = PDB.Model.Atom(CAind(resj)).Z;
    count = 1;
    count_conect = 1;
    MIconnectFile = ['cls',int2str(clsno),'_','path',int2str(i),'.pdb'];
    FID = fopen(MIconnectFile,'w');
    fprintf(FID,'ATOM  %5d O    NOD %s %3d    %8.3f%8.3f%8.3f  1.00%6.2f\n',...
        count,charset(pathstruc(index).cls),resi,xi,yi,zi,pathstruc(index).MI);
    count = count+1;
    %[dum1,path,dum2] = graphshortestpath(G,resi,resj);
    %[i/2 dum1]
    for j=2:length(pathstruc(index).path)
        xp = PDB.Model.Atom(CAind(pathstruc(index).path(j))).X;
        yp = PDB.Model.Atom(CAind(pathstruc(index).path(j))).Y;
        zp = PDB.Model.Atom(CAind(pathstruc(index).path(j))).Z;
        fprintf(FID,'ATOM  %5d O    PTH %s %3d    %8.3f%8.3f%8.3f  1.00%6.2f\n',...
            count,charset(pathstruc(index).cls),pathstruc(index).path(j),xp,yp,zp,...
            pathstruc(index).MI);
        count = count+1;
        %if j<length(path)
            conectrecord(count_conect) = cellstr( sprintf('CONECT%5d%5d\n',...
                count-2,count-1) );
            count_conect = count_conect+1;
        %end
    end
    fprintf(FID,'ATOM  %5d O    NOD %s %3d    %8.3f%8.3f%8.3f  1.00%6.2f\n',...
        count,charset(pathstruc(index).cls),resj,xj,yj,zj,pathstruc(index).MI);
    for j=1:length(conectrecord)
        fprintf(FID,'%s\n',char(conectrecord(j)));
    end
    fprintf(FID,'END\n');
    fclose(FID);
    count = count+1;
end