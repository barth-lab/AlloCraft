% calculate how many partners communicate with each dihedral at a given MI
% level

TMres = [];
if isGPCR 
    for i=1:size(helices,1)
        TMres = [TMres helices(i,1):helices(i,2)];
    end
else
    TMres = 0;
end

maxMI = max(MIraw(:));
minMI = min(MIraw(MIraw>0));

MIlist = logspace(log10(minMI),log10(maxMI),500);
MIcum = zeros(1,length(MIlist));
for i=1:length(MIlist)
    MIcum(i) = sum(MIraw(:)>MIlist(i));
end
MIcum = MIcum/(length(MIraw)^2)*100;

% point of inflecion
h = (maxMI-0)/500;
D1 = zeros(1,length(MIlist));
D2 = zeros(1,length(MIlist));
for i=1:length(MIlist)
    if(i==1)
        D1(i) = (MIcum(i+1)-MIcum(i))/h;
    elseif(i==length(MIlist))
        D1(i) = (MIcum(i)-MIcum(i-1))/h;
    else
        D1(i) = (MIcum(i+1)-MIcum(i-1))/2/h;
    end
end
for i=1:length(MIlist)
    if(i==1)
        D2(i) = (D1(i+1)-D1(i))/h;
    elseif(i==length(MIlist))
        D2(i) = (D1(i)-D1(i-1))/h;
    else
        D2(i) = (D1(i+1)-D1(i-1))/2/h;
    end
end

% get backbone & sidechain stats
% Create dihedral list
% get residue list
PDB = pdbread((fullfile(pathCalcdir,"protRenum.pdb")));
atomno = [PDB.Model.Atom(:).AtomSerNo];
resname = cellstr(char(PDB.Model.Atom(:).resName));
resno = [PDB.Model.Atom(:).resSeq];
atomname = char(PDB.Model.Atom(:).AtomName);
CAind = find( strncmp(cellstr(atomname),'CA',3) );
reslist = unique(resno);
atomnameStr = cellstr(atomname);
coordlist = [PDB.Model.Atom(:).X; PDB.Model.Atom(:).Y; PDB.Model.Atom(:).Z]';
for i=1:length(resind)
    dihresname(i) = resname(resind(i));
    if(dihtype(i)==1)
        dihname(i) = cellstr('phi');
    elseif(dihtype(i)==2)
        dihname(i) = cellstr('psi');
    else
        dihname(i) = cellstr('Chi');
    end
end

dihtype(dihtype==2) = 1;
% 1 = bb-bb, 2 = bb-sc, 3 = sc-sc
MIdistall = histc(MIraw(:),MIlist)/length(MIraw(:));
temp1 = repmat(dihtype,1,length(dihtype));
temp2 = temp1';
temp = MIraw(temp1 & temp2);
MIdist1 = histc(temp(:),MIlist)/length(temp(:));
temp = MIraw(temp1+temp2==1);
MIdist2 = histc(temp(:),MIlist)/length(temp(:));
temp = MIraw(~temp1 & ~temp2);
MIdist3 = histc(temp(:),MIlist)/length(temp(:));

figure
semilogx(MIlist,MIdistall,'-k')
hold
semilogx(MIlist,MIdist1,'-r')
semilogx(MIlist,MIdist2,'-','color',[0 0.5 0])
semilogx(MIlist,MIdist3,'-b')
set(gca,'fontsize',18)
xlabel('Inter torsion MI')
ylabel('Population density')
legend('All torsions','BB-BB','BB-SC','SC-SC')
title('Raw MI distribution')

natom = length(coordlist);
dismat = zeros(natom,natom);
for i=1:natom
    temp = coordlist - repmat(coordlist(i,:),natom,1);
    dismat(i,:) = sum(temp.^2,2)';
end
dismat = sqrt(dismat);

% MItable = repmat(' ',length(dihtype)*(length(dihtype)-1)/2,72); % Useless
MItable1 = zeros(length(dihtype)*(length(dihtype)-1)/2,11);
count = 1;
intertorsiondis = zeros(length(MIraw),length(MIraw));
for i=1:length(dihtype)-1
    for j=i+1:length(dihtype)
        if(dihtype(i) && dihtype(j))
            MItype = 1;
        elseif(dihtype(i)+dihtype(j)==1)
            MItype = 2;
        else
            MItype = 3;
        end
        ata1 = dihedral(i,1); ata2 = dihedral(i,2);
        atb1 = dihedral(j,1); atb2 = dihedral(j,2);
%         str1 = sprintf('%s%3d:%s -   %s%3d:%s',...
%             char(resname(ata1)),resno(ata1),(atomname(ata1,:)),...
%             char(resname(ata2)),resno(ata2),(atomname(ata2,:)));
%         str2 = sprintf('%s%3d:%s -   %s%3d:%s',...
%             char(resname(atb1)),resno(atb1),(atomname(atb1,:)),...
%             char(resname(atb2)),resno(atb2),(atomname(atb2,:)));
%         MItable(count,:) = sprintf('%s   %s   %d   %f',str1,str2,MItype,MI(i,j));
        res1 = resno(ata1); res2 = resno(atb1);
        % identify true allosteric MI (ignore pairs within same residue,
        % backbone torsions that are adjacent
        if(res1~=res2)
            if((MItype==1) && (abs(res1-res2)>1))
                isallo = 1;
            elseif(MItype>1)
                isallo = 1;
            else
                isallo = 0;
            end
        else
            isallo = 0;
        end
        %         intertorsiondis(i,j) = mean([dismat(ata1,atb1) dismat(ata1,atb2)...
        %             dismat(ata2,atb1) dismat(ata2,atb2)]);
        intertorsiondis(i,j) = norm(mean(coordlist([ata1,ata2],:))...
            - mean(coordlist([atb1,atb2],:)));
        intertorsiondis(j,i) = intertorsiondis(i,j);
        MItable1(count,:) = [i j ata1 ata2 atb1 atb2 MItype isallo intertorsiondis(i,j) MIraw(i,j) MI(i,j)];
        count = count+1;
    end
end
[~,sortind] = sortrows(MItable1,-10);

% collect allosteric MI
allosortind = sortind(MItable1(sortind,8)>0);

% calculate mean MI as function of distance
mindis = min(dismat(:)); maxdis = max(dismat(:));
dislist = mindis:(maxdis-mindis)/100:maxdis;
MI_dis_pop = zeros(1,length(dislist));
MI_dis_pop_bb = zeros(1,length(dislist));
MI_dis_pop_bs = zeros(1,length(dislist));
MI_dis_pop_ss = zeros(1,length(dislist));
ind = (MItable1(:,10)>0);
ind_bb = ind & (MItable1(:,7)==1);
ind_bs = ind & (MItable1(:,7)==2);
ind_ss = ind & (MItable1(:,7)==3);
ind_TM = (ismember(resind(MItable1(:,1)),TMres) & ismember(resind(MItable1(:,2)),TMres));
MIavg_bb_bb = mean(MItable1(ind_bb,10));
MIavg_bb_sc = mean(MItable1(ind_bs,10));
MIavg_sc_sc = mean(MItable1(ind_ss,10));
MIstd_bb_bb = std(MItable1(ind_bb,10));
MIstd_bb_sc = std(MItable1(ind_bs,10));
MIstd_sc_sc = std(MItable1(ind_ss,10));
% Used to remove average MI for each "type" from filtered MI
MIavg_bb_bb_TM = mean(MItable1(ind_bb & ind_TM,11));
MIavg_bb_sc_TM = mean(MItable1(ind_bs & ind_TM,11));
MIavg_sc_sc_TM = mean(MItable1(ind_ss & ind_TM,11));
MIstd_bb_bb_TM = std(MItable1(ind_bb & ind_TM,11));
MIstd_bb_sc_TM = std(MItable1(ind_bs & ind_TM,11));
MIstd_sc_sc_TM = std(MItable1(ind_ss & ind_TM,11));

MIavg_bb_bb_loop = mean(MItable1(ind_bb & ~ind_TM,11));
MIavg_bb_sc_loop = mean(MItable1(ind_bs & ~ind_TM,11));
MIavg_sc_sc_loop = mean(MItable1(ind_ss & ~ind_TM,11));
MIstd_bb_bb_loop = std(MItable1(ind_bb & ~ind_TM,11));
MIstd_bb_sc_loop = std(MItable1(ind_bs & ~ind_TM,11));
MIstd_sc_sc_loop = std(MItable1(ind_ss & ~ind_TM,11));

temp = MIraw(:);
temp1 = intertorsiondis(:);
for i=1:length(dislist)
    if(i<length(dislist))
        ind_dis = (MItable1(:,9)>=dislist(i)) & (MItable1(:,9)<dislist(i+1));
    else
        ind_dis = (MItable1(:,9)>=dislist(i)) & (MItable1(:,9)<inf);
    end
    MI_dis_pop(i) = mean(MItable1(ind & ind_dis,10));
    MI_dis_pop_bb(i) = mean(MItable1(ind_bb & ind_dis,10));
    MI_dis_pop_bs(i) = mean(MItable1(ind_bs & ind_dis,10));
    MI_dis_pop_ss(i) = mean(MItable1(ind_ss & ind_dis,10));
end
figure
hold
plot(dislist,MI_dis_pop,'-k')
plot(dislist,MI_dis_pop_bb,'-r')
plot(dislist,MI_dis_pop_bs,'-','color',[0 0.5 0])
plot(dislist,MI_dis_pop_ss,'-b')
set(gca,'fontsize',18)
xlabel('Inter-torsion distance,A')
ylabel('Mean MI')
legend('All torsions','BB-BB','BB-SC','SC-SC')
title('Raw MI distance distribution')
