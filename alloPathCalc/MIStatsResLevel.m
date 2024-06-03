%% Find residue distance, average MI at far/near distances, and summed MI per residue
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

% For a given residue, find average MI of other residues at close and far distances
MIfarnearratio = zeros(length(reslist),1);
MIperresidue = zeros(length(reslist),1);
for i=1:Nres
    avgMInear = 0;
    avgMIfar = 0;
    countnear = 0;
    countfar = 0;
    for j=1:Nres
        if dismat(i,j)<=disCutOff % Defaults to 10 A
            avgMInear = avgMInear+MIres(i,j);
            countnear = countnear+1;
        else
            avgMIfar = avgMIfar+MIres(i,j);
            countfar = countfar+1;
        end
        MIperresidue(i) = MIperresidue(i) + MIres(i,j);
    end
    avgMInear = avgMInear/countnear;
    avgMIfar = avgMIfar/countfar;
    if(avgMInear==0)
        avgMInear = 1;
    end
    MIfarnearratio(i) = avgMIfar/avgMInear;
end
MIperresidue = MIperresidue/length(reslist);

% MIfarnearratio and MIperresidue are not used beyond this point!!!
% Although they are kinda useful stats (KINDA)
writematrix(MIfarnearratio,fullfile(pathCalcdir,"MIfarnearratio.txt"));
writematrix(MIperresidue,fullfile(pathCalcdir,"MIperresidue.txt"));


 %% Now calculate MI distance distribution
clear disx MIx MIdist
disx = 0:1:ceil( max(dismat,[],'all'));
MIx = min(min(MIres)):0.01:max(max(MIres));
MIresavg = zeros(1,length(disx));
count = zeros(1,length(disx));

MIdist = zeros(length(disx),length(MIx));
for i=1:length(disx)-1
    MIvalues = MIres( (dismat>=disx(i)) & (dismat<disx(i+1)) );
    for j=1:length(MIx)-1
        MIdist(i,j) = sum( (MIvalues>=MIx(j)) & (MIvalues<MIx(j+1)) );
    end
    value = mean(MIvalues);
    if ~isnan(value)
        MIresavg(i) = value;
    end
end

DP = sum(MIdist,2)/sum(sum(MIdist));
MP = sum(MIdist,1)/sum(sum(MIdist));
MIdist = MIdist/sum(sum(MIdist));

MIdistnorm = zeros(length(disx),length(MIx));
for i=1:length(disx)
    for j=1:length(MIx)
        if DP(i) && MP(j)
            MIdistnorm(i,j) = MIdist(i,j)/DP(i)/MP(j);
        end
    end
end
writematrix([disx' MIresavg'],fullfile(pathCalcdir,"MI-dis.dat"),'Delimiter','tab')


figure; plot(disx,MIresavg,'LineWidth',1)
xlabel('Interresidue distance [A]')
ylabel('Average interresidue MI');
set(gca,'FontSize',16)
