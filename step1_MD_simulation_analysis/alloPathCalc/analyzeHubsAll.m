% Identify all residues that are involved in pathways
clear allhubsFreq allhubsFreqAll
allhubs = [pathstruc(I(1:Npath)).path];
allhubsMI = zeros(size(allhubs));
count = 1;
for i=1:Npath % Grab corresponding path pair MI for every hub in allhubs
    allhubsMI(count:count+pathstruc(I(i)).Npath-1) = pathstruc(I(i)).MI;
    count = count + pathstruc(I(i)).Npath;
end

allhubsList = unique(allhubs);
allhubsFreq = zeros(size(allhubsList));
for i=1:length(allhubsList)
    if MIWeightPaths
        ndxHere = allhubs==allhubsList(i);
        allhubsFreq(i) = sum(allhubsMI(ndxHere)); % scale counted pathways with pair MI
    else
        allhubsFreq(i) = sum(allhubs==allhubsList(i)); % count pathways
    end
end
allhubsData = [allhubsList; allhubsFreq]';
allhubsData = sortrows(allhubsData,-2);

writematrix(allhubsData,fullfile(pathCalcdir,"allhubResidues.txt"),'Delimiter','tab')

% Plot global hubs ion the 3D structure:
figure;
allohubs = allhubsData(1:20,1);
alloscores = allhubsData(1:20,2);

p3 = plot3(CAcoord(:,1),CAcoord(:,2),CAcoord(:,3),...
    'color',[0.8 0.5 0.],'Linewidth',2);

if exist("mainEntry",'var')
    resMain = mainEntry.chains{1}.resIds; % chain 1 is always main chain
    if includeLigPathwayCalc % Make list from both main Chain and ligand chain
        resLig = mainEntry.chains{3}.resIds;  % chain 3 is always ligand
        resText = [mainEntry.chains{1}.formatResidues(resMain,'BWonly',true); ...
            mainEntry.chains{3}.formatResidues(resLig)];
    else
        resText = mainEntry.chains{1}.formatResidues(resMain,'BWonly',true);
    end
    row = dataTipTextRow('Residue',resText);
else
    row = dataTipTextRow('Residue',1:Nres);
end
p3.DataTipTemplate.DataTipRows(end+1) = row;
axis('equal')
hold

s3 = scatter3(CAcoord(allohubs,1),CAcoord(allohubs,2),...
    CAcoord(allohubs,3),alloscores,alloscores,'filled');
hcb = colorbar;
hcb.Label.FontSize = 16;
hcb.Label.String = "Hubscore";

row = dataTipTextRow('Residue',allohubs);
s3.DataTipTemplate.DataTipRows(end+1) = row;
hold off
alpha(.5)
title('Top 20 global allosteric hubs','FontSize',16)

% Continue with GPI and BAI hubs:
allhubsFreqAll = zeros(Nres,1);
allhubsFreqAll(allhubsList) = allhubsFreq;

GPIhubs = zeros(Nres,1);
GPIhubsFreq = zeros(Nres,1);
% GPIres = dlmread('GPI_residues.txt');
GPIres1 = GPIres(GPIres>0);
BAIhubs = zeros(Nres,1);
BAIhubsFreq = zeros(Nres,1);
% BAIres = dlmread('BAI_residues.txt');
BAIres1 = BAIres(BAIres>0);
for i=1:Npath
    path = pathstruc(I(i)).path;
    if(sum(ismember(path,GPIres1))>0)       
        if MIWeightPaths
            GPIhubsFreq(path) = GPIhubsFreq(path) + pathstruc(I(i)).MI;
        else
            GPIhubsFreq(path) = GPIhubsFreq(path)+1; % count pathways
        end
    end
    if(sum(ismember(path,BAIres1))>0)    
        if MIWeightPaths
            BAIhubsFreq(path) = BAIhubsFreq(path) + pathstruc(I(i)).MI; 
        else
            BAIhubsFreq(path) = BAIhubsFreq(path)+1; % count pathways
        end
    end
end
writematrix([GPIhubsFreq BAIhubsFreq],fullfile(pathCalcdir,"allhubResiduesAll.txt"),'Delimiter','tab')
GPIhubs = find(GPIhubsFreq>0);
GPIhubsFreq = GPIhubsFreq(GPIhubsFreq>0);
GPIhubsMI = MIres(GPIhubs,GPIres1);
MIcutoff = pathstruc(I(Npath)).MI;
GPIhubsMI(GPIhubsMI<MIcutoff) = 0;
temp = zeros(length(GPIhubs),length(GPIres));
temp(:,GPIres>0) = GPIhubsMI;
GPIhubsMI = temp;

GPIhubsdata = [GPIhubs GPIhubsFreq GPIhubsMI];
writematrix(GPIhubsdata,fullfile(pathCalcdir,"GPIhubs.txt"),'Delimiter','tab')
