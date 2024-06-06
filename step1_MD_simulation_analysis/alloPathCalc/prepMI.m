randMIlim = 0; % This is obsolete with the new filterCorrectMI function (Thanks Simon!)
reSortRenum = load(fullfile(pathCalcdir,"reSort.txt"));
resind = reSortRenum(:,1);
dihtype = reSortRenum(:,2);
dihedral = reSortRenum(:,4:5);
excluded = load(fullfile(pathCalcdir,"excluded_res.txt"));

if isGPCR % Add more general secondary structure?
    % read helices definitions
    helices = load(fullfile(pathCalcdir,"helices.txt"));
end

excluded_res = [];
for i=1:size(excluded,1)
    excluded_res = [excluded_res excluded(i,1):excluded(i,2)];
end

MI = load(miName,'I');
MI = MI.I;
MI(MI<0) = 0;
MI(isnan(MI)) = 0;
MIraw = MI; % Save unfiltered MI for comparison with filtered



%% Filter MI according to a statistical significance criterion:

load(dihedralsPath);
% load(fullfile(pathCalcdir,'dihedrals.mat'));
MI = filterCorrectMI(dihedralsMat, MIraw, 'SaveDir', pathCalcdir,'SignificanceThreshold',0.05);
% [MI, binMiArray, iiArray, significanceArray]= filterCorrectMI(dihedralsMat, MIraw, 'SaveDir', pathCalcdir,'SignificanceThreshold',0.05);
%% Calculate MI stats, distribution vs distance and vs dihedral type

MIAnalysisDihLevel;

removeAverageMI;
% Outputs MIfilterAvg
%% Beyond this point, only MIres is used

%% Plot residue-wise MI:

MICell = {MIraw, MI, MIfilterAvg};
MINames = {'MIraw','Excess MI','Average filtered MI'};
Nres = length(reslist);
MIresHere = cell(3,1);
% compute residue summed MI

figure('Position', [10 100 1600 400])

tiledlayout(1,3)

for k=1:3
    MIhere = MICell{k};
    MIresHere{k} = zeros(Nres,Nres);

    for i=1:Nres-1
        for j=i+1:Nres
            temp = MIhere(resind==(i),resind==(j));
            temp = temp(:);
            MIresHere{k}(i,j) = sum(temp);
            MIresHere{k}(j,i) = MIresHere{k}(i,j);
        end
    end
    nexttile
    imagesc(MIresHere{k})
    if k==1
        ylabel('Residue')
    end
    xlabel('Residue')
    title(MINames{k})
    hcb = colorbar;
    hcb.Title.String = "KT";
    set(gca,'FontSize',16)
    if isGPCR % Add more general secondary structure?
        drawTMhelices(1:Nres,helices,1:Nres)
        drawTMhelices(1:Nres,helices,1:Nres,'Y')
    end

end
colormap turbo


MIres = MIresHere{3}; % Calculated from MIfilterAvg
clear MICell;
savefig(fullfile(pathCalcdir,"MIres"))
save(fullfile(pathCalcdir,"workspace.mat"),'MIres','-append');


% 
% Nres = length(reslist);
% MIres = zeros(Nres,Nres);
% % compute residue summed MI
% 
% for i=1:Nres-1
%     for j=i+1:Nres
%         temp = MIfilterAvg(resind==(i),resind==(j));
%         temp = temp(:);
%         MIres(i,j) = sum(temp);
%         MIres(j,i) = MIres(i,j);
%     end
% end


% Clear excluded residues by zeroing out their MI 
for i=1:Nres
    for j=1:Nres
        if ismember(i,excluded_res) || ismember(j,excluded_res)
            MIres(i,j) = 0;
        end
    end
end
