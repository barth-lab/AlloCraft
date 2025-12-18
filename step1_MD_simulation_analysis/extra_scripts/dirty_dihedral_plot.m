%% Load KLdiv data first!
% Temp dihedrals, please delete later
dihedralsRef = refSim.dihedralsMat(:,mainSim.klDivStuff.reSortndxRef);
dihedralsTest = mainSim.dihedralsMat(:,mainSim.klDivStuff.reSortndxHere);

dihCount = size(dihedralsRef, 2);

xDih = 1:dihCount;

% Just plot dihedrals
figure('Renderer', 'painters', 'Position', [10 10 1500 750]) ;
tiledlayout('flow')

res2study = [ 239 268];
lowKLDivLimit = 0; 
skipBB = true;

for i =1:length(res2study)
    resHere = res2study(i);
    % Get dihedrals belonging to residue:
    ndxHere = find(reSortCommon(:,1) == resHere);
    ndxRef = find(reSortCommonRef(:,1) == resHere);
    for j = 1:length(ndxHere)
        if kl1(ndxHere(j))< lowKLDivLimit % Skip low KLDiv dihedrals
            continue
        end
        dihType = reSortCommon(ndxHere(j),2);
        if skipBB && (dihType==1 || dihType==2)
            continue
        end
        nexttile
        binLimitsHere = binLimits(dihedralsRef(:,ndxRef(j)),dihedralsTest(:,ndxHere(j)));
        
        histogram(dihedralsTest(:,ndxHere(j)),50, 'BinLimits', binLimitsHere, 'Normalization', 'pdf')
        hold on 
        histogram(dihedralsRef(:,ndxRef(j)),50, 'BinLimits', binLimitsHere, 'Normalization','pdf')

        % Find the resname and dihType:
%         dihType = reSortCommon(ndxHere(j),2);
        if dihType == 1
            dihText = '\phi';
        elseif dihType == 2
            dihText = '\psi';
        else
            dihText = '\chi';
        end
    
        resMain = (mainEntry.pdb.resname(mainEntry.pdb.resseq==reSortCommon(ndxHere(j),1),:));
        resRef = (refEntry.pdb.resname(refEntry.pdb.resseq==reSortCommonRef(ndxRef(j),1),:));
        textMain = ['Test: ' resMain(1,:) num2str(reSortCommon(ndxHere(j),1)) ' ' dihText];
        textRef = ['Ref: ' resRef(1,:) num2str(reSortCommonRef(ndxRef(j),1)) ' ' dihText];
        legend(textMain,textRef,'Location','best')
        legend boxoff
        title(['KL1 = ' num2str(kl1(ndxHere(j)))])
    end
end
sgtitle(['Dihedral dist. of ' mainEntry.name '|' refEntry.name])

% savefig(['D1_dihedrals_' mainEntry.name '_' refEntry.name '_res' num2str(res2study) '.fig'])
%% Support functions:

function output = binLimits(xref,x) 
    % Creates common bin limits for two dihderal series while keeping them
    % in the [ 0 2*pi] limit
    minCommon = min(min(x), min(xref));
    maxCommon = max(max(x), max(xref));
    output = [max(minCommon, 0), min(maxCommon, 2*pi)];

end

