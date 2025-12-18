%
% OPEN ISSUES: 
% Preprocessing step: 
%   1- I may be missing important side chain contacts by
%   only taking CAs and CBs, possible solution: consider CG too
%   2- What's the best cutoff for preprocessing? 8 or 10A ?
% 
% Distribution comparison:
%   1- is WS distance enough with the sign addition? Seems like it is for now
%
% Attempting to do something similar to:
% Ensemble Difference Distance Matrix (eDDM) Analysis in Bio3D:
% http://thegrantlab.org/bio3d/articles/online/eddm_vignette/Bio3D_eddm.spin.html
% Grant, Barry J. et al., Protein Science 30.1 (2021).


% Ideas to increase usefulness of information present:
% Option to ignore loops (only changes in structured regions)
% Option to separate BB-BB, SC-BB, and SC-SC distance changes

options.preProc = "CA"; % all, CA, or auto
options.distCutPre  = 15; % in Angstroms, possibly try different distances
options.fracCutSim = 0.4; % Fraction of simulation for contact to be present to be considered significant
options.alpha = 0.05; % Cutoff p value for significance
options.smallwsLimit = 1; % wsd limit, below which values are ignored
options.nBlocks = 10;  % Number of blocks for bootstrapping significance testing of wsd
options.atomChoice = "backbone"; % backbone, sidechain, or heavy

% Add option to choose residues from main and ref systems? (Would be useful
% for when number of residues or starting residues are not identical 
%% First step: Preprocessing: which residue pairs to consider? 

% Another option for preprocessing, look at residues where a min delta(dis)
% occurs no matter their distance

contactMatPre = cell(2,1);
for entryHere = 1:2
    if entryHere == 1
        thisEntry = mainEntry;
        thisChain = mainChain;
    else
        thisEntry = refEntry;
        thisChain = refChain;
    end

    % Option 1: look for distances atom-wise (very slow)
    if strcmp(options.preProc,"all") 
    ['Preprocessing by searching for all heavy atoms within ' num2str(options.distCutPre) 'A']
     tic
    [~, ~, ~, ~, matTot] = proteinContacts(thisEntry.pdb, thisEntry.simulation.concatRuns('StartFrame', settings.frames2skip), ...
        thisChain.resIds, thisChain.name, thisChain.name);
     toc

    % Option 2: look for distances from CA 
    elseif strcmp(options.preProc,"CA")  || (strcmp(options.preProc,"auto") && strcmp(options.atomChoice ,"backbone"))
        ['Preprocessing by searching for CAs within ' num2str(options.distCutPre) 'A']

    indexCA_chainA = selectname(thisEntry.pdb.name,'CA') & selectname(thisEntry.pdb.chainid,thisChain.name);
    
    tic
    concatTraj = thisEntry.simulation.concatRuns('StartFrame', settings.frames2skip);
    nFrames = size(concatTraj,1);
    matTot = zeros(length(thisChain.resIds));
    for frameHere = 1:nFrames
        if mod(frameHere,100) ==0
            frameHere
        end
        [pair,list] = calcpairlist(concatTraj(frameHere,to3(indexCA_chainA)),options.distCutPre);
    
        for i = 1:size(pair,1)
            matTot(pair(i,1),pair(i,2)) =  matTot(pair(i,1),pair(i,2)) + 1;
        end
    end
    
    toc
    matTot = matTot /nFrames;
    %     title(['CA only, ' num2str(options.distCutPre ) ' A'])
    
    % Option 3:
    else  % What if I include CB too? Very similar results to heavy atom, while being
    % 10x faster! 

    % Maybe add CG too? Some side chain contacts are still missing
    if (strcmp(options.preProc,"auto") && strcmp(options.atomChoice ,"sidechain"))
        indexCACB_chainA = selectname(thisEntry.pdb.name, 'CB','CG') & selectname(thisEntry.pdb.chainid,thisChain.name);
        ['Preprocessing by searching for CBs and CGs within ' num2str(options.distCutPre) 'A']
    else % Auto and heavy atom: 
        indexCACB_chainA = selectname(thisEntry.pdb.name,'CA', 'CB','CG') & selectname(thisEntry.pdb.chainid,thisChain.name);
        ['Preprocessing by searching for CAs, CBs and CGs within ' num2str(options.distCutPre) 'A']
    end
    findindexCACB_chainA = find(indexCACB_chainA);
    tic
    concatTraj = thisEntry.simulation.concatRuns('StartFrame', settings.frames2skip);
    nFrames = size(concatTraj,1);
    matTot = zeros(length(thisChain.resIds));
    for frameHere = 1:nFrames
        matHere = zeros(length(thisChain.resIds));
        if mod(frameHere,100) ==0 % Debugging
            frameHere
        end
        [pair,list] = calcpairlist(concatTraj(frameHere,to3(indexCACB_chainA)),options.distCutPre );
    
        for i = 1:size(pair,1)
            res1 = thisEntry.pdb.resseq(findindexCACB_chainA(pair(i,1)));
            res2 = thisEntry.pdb.resseq(findindexCACB_chainA(pair(i,2)));
            
             % Ignore same residue contacts or if contact between res1 and res2
             % in this frame
            if res1 == res2 || matHere(res1,res2)
                continue
            else
              matHere(res1,res2) =  matHere(res1,res2) + 1;
            end
        end
        matTot = matTot + matHere;
    end
    toc
    matTot = matTot /nFrames;

    %      figure; imagesc(matTot>0.8)
    %     title(['CA and CB, ' num2str(options.distCutPre ) ' A'])
    end

    % Ignore residues with +-1 in sequence (-+2 in sequence can still interact
    % via side chains! Ex.: disulfide bond or salt bridging)
    matTot = triu(matTot,2);
    matTot = matTot > options.fracCutSim;  % Apply cutoff

    contactMatPre{entryHere} = matTot;
    clear concatTraj; % Clear HUGE variable (have mercy on the memory!)
end
%% Merge the contact maps from both systems: (union of both)

contactMatMerge = contactMatPre{1} + contactMatPre{2};
[a,b] = find(contactMatMerge);
pairs2Analyze = [a b];

['Preprocessing step done, ' num2str(size(pairs2Analyze,1)) ' pairs have been chosen for further analysis'] 
%% Second step: calculate distances for chosen residues:
tic
minDis = cell(2,1);
for entryHere = 1:2
    if entryHere == 1
        thisEntry = mainEntry;
        thisChain = mainChain;
    else
        thisEntry = refEntry;
        thisChain = refChain;
    end
    
    % Calculate minimum distance between any residue pair
    concatTraj = thisEntry.simulation.concatRuns('StartFrame', settings.frames2skip);
    nFrames = size(concatTraj,1);
    minDis{entryHere} = zeros(nFrames,size(pairs2Analyze,1));
    for i = 1:size(pairs2Analyze,1)
        if mod(i,100) ==0 % Debugging
            i
        end
        pairHere = pairs2Analyze(i,:);
        % Choice of index? backbone? sidechain? Heavy atom?

        if strcmp(options.atomChoice ,"heavy")
            extraSelection = ~selectname(thisEntry.pdb.name,'H*');
        elseif strcmp(options.atomChoice ,"backbone")
            extraSelection = selectname(thisEntry.pdb.name,'CA', 'C', 'N', 'O');
        elseif strcmp(options.atomChoice ,"sidechain")
            extraSelection = ~selectname(thisEntry.pdb.name,'CA', 'C', 'N', 'O') & ...
                ~selectname(thisEntry.pdb.name,'H*');   
            
        end
%         ndx1 = selectid(thisEntry.pdb.resseq,pairHere(1)) & extraSelection;
%         ndx2 = selectid(thisEntry.pdb.resseq,pairHere(2)) & extraSelection;

% This way of looking at resids may be better?
        ndx1 = selectid(thisEntry.pdb.resseq,thisChain.resIds(pairHere(1))) & extraSelection;
        ndx2 = selectid(thisEntry.pdb.resseq,thisChain.resIds(pairHere(2))) & extraSelection;
        % calcMinDis function seems to be working properly, I tested it
        % against VMD via inspection of selected BB-SC and SC-SC
        % interactions
        minDis{entryHere}(:,i) = calcMinDis(concatTraj,ndx1,ndx2); 
    end

end
clear concatTraj; % Clear HUGE variable (have mercy on the memory!)

toc

['Distances between pairs calculated, proceeding with the WS distance calculation now'] 


%% Third step: compare distance distributions between pairs in different systems:


% What metric can I use to compare two given distributions?
% Suggestion: Wasserstein distance:

ws = zeros(size(pairs2Analyze,1),1);
pVal = zeros(size(pairs2Analyze,1),1);

wsMap = zeros(length(refChain.resIds));

tic

for pairi = 1:size(pairs2Analyze,1)
    if mod(pairi,100) ==0 % Debugging
            pairi
    end
    x = minDis{1}(:,pairi);
    y = minDis{2}(:,pairi);

    [wsHere, pVal(pairi)] = ws_distance_pVal(x,y,'nBlocks', options.nBlocks);
    ws(pairi) = sign(mean(x) - mean(y)) * wsHere; % Add sign to wsd
    wsMap(pairs2Analyze(pairi,1),pairs2Analyze(pairi,2)) = ws(pairi);

end

wsFiltered = ws; 
wsFiltered(pVal >= options.alpha) = 0; % Zero out insignificant wsd
wsFiltered(abs(wsFiltered)<options.smallwsLimit) = 0; % zero out small wsd
toc

% Write WS and pairs to spreadsheet
writetable(table(pairs2Analyze(wsFiltered~=0,:),thisChain.resIds(pairs2Analyze(wsFiltered~=0,:)), wsFiltered(wsFiltered~=0), 'VariableNames',{'Pair Ndx','Pair ResID','Filtered WS'}), ...
    fullfile(md2pathdir, ['ws_' settings.mainName '_' settings.refName '_' char(options.atomChoice) '.xls']),'Sheet','Significant WS');

writetable(table(pairs2Analyze, thisChain.resIds(pairs2Analyze), ws, wsFiltered ,pVal, ...
    'VariableNames',{'Pair Ndx','Pair ResID','WS','Filtered WS','p value'}), ...
    fullfile(md2pathdir, ['ws_' settings.mainName '_' settings.refName '_' char(options.atomChoice) '.xls']),'Sheet','All WS');

%% Visualization section:

wsDraw = wsFiltered;
wsDraw(abs(wsDraw)<options.smallwsLimit) = NaN;
figure
s = scatter(pairs2Analyze(:,2),pairs2Analyze(:,1) , 10*abs(wsDraw),wsDraw,'s','filled');

if ~isempty(settings.helices) % Dirty fix to a recent change in Chain.resIds
%     drawTMhelices(pairs2Analyze(:,2), settings.helices,refChain.resIds,'Y');
%     drawTMhelices(pairs2Analyze(:,1), settings.helices,refChain.resIds);
    drawTMhelices(pairs2Analyze(:,2), settings.helices,sort(refChain.resIds),'Y');
    drawTMhelices(pairs2Analyze(:,1), settings.helices,sort(refChain.resIds));
end
% Add color as datatips
dtRows = [dataTipTextRow("Signed WS",wsDraw)];
s.DataTipTemplate.DataTipRows(end+1) = dtRows;

xlabel('Residue')
ylabel('Residue')
colorbar

colormapCustom('div','inputColor1', [0 0.4 0.4] ,'inputColor2', [0.8 0.4 0]);
% Side quest:
matVersion = version('-release'); % Get Matlab version
versionNbr = str2double(matVersion(1:4));
versionLtr = matVersion(5);

if versionNbr >= 2022 
    climTemp = clim(); % Make colormap centered at zero
    clim([-max(abs(climTemp)) max(abs(climTemp))]);
else
    climTemp = caxis(); % Make colormap centered at zero
    caxis([-max(abs(climTemp)) max(abs(climTemp))]);
end

title(['WS distance ' settings.mainName ' - ' settings.refName ', ' char(options.atomChoice) ' atoms'])
set(gca,'FontSize',16)

axis equal
axis square

figPath = fullfile(md2pathdir, "wsdis_eddm_" +  settings.mainName + "_" + settings.refName +options.atomChoice);
savefig( figPath + ".fig");
print2pdf(figPath);

%% Save the WS distance as the B-factor: (Not very useful)

bfactor = zeros(1, mainEntry.atomCount);

if settings.includeLigPathwayCalc % Include ligand in b-factor array
    chainSelect = selectname(mainEntry.pdb.chainid, settings.chains(Chains.receptor)) | selectname(mainEntry.pdb.chainid, settings.chains(Chains.ligand));
    selection = selectname( mainEntry.pdb.name,'CA') & chainSelect;
else
    selection = selectname( mainEntry.pdb.name,'CA') & selectname(mainEntry.pdb.chainid, settings.chains(Chains.receptor));
end
bfactor(selection) = sum(wsMap);

writepdb(fullfile(md2pathdir, "prot_wsDis.pdb"), mainEntry.pdb, [], 'default', bfactor);


%% Experiment with visualization 

% Only plot the TM parts 
isTM = zeros(length(mainChain.resIds),1);

for i = 1:size(settings.helices,1)
    isTM(settings.helices(i,1):settings.helices(i,2)) = 1;
end


wsDraw = wsFiltered;
wsDraw(abs(wsDraw)<options.smallwsLimit) = NaN;
isTMpair = ones(size(ws));
for pairi=1:length(pairs2Analyze)
    if ~isTM(pairs2Analyze(pairi,1)) || ~isTM(pairs2Analyze(pairi,2))
         wsDraw(pairi) = NaN;
        isTMpair(pairi) = 2;
    end
end

figure
s = scatter(pairs2Analyze(:,2),pairs2Analyze(:,1) , 10*abs(wsDraw),wsDraw,'s','filled');
drawTMhelices(pairs2Analyze(:,2), settings.helices,refChain.resIds,'Y');
drawTMhelices(pairs2Analyze(:,1), settings.helices,refChain.resIds);
% Add color as datatips
dtRows = [dataTipTextRow("Signed WS",wsDraw)];
s.DataTipTemplate.DataTipRows(end+1) = dtRows;

xlabel('Residue')
ylabel('Residue')
colorbar

colormapCustom('div','inputColor1', [0 0.4 0.4] ,'inputColor2', [0.8 0.4 0]);

if versionNbr >= 2022 
    climTemp = clim(); % Make colormap centered at zero
    clim([-max(abs(climTemp)) max(abs(climTemp))]);
else
    climTemp = caxis(); % Make colormap centered at zero
    caxis([-max(abs(climTemp)) max(abs(climTemp))]);
end

title(['WS distance ' settings.mainName ' - ' settings.refName ', ' char(options.atomChoice) ' atoms'])
set(gca,'FontSize',16)

axis equal
axis square

%%
% pairi = 203;
% 
% x = minDis{1}(:,pairi);
% y = minDis{2}(:,pairi);
% 
% figure;histogram(x,50)
% 
% hold on
% histogram(y,50)
% 
% % ws = ws_distance(x,y);
% 
% 
% title(['Pair number: ' num2str(pairi) ', ResPair: ' num2str(pairs2Analyze(pairi,:)) ', ws = ' num2str(ws(pairi))])
% 
% counts1 = histcounts(minDis{1}(:,1),50,'Normalization','pdf');
% counts2 = histcounts(minDis{2}(:,1),50,'Normalization','pdf');
% % wsp = ws_distance(counts1,counts2)
% % ws = ws_distance(x,y)
% % signrank(counts1,counts2)
% 
% %%
% temptrj = trj(:, to3(ndx1));
% tempMat  = reshape(trj(:, to3(ndx2)),nFrames,3,[]);
% 
% tempDis = zeros(nFrames,sum(ndx1),sum(ndx2));
% for natom = 1:(size(temptrj,2)/3)
%    tempDis(:,natom,:) = sqrt(sum(( temptrj(:,(3*(natom-1)+1):3*natom) -tempMat).^2, 2));
% end
% % Take the minimum along dimensions 2 (ndx1) and 3 (ndx2)
% temp1 = min(tempDis,[],2);
% minDis = min(temp1,[],3);
% 
% %%
% 
% X = rand(3,3);
% Y = rand(4,3);
% 
% 
% D1 = pdist2(X,Y);
% 
% % 
% %     0.5417    0.1210    0.1137    0.5836
% %     0.8219    0.5300    0.6283    1.0570
% %     0.6135    0.3236    0.3861    0.7670
% 
% D2 = vecnorm(X - Y,2,2);
%  %%
% 
%  test1 = 2+0.5*(randn(1000,1));
%  test2 = 4+0.5*(randn(1000,1));
% 
% 
% 
%  figure
%  subplot(1,2,1)
%  plot(test1)
%  hold on
%  plot(test2)
% 
%  plot(test2 - test1)
%  legend('1','2','DIff')
% 
%  subplot(1,2,2)
%  histogram(test1)
%  hold on 
%  histogram(test2)
% 
%   histogram(test2 - test1)
% 
%  [p1,h1] = signrank(test2 - test1,1)
%  [p2,h2] = signrank(test1,test2)
% 
