% Father directory for all the runs
metadir = 'D:\Telework_library\dopamine_phase_3';
databasePath = 'C:\Users\mahdi\Documents\md2pathDatabase\'; 

pdbName = 'prot.pdb';
% foldersToStudy = {'3-d2_bromo_WT','4-d2_bromo_T174M-C220L','15-d2_brc_F217M','17-d2_brc_I125N','36a-d2_brc_barr2_F217M','35-d2_brc_barr2_WT_forREAL'};
% numRunsSys = [10 5 5 5 5 4]; % number of runs for each system

% foldersToStudy = {'1-d2_dop_WT','2-d2_dop_T174M-C220L','14-d2_dop_F217M','16-d2_dop_I125N', ...
%     '36b-d2_dpa_barr2_F217M','37b-d2_dpa_barr2_I125N','34-d2_dpa_barr2_WT_for_Real'};

foldersToStudy = {'34-d2_dpa_barr2_WT_for_Real', ...
    '36b-d2_dpa_barr2_F217M','37b-d2_dpa_barr2_I125N'};
numRunsSys = [4 5 3]; % number of runs for each system
% thisSysLabel = {'WT','T5.54M-C6.47L','F6.44M','I4.46N','F6.44M barr','I4.46N barr', 'WT barr'};

thisSysLabel = {'WT barr','F6.44M barr','I4.46N barr'};

md2pathName = 'md2pathdev';
alloPathCalcName = 'alloPathCalc_Culled_data';
refNdx = 1; % Which system to compare to?

name = 'DD2R-barr2';
disCutoff = 10;

%% Make the meta-direcory
if isunix
    slash = '/';
    copy = 'cp';
elseif ispc
    slash = '\';
    copy = 'copy';
end
md2pathdirMeta = [metadir slash 'md2pathMeta' slash];
md2pathdir = [metadir slash 'md2pathMeta' slash name slash];

% Create meta directory
if ~exist(md2pathdirMeta, 'dir')
    mkdir(md2pathdirMeta)
    add2log(md2pathdirMeta, ['Creating md2path directory in ' metadir]);
end

% Create subdirectories for every bundle of systems studied
if ~exist(md2pathdir, 'dir')
    mkdir(md2pathdir)
    add2log(md2pathdir, ['Creating md2path directory in ' md2pathdirMeta]);
end

%% Load, fetch and align PDB files

metaLoadAlign
%% Plot average MI vs distance
figure
MIdis = cell(length(foldersToStudy),1);
for thisSys = 1:length(foldersToStudy)
%     mydir = ([metadir slash foldersToStudy{thisSys} slash md2pathName slash alloPathCalcName]); 
    mydir = fullfile(metadir,foldersToStudy{thisSys},md2pathName,alloPathCalcName);
    MIdis{thisSys} = importdata(fullfile(mydir,"MI-dis.dat"));
    plot(MIdis{thisSys}(:,1),MIdis{thisSys}(:,2),'LineWidth',2 - (thisSys/5))
%     scatter(MIdis{thisSys}(:,1),MIdis{thisSys}(:,2),10,'filled')
    hold on

end
legend(thisSysLabel)
title(name,'FontSize',20)
%% Calculate MIres, check effect of MI filtering:
MINames = {'MIraw','Excess MI','Average filtered MI'};
MIresCell = cell(length(foldersToStudy),1);

for thisSys = 1:length(foldersToStudy)
%     mydir = ([metadir slash foldersToStudy{thisSys} slash md2pathName slash alloPathCalcName]); 
    mydir = fullfile(metadir,foldersToStudy{thisSys},md2pathName,alloPathCalcName);

    reSortRenum = load(fullfile(mydir,"reSort.txt"));
    resind = reSortRenum(:,1);
    dihtype = reSortRenum(:,2);
    dihedral = reSortRenum(:,4:5);
    Nres = max(resind);
    
    MIraw = load(fullfile(mydir,"MI_.mat"));
    MIfiltered = load(fullfile(mydir,"MI_filtered.mat"));
    MICell = {MIraw.I, MIfiltered.excessMi};
    
    figure('Position', [10 100 1600 400])
    tiledlayout(1,2)


    for k=1:2
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
        colorbar
        colormap turbo

    %     if isGPCR % Add more general secondary structure?
    %         drawTMhelices(1:Nres,settings.helices,1:Nres)
    %         drawTMhelices(1:Nres,settings.helices,1:Nres,'Y')
    %     end
    end
    sgtitle(thisSysLabel{thisSys})
    MIresCell{thisSys} = MIresHere{2};

end

%% Distance filter MIres:
MIresDisCutCell = cell(length(foldersToStudy),1);

for thisSys = 1:length(foldersToStudy)
    MIres = MIresCell{thisSys};
%     mydir = ([metadir slash foldersToStudy{thisSys} slash md2pathName slash alloPathCalcName]); 
    mydir = fullfile(metadir,foldersToStudy{thisSys},md2pathName,alloPathCalcName); 
    [pdbHere,crdHere] = readpdb(fullfile(mydir,'protein_MI.pdb'));
    resList = unique(pdbHere.resseq);
    Nres = length(resList);
    ndxCA = selectname(pdbHere.name,'CA');
    disMat = calcdistancematrix(crdHere(to3(ndxCA))); % CA distance matrix

    MIresDisCut = zeros(size(MIres));

    for i=1:Nres-1
        for j=i+1:Nres
            if (disMat(i,j)<=disCutoff)
                MIresDisCut(i,j) = MIresCell{thisSys}(i,j);
    %             MIresDisCut(i,j) = MIresHere{2}(i,j);
                MIresDisCut(j,i) =  MIresDisCut(i,j);
            end
        end
    end
    MIresDisCutCell{thisSys} = MIresDisCut;
end

%% Visualize difference in MI (with reference state) via graph representation:


 
for thisSys = 3 %:length(foldersToStudy)
%     mydir = ([metadir slash foldersToStudy{thisSys} slash md2pathName slash alloPathCalcName]); 
    mydir = fullfile(metadir,foldersToStudy{thisSys},md2pathName,alloPathCalcName); 
    [pdbHere,crdHere] = readpdb(fullfile(mydir,'protein_MI.pdb'));
    CAcoord = pdbHere.xyz(ndxCA,:);

    GmatDiff = MIresDisCutCell{thisSys} - MIresDisCutCell{refNdx};
    GmatDiff(abs(GmatDiff)<0.5) =0;
    
    Gdiff = graph(GmatDiff);
    
    figure
    pHere = plot(Gdiff,'XData',CAcoord(:,1),'YData',CAcoord(:,2),'ZData',CAcoord(:,3),'MarkerSize',5);
    xlabel('X [A]');
    ylabel('Y [A]');
    zlabel('Z [A]');
    axis('equal')
    
%       pHere.NodeCData = [1:Nres];
    % Color edges by MI difference
     pHere.EdgeCData = Gdiff.Edges.Weight; 
    % Make colorbar symmetric
    caxisTemp = caxis;
    caxis([- max(abs(caxisTemp))  max(abs(caxisTemp))]);

     colormapCustom('div','inputColor1', [0 0.4 0.4] ,'inputColor2', [0.8 0.4 0]);
     colorbar
    title([thisSysLabel{thisSys} ' vs ' thisSysLabel{refNdx}])
end


%

% Remember: it's test - WT
MIDiffCut = 5;
resDiff = sum(GmatDiff);
[a,b] = sort(resDiff);   

% MI gains
highlight(pHere, b(a>MIDiffCut) ,'NodeColor',[0.8 0.4 0],'MarkerSize',10)

highlight(pHere, b(a< - MIDiffCut),'NodeColor', [0 0.4 0.4],'MarkerSize',10)


refChain.formatResidues(b(a>MIDiffCut),'SecondaryEntry',database.entries{3})
refChain.formatResidues( b(a< - MIDiffCut),'SecondaryEntry',database.entries{3})

 % Add the legend in a sneaky way
h = zeros(2, 1);
hold on
h(1) = scatter(NaN,NaN,150,'MarkerEdgeColor',[0.8 0.4 0],...
              'MarkerFaceColor',[0.8 0.4 0], ...
              'LineWidth',1.5);
h(2) = scatter(NaN,NaN,150,'MarkerEdgeColor',[0 0.4 0.4],...
              'MarkerFaceColor',[0 0.4 0.4], ...
              'LineWidth',1.5);

legend(h, [thisSysLabel{thisSys} ' focused residues'], [thisSysLabel{refNdx} ' focused residues'])

%  set(gca,'zdir','reverse')
%  set(gca,'ydir','reverse')
figName = md2pathdir;
%%

resFocus = Gdiff.Edges( logical(sum(Gdiff.Edges.EndNodes ==253,2)),:);

refChain.formatResidues(resFocus.EndNodes(:,1))
refChain.formatResidues(resFocus.EndNodes(:,2))


resFocus.Weight


%% Compare MI between different domains (needs domain definition, of course)

