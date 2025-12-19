%% Load shit in:
md2pathdir = 'D:\Telework_library\dopamine_phase_3\35-d2_brc_barr2_WT_forREAL\md2pathdev\';
refdir = 'D:\Telework_library\dopamine_phase_3\34-d2_dpa_barr2_WT_for_Real\md2pathdev\';
name = cell(2,1);
tempStruc = load(fullfile(md2pathdir,'workspace.mat'));
rmsf1 = tempStruc.RMSF;
name{1} = 'D2-BRC-WT-barr2';

tempStruc = load(fullfile(refdir,'workspace.mat'));
rmsf2 = tempStruc.RMSF;
name{2} = 'D2-DA-WT-barr2';


%% Plot main figure
errStrideRMSF = 3;
figure
RMSF_mean_cell = cell(2,1);
for i = 1:2
    if i ==1
        RMSF = rmsf1;
    else
%         RMSF = rmsf_6vms(3:end,:);
        RMSF = rmsf2;
    end
    RMSF_mean_CA = mean(RMSF,2);
    RMSF_std_CA = std(RMSF,0,2);
    
    plot(RMSF_mean_CA,'LineWidth',1)
    hold on
    errorbar(1:errStrideRMSF:length(RMSF_mean_CA),RMSF_mean_CA(1:errStrideRMSF:length(RMSF_mean_CA)),RMSF_std_CA(1:errStrideRMSF:length(RMSF_mean_CA))/sqrt(size(RMSF,2)) ...
    ,'o','MarkerSize',5)
    legend_entries{2*i-1} =name{i};
    legend_entries{2*i} = '';

    RMSF_mean_cell{i} = RMSF_mean_CA;
end
% Draw helices if helices are defined
if ~isempty(settings.helices)
    drawTMhelices(RMSF_mean_CA, settings.helices, mainChain.resIds)
end
xlabel('Residue'); ylabel('RMSF [Angstrom]')

legend(legend_entries,'location','best');
legend boxoff
title('RMSF of C\alpha atoms - Receptor')
set(gca,'FontSize',16)
%% Add ligand binding residues:

bsRes.main = importdata(fullfile(md2pathdir,'BS_residues.txt'));
bsRes.ref = importdata(fullfile(refdir,'BS_residues.txt'));

hold on
% scatter(bsRes.main,RMSF_mean_CA(receptorResIdsNdx(bsRes.main)), 'LineWidth',1.5)
scatter(bsRes.main,RMSF_mean_cell{1}(bsRes.main),'k', 'LineWidth',1.5,'DisplayName',[name{1} ' binding'])

scatter(bsRes.ref,RMSF_mean_cell{1}(bsRes.ref), 25,'r', 'filled','DisplayName',[name{2} ' binding'])
%% On another note: MIperresidue?
alloPathName = 'alloPathCalc_Culled_data'; % alloPathCalc or alloPathCalc_Culled_data

miPerRes1 = importdata(fullfile(md2pathdir,alloPathName,'MIperresidue.txt'));
miPerRes2 = importdata(fullfile(refdir,alloPathName,'MIperresidue.txt'));
% 3rd system too?
% extradir = 'D:\Telework_library\dopamine_phase_3\21-d2_ris_6cm4\md2pathdev\';
% name{3} = 'RIS-6CM4';

figure; 
s1 = scatter( mainChain.resIds,miPerRes1,25,'^','filled'); hold on; 
s2 = scatter( mainChain.resIds,miPerRes2,25,'filled');
resMain = mainEntry.chains{1}.resIds; % chain 1 is always main chain
resText = mainEntry.chains{1}.formatResidues(resMain,'BWonly',true);
row = dataTipTextRow('Residue',resText);
s1.DataTipTemplate.DataTipRows(end+1) = row;
s2.DataTipTemplate.DataTipRows(end+1) = row;

if ~isempty(extradir)
    miPerRes3 = importdata(fullfile(extradir,alloPathName,'MIperresidue.txt'));
    scatter( mainChain.resIds,miPerRes3,25,'v','filled')
end

drawTMhelices(miPerRes1, settings.helices, mainChain.resIds)

xlabel('Residue i'); ylabel('\Sigma_j(MI(i,j)/Nres')

legend(name,'location','best');
legend boxoff
title('Average summed MI per residue')

set(gca,'FontSize',16)

% scatter(settings.highlightRes, withPam(settings.highlightRes),60,'x', 'LineWidth',1.5);
% legend('Without PAM','With PAM',settings.highlightText)
% 
