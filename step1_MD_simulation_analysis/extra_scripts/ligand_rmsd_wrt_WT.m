% DA
% load(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\dirty matlab scripts','D2_DA_ligand_RMSD.mat'));

rmsd_all = calcrmsd(tempLigTraj(1,:),tempLigTraj);


rmsd_stats = zeros(length(foldersToStudy),2);
rmsf_stats = zeros(length(foldersToStudy),2);

figure
for thisSys = 1:length(foldersToStudy)
    numRuns = numRunsSys(thisSys);
    rmsd_here = rmsd_all(CSys(:,1)==thisSys);
    plot(rmsd_here);
    hold on
    rmsd_stats(thisSys,:) = [mean(rmsd_here) std(rmsd_here)/sqrt(numRuns)];
end
xlabel('Frame')
ylabel('RMSD [A]')
legend(thisSysLabel)
title('RSMD to DA WT input model')
set(gca,'FontSize',20)

% RMSF
figure
for thisSys = 1:length(foldersToStudy)
    numRuns = numRunsSys(thisSys);
    rmsf_here = calcrmsf(tempLigTraj(CSys(:,1)==thisSys,:),[],numRuns);
    errorbar(1:(size(tempLigTraj,2)/3), rmsf_here(:, 1), rmsf_here(:, 2));
    hold on
    rmsf_stats(thisSys,:) = [mean(rmsf_here(:, 1)) mean(rmsf_here(:, 2))];
end
xlabel('Frame')
ylabel('RMSF [A]')
legend(thisSysLabel)
title('RSMF of DA')
set(gca,'FontSize',20)

tab1 = table(thisSysLabel',rmsd_stats,rmsf_stats);

%% BRC
load(fullfile('D:\Telework_library\dopamine_phase_3\a-analysis\dirty matlab scripts','D2_BRC_ligand_RMSD.mat'));

rmsd_all = calcrmsd(tempLigTraj(1,:),tempLigTraj);
rmsd_stats = zeros(length(foldersToStudy),2);
rmsf_stats = zeros(length(foldersToStudy),2);

% RMSD
figure
for thisSys = 1:length(foldersToStudy)
    numRuns = numRunsSys(thisSys);
    rmsd_here = rmsd_all(CSys(:,1)==thisSys);
    plot(rmsd_here);
    hold on
    rmsd_stats(thisSys,:) = [mean(rmsd_here) std(rmsd_here)/sqrt(numRuns)];
end
xlabel('Frame')
ylabel('RMSD [A]')
legend(thisSysLabel)
title('RSMD to BRC WT input model')
set(gca,'FontSize',20)

% RMSF
figure
for thisSys = 1:length(foldersToStudy)
    numRuns = numRunsSys(thisSys);
    rmsf_here = calcrmsf(tempLigTraj(CSys(:,1)==thisSys,:),[],numRuns);
    errorbar(1:(size(tempLigTraj,2)/3), rmsf_here(:, 1), rmsf_here(:, 2));
    hold on
    rmsf_stats(thisSys,:) = [mean(rmsf_here(:, 1)) mean(rmsf_here(:, 2))];
end
xlabel('Frame')
ylabel('RMSF [A]')
legend(thisSysLabel)
title('RSMF of BRC')
set(gca,'FontSize',20)

tab2 = table(thisSysLabel',rmsd_stats,rmsf_stats);

% Put it together
tab = [tab1; tab2];

% writetable(tab,'ligand_rmsd_stats.xlsx')