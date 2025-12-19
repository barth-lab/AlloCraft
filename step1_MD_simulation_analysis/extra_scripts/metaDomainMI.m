% Meta domainMI
interTMDomainMISumCell = cell(length(settings.foldersToStudy),1);
interTMDomainMIMeanCell = cell(length(settings.foldersToStudy),1);
interDomainMISumCell = cell(length(settings.foldersToStudy),1);
interDomainMIMeanCell = cell(length(settings.foldersToStudy),1);
% [domain.name]
% [tmDomain.name]
doPlot = 0; % Plot heatmaps?
% xlsname = l;
for thisSys=1:length(settings.foldersToStudy)
    pathCalcdir = fullfile(settings.metadir,settings.foldersToStudy{thisSys},md2pathName,alloPathCalcName);
    name = settings.thisSysLabel{thisSys};
    % Run the script
    domainMI;
    % Save shit:
    interTMDomainMISumCell{thisSys} = interTMDomainMISum;
    interTMDomainMIMeanCell{thisSys} = interTMDomainMIMean;
    interDomainMISumCell{thisSys} = interDomainMISum;
    interDomainMIMeanCell{thisSys} = interDomainMIMean;
end
%% Plot shit?
% Which plot type to use?

figure;

for i = 1:3
% Scatter3
[X,Y] = meshgrid(1:1:length(interTMDomainMISumCell{i}));
X = triu(X);
Y = triu(Y);
data = triu(interTMDomainMISumCell{i});
% figure; scatter3([X(:)],[Y(:)],[interTMDomainMISumCell{1}(:)])

scatter3([X(:)],[Y(:)],[data(:)],20,'filled')
hold on
end
%% Mesh3:
[X,Y] = meshgrid(1:1:length(interTMDomainMISumCell{1}));
figure
mesh(X,Y,interTMDomainMISumCell{1},interTMDomainMISumCell{1},'EdgeColor','none','FaceColor','flat')


%% Try another way:
interTMDomainMISum3D = zeros(length(interTMDomainMISumCell{1}),length(interTMDomainMISumCell{1}),length(interTMDomainMISumCell));
interTMDomainMIMean3D = zeros(length(interTMDomainMISumCell{1}),length(interTMDomainMISumCell{1}),length(interTMDomainMISumCell));
for i =1:length(interTMDomainMISumCell)
    interTMDomainMISum3D(:,:,i) = interTMDomainMISumCell{i};
    interTMDomainMIMean3D(:,:,i) = interTMDomainMIMeanCell{i};
end

interDomainMISum3D = zeros(length(interDomainMISumCell{1}),length(interDomainMISumCell{1}),length(interDomainMISumCell));
interDomainMIMean3D = zeros(length(interDomainMISumCell{1}),length(interDomainMISumCell{1}),length(interDomainMISumCell));
for i =1:length(interDomainMISumCell)
    interDomainMISum3D(:,:,i) = interDomainMISumCell{i};
    interDomainMIMean3D(:,:,i) = interDomainMIMeanCell{i};
end

% Which variables to care about? Assumption: go with highest variance
% variables

stdInterDomainMISum = std(interDomainMISum3D,[],3);
stdInterTMDomainMISum = std(interTMDomainMISum3D,[],3);

figure; 
tiledlayout('flow')
nexttile
heatmap([domain.name],[domain.name],std(interDomainMISum3D,[],3))
title('\sigma(Interdomain MI)')
nexttile
heatmap([domain.name],[domain.name],std(interDomainMIMean3D,[],3))
title('\sigma(Mean interdomain MI)')
nexttile
heatmap([tmDomain.name],[tmDomain.name],std(interTMDomainMISum3D,[],3))
title('\sigma(Interdomain MI)')
nexttile
heatmap([tmDomain.name],[tmDomain.name],std(interTMDomainMIMean3D,[],3))
title('\sigma(Mean interdomain MI)')

sgtitle('Standard deviation over all studied systems','FontSize', 20)
%% x-axis are the systems, different plots are domains

pairs2study = [8 11 11 13;8 11 13 13];
i = 13;
j = 13;
shapes = {'o','s','+','x','^','d'};
counter = 1;
figure;
for pairHere = pairs2study
    i = pairHere(1);
    j = pairHere(2);
    nameHere = "MI(" + tmDomain(i).name + "," + tmDomain(j).name + ")";
    scatter(1:size(interTMDomainMISum3D,3),reshape(interTMDomainMISum3D(i,j,:),[],1,1),20, ...
      shapes{counter},'LineWidth',1.5,'DisplayName',nameHere)
    hold on
    counter = counter + 1;
end
% title("MI(" + tmDomain(i).name + "," + tmDomain(j).name + ")")
ylabel('MI')
xticks(1:1:length(settings.thisSysLabel))
xticklabels(settings.thisSysLabel)
legend('-DynamicLegend');
legend boxoff
set(gca,'FontSize',16)

%% x-axis are domains, different plots are systems

% Set up the parameters for the plot
stdCutoff = 15;
inactiveNdx = find(strcmp('RIS-inactive',settings.thisSysLabel));

MIHere = interDomainMISum3D;
stdHere = std(MIHere,[],3);
domainHere = domain;

xtickNames = [];
figure; 
counter = 1;
for i = 1:length(domainHere)
    for j = i:length(domainHere)
        if stdHere(i,j) < stdCutoff % If std is too low, don't bother plotting
            continue
        end
        s3 = scatter(counter*ones(size(MIHere,3),1),reshape(MIHere(i,j,:),[],1,1),20, ...
      'filled');
        hold on
        row = dataTipTextRow('System',settings.thisSysLabel);
        s3.DataTipTemplate.DataTipRows(end+1) = row;
        scatter(counter,MIHere(i,j,inactiveNdx),'rx','LineWidth',1.5)

        xtickNames{counter} = "MI(" + domainHere(i).name + "," + domainHere(j).name + ")";
        counter = counter + 1;
    end
end
ylabel('MI')
xticks(1:1:length(xtickNames))
xticklabels(xtickNames)
set(gca,'FontSize',16)
title('Inter-domain MI, summed','FontSize',20)

%% x-axis are variants, different plots are domains

% Set up the parameters for the plot
stdCutoff = 15;
inactiveNdx = find(strcmp('RIS-inactive',settings.thisSysLabel));

MIHere = interTMDomainMISum3D;
stdHere = std(MIHere,[],3);
domainHere = tmDomain;

xtickNames = [];
figure; 
counter = 1;
for i = 1:length(domainHere)
    for j = i:length(domainHere)
        if stdHere(i,j) < stdCutoff % If std is too low, don't bother plotting
            continue
        end
        s3 = scatter(1:size(MIHere,3),reshape(MIHere(i,j,:),[],1,1),20, ...
      'filled');
        hold on
        row = dataTipTextRow('Domain',repmat(domainHere(i).name + "-" + domainHere(j).name,size(MIHere,3),1));
        s3.DataTipTemplate.DataTipRows(end+1) = row;
%         scatter(counter,MIHere(i,j,inactiveNdx),'rx','LineWidth',1.5)

%         xtickNames{counter} = "MI(" + domainHere(i).name + "," + domainHere(j).name + ")";
%         counter = counter + 1;
    end
end
ylabel('MI')
% xticks(1:1:length(xtickNames))
% xticklabels(xtickNames)
xticks(1:1:length(settings.thisSysLabel))
xticklabels(settings.thisSysLabel)
set(gca,'FontSize',16)
title('Inter-domain MI, summed','FontSize',20)
