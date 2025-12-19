% Example of using GLM
% A set of car weights
weight = [2100 2300 2500 2700 2900 3100 3300 3500 3700 3900 4100 4300]';
% The number of cars tested at each weight
tested = [48 42 31 34 31 21 23 23 21 16 17 21]';
% The number of cars failing the test at each weight
failed = [1 2 0 3 8 8 14 17 19 15 17 21]';
% The proportion of cars failing for each weight
proportion = failed ./ tested;

figure
plot(weight,proportion,'s')
xlabel('Weight'); ylabel('Proportion');

%
[logitCoef,dev] = glmfit(weight,[failed tested],'binomial','logit');
logitFit = glmval(logitCoef,weight,'logit');
figure
plot(weight,proportion,'bs', weight,logitFit,'r-');
xlabel('Weight'); ylabel('Proportion');

[logitCoef2,dev2] = glmfit([weight weight.^2],[failed tested],'binomial','logit');
pval = 1 - chi2cdf(dev-dev2,1)

%% Prep my data:
% Experimental data:
tabExp = readtable('D:\Telework_library\dopamine_phase_3\a-analysis\Experimental data\experimental_data_labeled.xlsx');

%% All my other shit:

MIPerRes; % Loaded from simple_pca.m
thisSysLabel;
ndxExpSim = zeros(length(tabExp.Label),1);
% look for the same labels in experimental data and variable 
for i = 1:length(tabExp.Label)

% strcmp(string(thisSysLabel),string(tabExp.Label(i)))
    ndxHere = find(strcmp(string(thisSysLabel),string(tabExp.Label(i))));
    if ndxHere & (~isnan(tabExp.LogEC50(i)) | ~isnan(tabExp.Efficacy(i)))
        ndxExpSim(i,1) = i;
        ndxExpSim(i,2) = ndxHere;
    end
end

ndxCommon=ndxExpSim(any(ndxExpSim,2),any(ndxExpSim,1)); % Remove rows with zeros
%%
MIPerResFiltered=MIPerRes(any(MIPerRes,2),any(MIPerRes,1));
X = sum(MIPerResFiltered(:,ndxCommon(:,2)))';
Y = [tabExp.Efficacy(ndxCommon(:,1)) -tabExp.LogEC50(ndxCommon(:,1))];

[Coef,dev] = glmfit(X,Y);

logitFit = glmval(Coef,MIPerRes(:,ndxCommon(:,2)));
figure
plot(weight,proportion,'bs', weight,logitFit,'r-');
xlabel('Weight'); ylabel('Proportion');

%% Try simpler regression
X = [ones(length(ndxCommon),1) sum(MIPerResFiltered(:,ndxCommon(:,2)))'];
% y = log(Y(:,1)).*Y(:,2);

[b,bint,r,rint,stats] = regress(Y(:,2),X);

yFit = b(1) + b(2)*linspace(min(X(:,2)),max(X(:,2)),100);


%%
ligNames = ["DA" "BRC"];
effectorNames = ["Gi" "Barr2"];
markers ={'o','d','s','^'};
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];

xRange = range(X(:,2));
yRange = range( Y(:,1));

thisSysLabel     
tabExp.Label
ndxCommon

figure
tiledlayout('flow')
for ligHere = 1:length(ligNames)
        ligNdx = contains(thisSysLabel(ndxCommon(:,2)),ligNames(ligHere));
        for effectorHere = 1:length(effectorNames)
            % Plot the stuff
            nexttile
            effectorNdx = contains(thisSysLabel(ndxCommon(:,2)),effectorNames(effectorHere));

            Xhere = sum(MIPerResFiltered(:,ndxCommon(ligNdx&effectorNdx,2)))';
            Yhere = tabExp.Efficacy(ndxCommon(ligNdx&effectorNdx,1)) ;

            s = scatter(Xhere, Yhere,25, ...
               colors(effectorHere,:),'Marker',markers{ligHere},'LineWidth',1.5, ...
               'DisplayName', ligNames(ligHere)+"-"+effectorNames(effectorHere));
            hold on
            row = dataTipTextRow('Label',thisSysLabel(ndxCommon(ligNdx&effectorNdx,2)));
            s.DataTipTemplate.DataTipRows(end+1) = row;
            text(Xhere(:,1)+xRange/50,Yhere(:,1)+yRange/50, thisSysLabel(ndxCommon(ligNdx&effectorNdx,2)));   
            
            % Regression:
            [b,bint,r,rint,stats] = regress(Yhere,[ones(size(Xhere)) Xhere]);
            yFit = b(1) + b(2)*linspace(min(Xhere),max(Xhere),100);
            
            plot(linspace(min(Xhere),max(Xhere),100),yFit,'--r')
            legend('Data',['Linear fit, R^2=' num2str(stats(1))])
            legend boxoff
            % Labels and formatting
            xlabel('\Sigma(MI)')
            ylabel('Efficacy (normalized to WT)')
            title( ligNames(ligHere)+"-"+effectorNames(effectorHere))
            set(gca,'FontSize',12)
        end
end





% text(X(:,2)+xRange/50,Y(:,1)+yRange/50, thisSysLabel(ndxCommon(:,2)));   
% set(gca,'FontSize',16)