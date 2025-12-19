clear MIres;
load(fullfile(pathCalcdir,"workspace.mat"),"MIres","isGPCR");

 if ~exist('MIres','var') % Recalculate MIres fif it doesn't exist
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
    % Load excessMI:
    load(fullfile(pathCalcdir,"MI_.mat"));
    MIraw = I;
    load(fullfile(pathCalcdir,"MI_filtered.mat"));
    MI = excessMi;
   
    MIAnalysisDihLevel;
    removeAverageMI

    Nres = length(reslist);
    MIres = zeros(Nres);
     for i=1:Nres-1
        for j=i+1:Nres
            temp = MIfilterAvg(resind==(i),resind==(j));
            temp = temp(:);
            MIres(i,j) = sum(temp);
            MIres(j,i) = MIres(i,j);
        end
    end
    save(fullfile(pathCalcdir,"workspace.mat"),'MIres','-append');
 end
helices = load(fullfile(pathCalcdir,"helices.txt"));
BSres = load(fullfile(pathCalcdir,"BS_residues.txt"));
GPIres = load(fullfile(pathCalcdir,"GPI_residues.txt"));
BAIres = load(fullfile(pathCalcdir,"BAI_residues.txt"));
PDB = pdbread((fullfile(pathCalcdir,"protRenum.pdb")));
atomno = [PDB.Model.Atom(:).AtomSerNo];
resname = cellstr(char(PDB.Model.Atom(:).resName));
resno = [PDB.Model.Atom(:).resSeq];
reslist = unique(resno);
csts = getconstants;

% determine domain residues: we divide the GPCr in two ways, EC to IC
% domains, and transmembrane helices and loops (tmDomain)
if isGPCR
    ECres = [];
    ICres = [];
    ECLres = [];
    counter = 1;
    eclCounter = 1;
    iclCounter = 1;
    clear tmDomain
    if helices(1,1) > 3 % If Nterm is longer than 3 residues
        tmDomain(counter).name = "Nterm";
        tmDomain(counter).residues = reslist(1:helices(1,1)-1);
        counter = counter + 1;
    end

    for i=1:7
        tmDomain(counter).name = "TM" + num2str(i);
        tmDomain(counter).residues = reslist(helices(i,1):helices(i,2));
        counter = counter + 1;

        halfTMlength = int32(0.5*(helices(i,2)-helices(i,1)+1));
        Nhalf = reslist( helices(i,1):helices(i,1)+halfTMlength );
        Chalf = reslist( helices(i,1)+halfTMlength+1:helices(i,2) );
        if ~rem(i,2)
            ECres = [ECres Chalf];
            ICres = [ICres Nhalf];
            ECLres = [ECLres reslist(helices(i,2)+1:helices(i+1,1)-1)];

            % Even helices: ECLs:
            tmDomain(counter).name = "ECL" + num2str(eclCounter);
            tmDomain(counter).residues = reslist(helices(i,2)+1:helices(i+1,1)-1);
            eclCounter = eclCounter + 1;
            counter = counter + 1;
        else
            ECres = [ECres Nhalf];
            ICres = [ICres Chalf];

            % Odd helices: ICLs or Cterm
            if i ~= 7 
            tmDomain(counter).name = "ICL" + num2str(iclCounter);
            tmDomain(counter).residues = reslist(helices(i,2)+1:helices(i+1,1)-1);
            iclCounter = iclCounter + 1;
            counter = counter + 1;
            else
            tmDomain(counter).name = "Cterm";
            tmDomain(counter).residues = reslist(helices(i,2):end);
            end
        end
    end
    clear domain
    domain(1).residues = unique(ECres);
    domain(1).name = "EC";
    domain(2).residues = unique(ICres);
    domain(2).name = "IC";
    domain(3).residues = unique(BSres);
    domain(3).name = "BS";
    domain(4).residues = unique(GPIres);
    domain(4).name = "GPI";
    domain(5).residues = unique(ECLres);
    domain(5).name = "ECL";
    
    interDomainMIMean = zeros(length(domain),length(domain));
    interDomainMISum = zeros(length(domain),length(domain));

    for i=1:length(domain)
        for j=i:length(domain)
            resi = domain(i).residues;
            resj = domain(j).residues;
            interDomainMIMean(i,j) = mean(mean(MIres(resi,resj)));            
            interDomainMIMean(j,i) = interDomainMIMean(i,j);
            interDomainMISum(i,j) = sum(MIres(resi,resj),'all');
            interDomainMISum(j,i) = interDomainMISum(i,j);
        end
    end


    interTMDomainMIMean = zeros(length(tmDomain),length(tmDomain));
    interTMDomainMISum = zeros(length(tmDomain),length(tmDomain));
    for i=1:length(tmDomain)
        for j=i:length(tmDomain)
            resi = tmDomain(i).residues;
            resj = tmDomain(j).residues;
            interTMDomainMIMean(i,j) = mean(mean(MIres(resi,resj)));            
            interTMDomainMIMean(j,i) = interTMDomainMIMean(i,j);
            interTMDomainMISum(i,j) = sum(MIres(resi,resj),'all');
            interTMDomainMISum(j,i) = interTMDomainMISum(i,j);
        end
    end
end

if doPlot 
% Plot heat map with MI:
% figure('Position',[163 270 1400 1000]); 
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
heatmap([domain.name],[domain.name],interDomainMISum)
title('Interdomain MI')
axs = struct(gca); %ignore warning that this should be avoided
cb = axs.Colorbar;
cb.Label.String = '\Sigma_{i,j}MI(i,j)';
cb.Label.FontSize = 16; 
subplot(2,2,2)
heatmap([domain.name],[domain.name],interDomainMIMean)
title('Mean interdomain MI')
axs = struct(gca); %ignore warning that this should be avoided
cb = axs.Colorbar;
cb.Label.String = 'mean[MI(i,j)]';
cb.Label.FontSize = 16; 
subplot(2,2,3)
heatmap([tmDomain.name],[tmDomain.name],interTMDomainMISum)
title('Interdomain MI')
axs = struct(gca); %ignore warning that this should be avoided
cb = axs.Colorbar;
cb.Label.String = '\Sigma_{i,j}MI(i,j)';
cb.Label.FontSize = 16; 
subplot(2,2,4)
heatmap([tmDomain.name],[tmDomain.name],interTMDomainMIMean)
title('Mean interdomain MI')
axs = struct(gca); %ignore warning that this should be avoided
cb = axs.Colorbar;
cb.Label.String = 'mean[MI(i,j)]';
cb.Label.FontSize = 16; 
colormap summer
sgtitle(name,'FontSize', 20)

%% Save the figure
savefig(fullfile(pathCalcdir,"domainMI"))
print2pdf(fullfile(pathCalcdir,"domainMI"))
end