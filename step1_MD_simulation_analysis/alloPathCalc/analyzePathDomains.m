% reslist = 1:Nres;
BSres = load(fullfile(pathCalcdir,"BS_residues.txt"));
GPIres = load(fullfile(pathCalcdir,"GPI_residues.txt"));
BAIres = load(fullfile(pathCalcdir,"BAI_residues.txt"));

% find all channels that include any BS residue
for i=1:length(channelstruc)
    % find all pathways belonging to cluster i
    clslist = [pathstruc(I(1:Npath)).cls]==i;
    clsres = unique([pathstruc(I(clslist)).path]);
    BScount = sum(ismember(BSres,clsres));
    if(BScount>=0.5*length(BSres))
        i
    end
end

for i=1:Npath
    res_terminus(i,1) = pathstruc(I(i)).path(1);
    res_terminus(i,2) = pathstruc(I(i)).path(end);
end
% determine domain residues
if isGPCR
    ECres = [];
    ICres = [];
    ECLres = [];
    for i=1:7
        halfTMlength = int32(0.5*(helices(i,2)-helices(i,1)+1));
        Nhalf = reslist( helices(i,1):helices(i,1)+halfTMlength );
        Chalf = reslist( helices(i,1)+halfTMlength+1:helices(i,2) );
        if ~rem(i,2)
            ECres = [ECres Chalf];
            ICres = [ICres Nhalf];
            ECLres = [ECLres reslist(helices(i,2)+1:helices(i+1,1)-1)];
        else
            ECres = [ECres Nhalf];
            ICres = [ICres Chalf];
        end
    end
    clear domain
    domain(1).residues = unique(ECres);
    domain(1).name = 'EC';
    domain(2).residues = unique(ICres);
    domain(2).name = 'IC';
    domain(3).residues = unique(BSres);
    domain(3).name = 'BS';
    domain(4).residues = unique(GPIres);
    domain(4).name = 'GPI';
    domain(5).residues = unique(ECLres);
    domain(5).name = 'ECL';
    
    interDomainMI = zeros(length(domain),length(domain));
    domainsfile = fopen(fullfile(pathCalcdir,"domains_MI.txt"),'w');
    outstring = [];
    for i=1:length(domain)
        for j=i:length(domain)
            resi = domain(i).residues;
            resj = domain(j).residues;
            interDomainMI(i,j) = mean(mean(MIres(resi,resj)));
            interDomainMI(j,i) = interDomainMI(i,j);
            writestring = sprintf('%s-%s\t%f\n',...
                domain(i).name,domain(j).name,interDomainMI(i,j));
            outstring = [outstring writestring];
            fprintf(domainsfile,'%s',writestring);
        end
    end
    fclose(domainsfile); 
    outstringMI = outstring;
    
    % count all pathways in the top X belonging to different domains
    interDomainPathcount = zeros(length(domain),length(domain));
    domainsfile = fopen(fullfile(pathCalcdir,"domains_pathcount.txt"),'w');
    outstring = [];
    for i=1:length(domain)
        for j=i:length(domain)
            resi = domain(i).residues;
            resj = domain(j).residues;
            interDomainPathcount(i,j)...
                = sum( (ismember(res_terminus(:,1),resi) & ismember(res_terminus(:,2),resj))...
                | (ismember(res_terminus(:,2),resi) & ismember(res_terminus(:,1),resj)) );
            interDomainPathcount(j,i) = interDomainPathcount(i,j);
            writestring = sprintf('%s-%s\t%f\n',...
                domain(i).name,domain(j).name,interDomainPathcount(i,j));
            outstring = [outstring writestring];
            fprintf(domainsfile,'%s',writestring);
        end
    end
    fclose(domainsfile); 
    outstringPathcount = outstring;
end

% Collect all pathways leading to GPI and BAI
count = 1;
ResOriginFreq = zeros(Nres,2);
for i=1:Npath
    if (sum(ismember(res_terminus(i,:),GPIres))...
            +sum(ismember(res_terminus(i,:),BAIres)))
        GPI_BAI_pathways(count) = i;
        if ismember(res_terminus(i,1),GPIres)
            ResOriginFreq(res_terminus(i,2),1) = ResOriginFreq(res_terminus(i,2),1)+1;
        end
        if ismember(res_terminus(i,2),GPIres)
            ResOriginFreq(res_terminus(i,1),1) = ResOriginFreq(res_terminus(i,1),1)+1;
        end
        if ismember(res_terminus(i,1),BAIres)
            ResOriginFreq(res_terminus(i,2),2) = ResOriginFreq(res_terminus(i,2),2)+1;
        end
        if ismember(res_terminus(i,2),BAIres)
            ResOriginFreq(res_terminus(i,1),2) = ResOriginFreq(res_terminus(i,1),2)+1;
        end
        count = count+1;
    end
end
[~,index] = sortrows(ResOriginFreq(:,1)+ResOriginFreq(:,2),-1);
writematrix([reslist(index)',ResOriginFreq(index,:)],fullfile(pathCalcdir,"res_origin_freq.txt"),'Delimiter','tab');

% Identify intermediate residues and rank them
ResIntermFreq = zeros(Nres,2);
for i=1:Npath
    if (sum(ismember(res_terminus(i,:),GPIres))...
            +sum(ismember(res_terminus(i,:),BAIres)))
        if sum(ismember(res_terminus(i,:),GPIres))
            ResIntermFreq(pathstruc(I(i)).path(2:end-1),1)...
                = ResIntermFreq(pathstruc(I(i)).path(2:end-1),1)+1;
        end
        if sum(ismember(res_terminus(i,:),BAIres))
            ResIntermFreq(pathstruc(I(i)).path(2:end-1),2)...
                = ResIntermFreq(pathstruc(I(i)).path(2:end-1),2)+1;
        end
    end
end
[~,index] = sortrows(ResIntermFreq(:,1)+ResIntermFreq(:,2),-1);
writematrix([reslist(index)',ResIntermFreq(index,:)],fullfile(pathCalcdir,"res_interm_freq.txt"),'Delimiter','tab');


if exist("mainEntry",'var') % Add reference and BW numbering for convenience
    % Visualize density of pathways from binding site residues and 
    density = plotAlloDensityFun(PDB,pathstruc,BSres,'cut',[0.1 0.2 0.3 0.5],'MIFractionCutoff',MIFractionCutoff, ...
        'title', 'BS residues','mainEntry',mainEntry,'includeLigPathwayCalc',includeLigPathwayCalc);
    savefig(fullfile(pathCalcdir,"alloDensity_bs"))
    % Pathways that go all the way from BS to GPI are scarce, so the cutoff is
    % lower than default
    plotAlloDensityFun(PDB,pathstruc,BSres,'MIFractionCutoff',MIFractionCutoff,'optRes',GPIres,'cut',[0.1 0.2 0.3 0.5] ...
        ,'title', 'BS2GP','mainEntry',mainEntry,'includeLigPathwayCalc',includeLigPathwayCalc);
    savefig(fullfile(pathCalcdir,"alloDensity_bs2gp"))
else
    % Visualize density of pathways from binding site residues and 
    density = plotAlloDensityFun(PDB,pathstruc,BSres,'cut',[0.1 0.2 0.3 0.5],'MIFractionCutoff',MIFractionCutoff,'title', 'BS residues');
    savefig(fullfile(pathCalcdir,"alloDensity_bs"))
    % Pathways that go all the way from BS to GPI are scarce, so the cutoff is
    % lower than default
    plotAlloDensityFun(PDB,pathstruc,BSres,'MIFractionCutoff',MIFractionCutoff,'optRes',GPIres,'cut',[0.1 0.2 0.3 0.5],'title', 'BS2GP');
    savefig(fullfile(pathCalcdir,"alloDensity_bs2gp"))
end
 save(fullfile(pathCalcdir,"workspace.mat"),'density','-append')
% If protein is inversed along z, use this command to make it upright again
% set(gca,'zdir','reverse')