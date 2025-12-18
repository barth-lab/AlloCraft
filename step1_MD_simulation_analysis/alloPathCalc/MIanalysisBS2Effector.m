clear alloMI BS_GPI GPI
nMI = length(MIres);
MIsd = std(MIres(:));
count = 1;
for i=1:nMI-1
    for j=i+1:nMI
        if( (MIres(i,j)>0) && (dismat(i,j)>disCutOff) )
            %MIzscore = (MIres(i,j)-MIavg)/MIsd;
            ptile = sum(MIres(:)>MIres(i,j) & dismat(:)>disCutOff)/sum(MIres(:)>0 & dismat(:)>disCutOff)*100;
            alloMI(count,:) = [i,j,ptile MIres(i,j)];
            count  = count+1;
        end
    end
end

% calculate MI between BI and GPI
count = 1;
count1 = 1;
for i=1:length(alloMI)
    if( (ismember(alloMI(i,1),BSres) && ismember(alloMI(i,2),GPIres)) ||...
            (ismember(alloMI(i,2),BSres) && ismember(alloMI(i,1),GPIres)) )
        BS_GPI(count,:) = alloMI(i,:);
        count = count+1;
    end
    if( (~ismember(alloMI(i,1),GPIres) && ismember(alloMI(i,2),GPIres)) ||...
            (~ismember(alloMI(i,2),GPIres) && ismember(alloMI(i,1),GPIres)) )
        GPI(count1,:) = alloMI(i,:);
        count1 = count1+1;
    end
end

% write pdb with per residue allosteric MI as temp factor
% Allosteric MI is considered to be the sum of MI with residues beyond
% disCutOff
PDBout = PDB;
MIperResAllo = zeros(Nres,1);
for i=1:Nres
    MIperResAllo(i) = sum(MIres(i,:).*(dismat(i,:)>disCutOff));
    tf = MIperResAllo(i) / sum((dismat(i,:)>disCutOff)); % Normalize by number of "far" residues
    if(isnan(tf))
        tf = 0;
    end
    [PDBout.Model.Atom(resno==i).tempFactor] = deal(tf); 

end
pdbwrite(fullfile(pathCalcdir,"protein_MI.pdb"),PDBout)
netGPIMI = sum(GPI(:,4));

% Run the analyzeHubsAll script:
analyzeHubsAll;

if isGPCR
    % calculate MI and number of pathways between BS and TM5/TM6 vs TM7
    TM5ICres = helices(5,1):helices(5,2); 
    TM6ICres = helices(6,1):helices(6,2); 
    TM7ICres = helices(7,1):helices(7,2); 
    temp = MIres(BSres,TM5ICres);
    MI_BS_TM5IC = mean(temp(:));
    temp = MIres(BSres,TM6ICres);
    MI_BS_TM6IC = mean(temp(:));
    temp = MIres(BSres,TM7ICres);
    MI_BS_TM7IC = mean(temp(:));
    
    % Number of pathways from BS to IC ends
    count_pth_TM5IC = 0;
    count_pth_TM6IC = 0;
    count_pth_TM7IC = 0;
    count_pth_BS_TM5IC = 0;
    count_pth_BS_TM6IC = 0;
    count_pth_BS_TM7IC = 0;
    for i=1:length(pathstruc)
        pth = pathstruc(i).path;
        if(~isempty(pth))
            if(sum(ismember(TM5ICres,pth))>0)
                count_pth_TM5IC = count_pth_TM5IC+1;
            end
            if(sum(ismember(TM6ICres,pth))>0)
                count_pth_TM6IC = count_pth_TM6IC+1;
            end
            if(sum(ismember(TM7ICres,pth))>0)
                count_pth_TM7IC = count_pth_TM7IC+1;
            end
            %         if(ismember(pth(1),BSres) || ismember(pth(end),BSres))
            %             if(sum(ismember(TM5ICres,pth))>0)
            %                 count_pth_BS_TM5IC = count_pth_BS_TM5IC+1;
            %             end
            %             if(sum(ismember(TM6ICres,pth))>0)
            %                 count_pth_BS_TM6IC = count_pth_BS_TM6IC+1;
            %             end
            %             if(sum(ismember(TM7ICres,pth))>0)
            %                 count_pth_BS_TM7IC = count_pth_BS_TM7IC+1;
            %             end
            %         end
            if(ismember(pth(1),BSres) || ismember(pth(end),BSres))
                if(ismember(pth(1),TM5ICres) || ismember(pth(end),TM5ICres))
                    count_pth_BS_TM5IC = count_pth_BS_TM5IC+1;
                end
                if(ismember(pth(1),TM6ICres) || ismember(pth(end),TM6ICres))
                    count_pth_BS_TM6IC = count_pth_BS_TM6IC+1;
                end
                if(ismember(pth(1),TM7ICres) || ismember(pth(end),TM7ICres))
                    count_pth_BS_TM7IC = count_pth_BS_TM7IC+1;
                end
            end
        end
    end
    % dlmwrite('Gs_barr_pathways.dat',[MI_BS_TM5IC MI_BS_TM6IC MI_BS_TM7IC...
    %     count_pth_BS_TM5IC count_pth_BS_TM6IC count_pth_BS_TM7IC...
    %     count_pth_TM5IC count_pth_TM6IC count_pth_TM7IC],' ');
    writematrix([MI_BS_TM5IC MI_BS_TM6IC MI_BS_TM7IC...
        count_pth_BS_TM5IC count_pth_BS_TM6IC count_pth_BS_TM7IC...
        count_pth_TM5IC count_pth_TM6IC count_pth_TM7IC], fullfile(pathCalcdir,"Gs_barr_pathways.dat"), ...
        'Delimiter'," ");
end

count = 1;
pathstruc1 = [];
for i=1:length(BS_GPI)
    res1 = BS_GPI(i,1);
    res2 = BS_GPI(i,2);
    
    % find all MI pairs involving res1
    temp = alloMI(logical(sum(ismember(alloMI(:,1:2),res1),2)),1:2);
    temp = unique(temp(:));
    temp = temp(~ismember(temp,BS_GPI(i,1:2)));
    
    % find if any of these residues also show high MI with res2
    temp = logical(sum(ismember(alloMI(:,1:2),temp),2)) & ...
        logical(sum(ismember(alloMI(:,1:2),res2),2));
    temp = alloMI(temp,1:2);
    temp = unique(temp(:));
    temp = temp(temp~=res2);
    path = [[res1; temp; res2]];% MIres([res1; temp; res2],[res1 res2])];
    if(isempty(temp))
        path = [res1 res2];
        pathstruc1(count).path = path;
        pathstruc1(count).Npath = length(path);
        pathstruc1(count).MI = MIres(res1,res2);
        pathstruc1(count).meanpathMI = 0;
        pathstruc1(count).cls = i;
    else
        for j=1:length(temp)
            path = [res1 temp(j) res2];
            pathstruc1(count).path = path;
            pathstruc1(count).Npath = length(path);
            pathstruc1(count).MI = MIres(res1,res2);
            pathstruc1(count).meanpathMI = 0;
            pathstruc1(count).cls = i;
        end
    end
    count = count+1;
end
writepathwayClusters(PDB,pathstruc1,1:length(pathstruc1),length(pathstruc1),fullfile(char(pathCalcdir),'BS_GPI'))
