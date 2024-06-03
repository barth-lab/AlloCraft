%==================================================================================================================
%==================================================================================================================
% Those are defined in the md2path script
% all_BS = [9 83 172 175 176 183 253 259 260 261 264 268]; %to be filled in before running
% all_EC = [1:24 78:99 151:188 246:268];  %to be filled in before running

%==================================================================================================================
%==================================================================================================================

%Initialization begins

EC_BS_BAI_path = [];
EC_BS_GPI_path = [];

ec_bs_gpi_hubs = zeros(Nres,1);
ec_bs_bai_hubs = zeros(Nres,1);

%Initialization ends

%topcutoff20 = round(sum([pathstruc(:).MI]>0)*0.20);
%topcutoff05 = round(sum([pathstruc(:).MI]>0)*0.05);
%topcutoff15 = round(sum([pathstruc(:).MI]>0)*0.15);
% topcutoff10 = round(sum([pathstruc(:).MI]>0)*0.100);



for i=1:Npath
    path = pathstruc(I(i)).path;
    if( (sum(ismember(path([1,end]),all_EC))>0) &&...
        (sum(ismember(path([1,end]),BAIres))>0) &&...
        (sum(ismember(path,all_BS))>0) )
        EC_BS_BAI_path = [I(i) EC_BS_BAI_path]; %Pathways passing from EC through BS to BAI
		for j=1:length(path)
			ec_bs_bai_hubs(pathstruc(I(i)).path(j))=ec_bs_bai_hubs(pathstruc(I(i)).path(j))+1;
		end
    end
    if( (sum(ismember(path([1,end]),all_EC))>0) &&...
        (sum(ismember(path([1,end]),GPIres))>0) &&...
        (sum(ismember(path,all_BS))>0) )
        EC_BS_GPI_path = [I(i) EC_BS_GPI_path]; %Pathways passing from EC through BS to GPI
		for j=1:length(path)
			ec_bs_gpi_hubs(pathstruc(I(i)).path(j))=ec_bs_gpi_hubs(pathstruc(I(i)).path(j))+1;
		end
    end
end

% dlmwrite('gpi_bai_hubs.txt',[ec_bs_gpi_hubs ec_bs_bai_hubs],' ') %residue-wise hubscore for communication from ec region, through lig binding site to g-ptn interface and b-arr interface
% dlmwrite('lengths.txt',[length(EC_BS_BAI_path) length(EC_BS_GPI_path)],' ') %numerator and denominator to calculate pathway ratio for arr bias.
writematrix([ec_bs_gpi_hubs ec_bs_bai_hubs], fullfile(pathCalcdir,"gpi_bai_hubs.txt"), ...
    'Delimiter'," ")
writematrix([length(EC_BS_BAI_path) length(EC_BS_GPI_path)],fullfile(pathCalcdir,"lengths.txt"), ...
    'Delimiter'," ")

%==================================================================================================================
