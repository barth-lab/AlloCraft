% Load the hyperDB first
dataTableMI = hdb.dbFetchData("DA","hyperEntryNdx",[1 3]);


dataTableKL = hdb.dbFetchData("DA","hyperEntryNdx",[1 7]);


%%
% BRC I4.46N
% wtName = "BRC-WT-Gi";
% sysName = "BRC-I4.46N-Gi";
% impRes = ["2.41" "2.46"  "3.29" "4.46" "7.39" "7.40" "7.43" "7.53"] ;

% DA I4.46N
wtName = "DA-WT-Gi";
sysName = "DA-I4.46N-Gi";
impRes = ["2.41" "2.45"  "3.42" "4.46" "6.32" "6.37" "6.44" "6.48" "6.58" "7.43" "7.53"] ;


% strcmp(database.residues{1}.Label,"1.52")
resIdsSim = zeros(length(impRes),1);

for i = 1:length(impRes)

resIdsSim(i) = database.residues{1}.output1((strcmp(database.residues{1}.Label,impRes(i))));

end



tempKL = dataTableKL(strcmp(dataTableKL.Label,sysName),"KLDiv wrt WT").Variables;
KL1Here = tempKL(resIdsSim);


tempMI = dataTableMI(strcmp(dataTableMI.Label,sysName),"MI per residue").Variables;
tempMIWT = dataTableMI(strcmp(dataTableMI.Label,wtName),"MI per residue").Variables;
MIHereTest = tempMI(resIdsSim);
MIHereWT = tempMIWT(resIdsSim);


tabtabtab = table(impRes',KL1Here',MIHereTest',MIHereWT',(MIHereTest - MIHereWT)','VariableNames',["Residue" "KL1" "MI" "WT MI" "delta(MI)"])

writetable(tabtabtab, ...
    fullfile('C:\Users\mahdi\Desktop\PhD EPFL\Manuscripts\D2 paper 2022\New submission 2023\figure elements\Supplementaries',sprintf("table_KL1_MI_%s.xlsx",sysName)));
