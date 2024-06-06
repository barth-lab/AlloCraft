library(bio3d)

# Prepare input:
# Define the WT file identifier (modify this as per your WT file naming convention)
wt_identifier <- "WT"
nmodes = 20
ncore = 1 # Change this as needed

# Takes two file lists and two sets of hubs as input and calculates dd(sum(sum(DCCM))) over allosteric hubs
file_list_state1 = read.table('file_list_ina.txt')
file_list_state2 = read.table('file_list_act.txt')
# state1_name = "inactive"
# state2_name = "active"

#D2 Mahdi hubs, DA
# hubs_state1 <- c(249,33,252,107,183,95, 244,197,193,100,31,253,256,177,181,159,109,34,115,184)
# hubs_state2 <- c(261,38,264,112,188,100,256,209,205,105,36,265,268,182,186,164,114,39,120,189)
hubs_state1 <-scan("hubs_state1.txt", what = numeric(), sep = " ")
hubs_state2 <-scan("hubs_state2.txt", what = numeric(), sep = " ")

# Compute functions
compute_ssDCCM <- function(pdb, hubs, nmodes = 20, ncore = 1) {
    mode_pdb_data <- nma(pdb)
    dccm_pdb_data <- dccm(mode_pdb_data, nmodes=nmodes, ncore=ncore)
    
    ssDCCM = 0.0

    for(x in hubs) { 

        for(y in hubs) { 

            if (y>x) { 

                ssDCCM = ssDCCM + dccm_pdb_data[x,y]

                # Potentially use absolute value so +ve and -ve corrs don't cancel out
                # ssDCCM = ssDCCM + abs(dccm_pdb_data[x,y])

            } 
        }
    }
    return(ssDCCM)
}

compute_dssDCCM_list <- function(file_list, hubs, wt_identifier, nmodes = 20, ncore = 1){
    ssDCCM_list = c()
    for (file in file_list[['V1']]) {

        print(file)

        # Is the file a WT?
        # Check if the file corresponds to the WT
        isWT <- grepl(wt_identifier, file, ignore.case = TRUE)
        
        pdb_path = paste('./' , trimws(file) , sep='')
        
        pdb_data = read.pdb(pdb_path)
        
        ssDCCM <- compute_ssDCCM(pdb_data, hubs, nmodes, ncore)

        ssDCCM_list = append(ssDCCM_list, ssDCCM)
        if (isWT) {
            ssDCCMWT = ssDCCM
        }
    }
    return(list(ssDCCM_list = ssDCCM_list, ssDCCMWT = ssDCCMWT))
}

# "Main" part of the script
listHere <- compute_dssDCCM_list(file_list_state1, hubs_state1, wt_identifier, nmodes, ncore)
ssDCCM_list_state1 <- listHere$ssDCCM_list
ssDCCMWT_state1 <- listHere$ssDCCMWT

listHere <- compute_dssDCCM_list(file_list_state2, hubs_state2, wt_identifier, nmodes, ncore)
ssDCCM_list_state2 <- listHere$ssDCCM_list
ssDCCMWT_state2 <- listHere$ssDCCMWT
dssDCCM_list_state1 = ssDCCM_list_state1 - ssDCCMWT_state1
dssDCCM_list_state2 = ssDCCM_list_state2 - ssDCCMWT_state2

write.table(dssDCCM_list_state2 - dssDCCM_list_state1, 'ddssdccm_raw.txt', sep = '\n', row.names = FALSE, col.names = FALSE)

# Combine file names and numbers into a data frame
data <- data.frame(file_list_state1 = file_list_state1, file_list_state2 = file_list_state2, ssDCCM_list_state1 = ssDCCM_list_state1, 
 ssDCCM_list_state2 = ssDCCM_list_state2, dssDCCM_list_state1 = dssDCCM_list_state1, dssDCCM_list_state2 = dssDCCM_list_state2,
 ddssDCCM = dssDCCM_list_state2 - dssDCCM_list_state1)

output_file <- "output_ddssdccm_all.txt"
write.table(data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE)