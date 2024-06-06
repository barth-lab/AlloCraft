library(bio3d)

file_list = read.table('file_list.txt')

# Sample hubs from 5ht1b simultions
# hubs <- c(126, 268, 123, 122, 278, 275, 273, 127, 47, 282, 54, 121, 266, 270, 49, 131, 279, 43, 137, 226)
hubs <-scan("hubs.txt", what = numeric(), sep = " ")

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


# Define the WT file identifier (modify this as per your WT file naming convention)
wt_identifier <- "WT"
nmodes = 20
ncore = 1

listHere <- compute_dssDCCM_list(file_list, hubs, wt_identifier, nmodes, ncore)
ssDCCM_list <- listHere$ssDCCM_list
ssDCCMWT <- listHere$ssDCCMWT

# # write.table(ssDCCM_list, '5ht1b_dmut_20hubs_bio3d_dccm.txt', sep='\n')
write.table(ssDCCM_list, 'ssdccm_raw.txt', sep = '\n', row.names = FALSE, col.names = FALSE)

# # write.table(ssDCCM_list - ssDCCMWT, '5ht1b_dmut_20hubs_bio3d_dccm.txt', sep='\n')
write.table(ssDCCM_list - ssDCCMWT, 'dssdccm_raw.txt', sep = '\n', row.names = FALSE, col.names = FALSE)

# Combine file names and numbers into a data frame
data <- data.frame(File_Name = file_list, ssDCCM_list = ssDCCM_list, dssDCCM_list = ssDCCM_list - ssDCCMWT)
output_file <- "output_dssdccm_all.txt"
write.table(data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE)