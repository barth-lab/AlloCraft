# Define the special identifier
wt_identifier <- "WT"

# Get the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the filename is provided
if (length(args) == 0) {
  stop("Please provide the filename of the file list as a command-line argument.")
}

# Read the filename from the command-line arguments
file_list_filename <- args[1]

# Read the list of filenames from file_list.txt
# file_list <- readLines("file_list_ina.txt")
file_list <- readLines(file_list_filename)

# Initialize an empty character vector to store the output
output <- character()

# Variable to store the energy from the special identifier file
wt_energy <- NULL

# First, find the special identifier's energy
for (i in file_list) {
  isWT <- grepl(wt_identifier, i, ignore.case = TRUE)
  if (isWT) {
    # Read the file content
    lines <- readLines(i)
    
    # Find lines containing ' bk_tot'
    bk_tot_lines <- grep(' bk_tot', lines, value = TRUE)
    
    # Extract the first matching line and parse the energy
    if (length(bk_tot_lines) > 0) {
      first_line <- bk_tot_lines[1]
      wt_energy <- as.numeric(strsplit(first_line, "\\s+")[[1]][3])  # Assuming the energy is the third field
    }
    sprintf("WT energy is %f",wt_energy)
    break
  }
}

# Check if the special identifier energy was found
if (is.null(wt_energy)) {
  stop("WT energy not found!")
}

# Loop over each file
for (i in file_list) {
  # Print the file name
  cat(i, "\n")
  
  # Read the file content
  lines <- readLines(i)
  
  # Find lines containing ' bk_tot'
  bk_tot_lines <- grep(' bk_tot', lines, value = TRUE)
  
  # Extract the first matching line
  if (length(bk_tot_lines) > 0) {
    first_line <- bk_tot_lines[1]
    energy <- as.numeric(strsplit(first_line, "\\s+")[[1]][3])  #  Assuming the energy is the third field
    
    # Calculate the energy difference
    energy_diff <- energy - wt_energy
    
    # Create the output line by pasting the file name, the first matching line, and the energy difference
    output_line <- paste(i, first_line, energy_diff, sep = "    ")
    
    # Append the output line to the output vector
    output <- c(output, output_line)
  }
}

# Write the output to energies.txt
writeLines(output, sprintf("energies_%s",file_list_filename))
