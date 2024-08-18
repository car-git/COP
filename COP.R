# Modify as needed:
file_path <- "/Users/carlosmh1/Desktop/Barley.txt"
file_path <- "/Users/carlosmh1/Desktop/kemp.txt"

# Read the file, skipping lines that start with "#" and using tab as the delimiter
Data <- read.table(file_path, header = FALSE, sep = "\t", comment.char = "#")

# Assign column names
colnames(Data) <- c("STRING", "PARENT_1", "PARENT_2")

# Display the Data
print(Data)

# Extract all individuals from the Data table
all_individuals <- unlist(Data)

# Get unique individuals and sort them
strains <- sort(unique(all_individuals))

# Display the strains list
print(strains)

# Count the number of unique strains
n <- length(strains)

# Create the square matrix P of size n filled with zeros
P <- matrix(0, nrow = n, ncol = n)

# Set the row and column names of P to the strain names
rownames(P) <- strains
colnames(P) <- strains

# Populate the matrix P
for (i in 1:n) {
  # Get the current strain name
  current_strain <- strains[i]
  
  # Find rows in Data where current_strain is a parent
  parent_indices <- which(Data$PARENT_1 == current_strain | Data$PARENT_2 == current_strain)
  
  # For each found row, set P(i, j) = 1 where j is the index of the strain that has i as a parent
  for (j in parent_indices) {
    # Find the index of the strain that has current_strain as a parent
    child_strain_index <- which(strains == Data$STRING[j])
    
    # Set P(i, child_strain_index) = 1
    P[i, child_strain_index] <- 1
  }
}

# Display the populated matrix P with row and column names
print(P)

# Calculate Q matrix
Q <- solve(diag(nrow(P)) - P / 2)

# Create matrix E with the same row and column names as P
E <- matrix(0, nrow = nrow(P), ncol = ncol(P))
rownames(E) <- rownames(P)
colnames(E) <- colnames(P)

# Extract the strains in the first column "STRING"
strains_in_string <- Data$STRING

# Find strains in the "strains" list that are not in the "STRING" column
strains_not_in_string <- setdiff(strains, strains_in_string)

# Iterate through each strain in strains_not_in_string
for (strain in strains_not_in_string) {
  # Find the index of the strain in the strains list
  strain_index <- which(strains == strain)
  
  # Set the corresponding diagonal element in E to 1
  E[strain_index, strain_index] <- 1
}

# Display the updated matrix E
print(E)


# Calculate matrices F and cop
F   <- E %*% Q
cop <- (1/2) * t(F) %*% F

# If terminal ancestors are homozygous use this:
# cop <- 2 * cop

# Extract the directory from the file path
directory <- dirname(file_path)

# Construct the full path for saving the sorted_cop matrix (modify as needed):
output_file <- file.path(directory, "cop_matrix.csv")

# Save the sorted_cop matrix to a CSV file in the extracted directory
write.csv(cop, output_file, row.names = TRUE)

# Print the path to the saved file
print(paste("File saved at:", output_file))
