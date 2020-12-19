#________________________________________________________________________
#
# BIND MATRIX
#
# Function to bind IF matrices of all chromosomes for one viewpoint
#________________________________________________________________________


Bind_Matrices <- function(file = "Compiled Matrices.txt", 
                          viewpt, 
                          resolution){
  
# ...generates matrix of IF values for all chromosomes of one viewpoint
#
# Args:
#     compiled.filenames: Text file containing filenames of generated matrices for each chromosome
# Returns:
#     matrix containing values for all chromosomes for one viewpoint
  
  filenames <- read.table(file)
  load(paste(filenames[1,viewpt+1]))
  matrix1 <- matrix
  
  for(i in 2:23){
    load(paste(filenames[i,viewpt+1]))
    matrix2 <- matrix
    
    matrix1 <- rbind(matrix1, matrix2)
    rm(matrix2)
  }
  
  matrix <- matrix1
  
  save(matrix, file=paste("Compiled Matrices_", viewpt, "_", resolution, ".Rdata",sep = ""))
  print("Matrix is saved")
}






