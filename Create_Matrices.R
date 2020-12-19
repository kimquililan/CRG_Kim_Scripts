#________________________________________________________________________
#
# CREATE MATRIX
#
# Function to create a matrix from filenames of IF values for different 
# samples with different viewpoints
#________________________________________________________________________


Create_Matrix<- function(compiled.filenames,
                        chrom,
                        resolution = "1Mb",
                        viewpt = 0) {

  #...generates matrix of values for all bins on one chromosome 
  #
  #Args:
  #   compiled.filenames: text file containing filenames (.abc) 
  #     of all samples for all viewpoints for specified chromosome;
  #     each column is a viewpoint
  #   chrom: chromosome number from 1 to 23
  #   resolution: one of the following (1Mb, 250kb, 500kb)
  #   viewpt: viewpoint from 0 to 2
  #
  #Returns
  #   matrix of IF values for one chromosome for one viewpoint
  
  if (chrom == 23) chromosome <- paste("X") else chromosome <- chrom
  
  
  #make the matrix
  sample_order <- read.table("Sample_order.txt")
  chrom.bins <- read.table(paste0("chrom.bins.",resolution,".txt"))
  nrows <- chrom.bins[chrom,2]
  matrix <- matrix(0, nrow=nrows,ncol=34)
  matrix[,1] <- viewpt
  matrix[,2] <- 0:(nrows-1)
  matrix[,34] <- chrom
  colnames(matrix)<-c("row","col", paste(sample_order[,1],sep=""), "chrom")
  
  #read the files
  chr14_n <- read.table(paste(compiled.filenames))
  
  for(i in 1:31){
    current.file <- read.table(as.character(chr14_n[i,1]))
    current.selection <- current.file[which(current.file[,1]==viewpt),]
    idx <- match(matrix[,2],current.selection[,2])
    matrix[,(2+i)] <- current.selection[idx,3]
    
    #Switch back NAs to zeros
    rows.w.nas <- which(is.na(matrix[,(2+i)]))
    matrix[rows.w.nas,(2+i)] <- 0
    rm(rows.w.nas)
    
    #Determine the badcols and input NA in matrix
    current.file <- read.table(as.character(chr14_n[i,1]), comment.char = "", fill=TRUE)
    badcols <- as.numeric(unlist(strsplit(as.character(current.file[5,3]),",")))
    rows.w.nas <- match(badcols,matrix[,2])
    matrix[rows.w.nas,(2+i)] <- NA
    
    #clean
    rm(badcols)
    rm(current.file)
    rm(current.selection)
    rm(rows.w.nas)
    
  }
  
  #save the file
  save(matrix, file=paste0("chr14_",chromosome, "_matrix_", viewpt,"_", resolution, ".Rdata"))
}

#-------- Loop to create all matrices for all chromosomes for one viewpoint --------

  for (j in 1:23){
    if (j == 23) chromosome <- paste("X") else chromosome <- j
    Create_Matrix(paste("chr14_",chromosome,"_filenames.txt",sep=""),j,"1Mb",0)
    print(paste("chr14_",j, "_matrix_", 0,"_1Mb.Rdata is created", sep=""))
  }

#-----for multiple viewpoints-------
for(i in 0:2){
  for (j in 1:23){
    if (j == 23) chromosome <- paste("X") else chromosome <- j
    Create_Matrix(paste("chr14_",chromosome,"_filenames.txt",sep=""),j,"250kb",i)
    print(paste("chr14_",j, "_matrix_", i, "_250kb.Rdata is created", sep=""))
  }
}


