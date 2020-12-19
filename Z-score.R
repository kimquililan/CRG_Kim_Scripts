#_______________________________________________________________________
#
# Z-score
#
# Calculates the Row z-score 
# 
#_______________________________________________________________________

#Needed for analysis:
#   selected.bins = list of sample names; specified below
#   signif.rows = vector of bin numbers
#
#Returns:
#   output: which can be used to create a heatmap


setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/Normalized")
load("Compiled Matrices_0_1Mb_norm.Rdata")
selected.bins <- c(control,experimental) #experimental may change to either MCL or CLL set
matrix[,2] <- 0:(dim(matrix)[1]-1) 
output <- matrix(0, nrow=length(signif.rows),ncol=length(selected.bins)+1)
colnames(output)<- c("bin",colnames(matrix[,selected.bins]))



for(j in 1:length(signif.rows)){
  
  x <- matrix[signif.rows[j]+1,selected.bins]
  output[j,1] <- matrix[signif.rows[j]+1,2]
  mean.sample <- mean(as.numeric(x),na.rm=TRUE)
  sd.sample <- sd(as.numeric(x), na.rm=TRUE)
  
  x.z <- rep(0,length(selected.bins))
  
  for(i in 1:length(selected.bins)){
    x.z[i] <- (x[i]-mean.sample)/sd.sample
  }
  
  output[j,selected.bins] <- x.z
  rm(x)
}

