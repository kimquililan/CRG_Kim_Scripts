#_______________________________________________________________________
#
# Permutation Test for Telomeres
#
# Performs a permutation test on specific bins
# 
#_______________________________________________________________________

#Needed for analysis
#   t.matrix = matrix of distances from the telomere of specific bins; must contain columns "chr"


#count signif bins per chromosome and store in a vector
length.signif.bins <- rep(0,23)

for(i in 1:23){
  length.signif.bins[i] <- length(which(t.matrix[,2]==i))
}


#Load IF Matrix and exclude chr11,14,X
setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/Raw")
load("Compiled Matrices_0_1Mb_raw.Rdata")
start <- matrix[,2]*1e6
end <- (matrix[,2]+1)*1e6
matrix <- cbind(matrix,start,end)
rows.wo.chr11.14.X <- which(!matrix[,34]==11 & !matrix[,34]==14 & !matrix[,34]==23)
matrix <- matrix[rows.wo.chr11.14.X,]

#Function to random sample for each chromosome
randomTest <- function(chr){
  rows <- which(matrix[,34]==chr)
  random_rows <- matrix[sample(rows,length.signif.bins[chr],replace=FALSE),]
  chrom.bins <- read.table("chrom.bins.1Mb.txt")
  chr.length <- chrom.bins[chr,2]
  tmp.matrix <- matrix(0, nrow=length.signif.bins[chr], ncol=4)
  colnames(tmp.matrix) <- c("chr","start","end","telomere.dist")
  
  tmp.matrix[,1] <- chr
  
  #Calculate telomere distance
  {if(length.signif.bins[chr]>1){
    for (i in 1:length.signif.bins[chr]){
      bin.start <- random_rows[i,35]
      tmp.matrix[i,2] <-bin.start
      bin.end <- random_rows[i,36]
      tmp.matrix[i,3] <-bin.end
      dist1 = (bin.start+bin.end)/2
      dist2 = abs(chr.length*1e06-dist1)
      dist = min(dist1,dist2)
      tmp.matrix[i,4] <- dist
    }
    telomere.dist <- tmp.matrix[,4]
  }
    else{
      bin.start <- random_rows[35]
      tmp.matrix[1,2] <-bin.start
      bin.end <- random_rows[36]
      tmp.matrix[1,3] <-bin.end
      dist1 = (bin.start+bin.end)/2
      dist2 = abs(chr.length*1e06-dist1)
      dist = min(dist1,dist2)
      tmp.matrix[1,4] <- dist
      telomere.dist <- tmp.matrix[1,4]
    }}
  
}

#Running the function
n=1e4
cycle <- matrix(0, nrow=sum(length.signif.bins), ncol=n)

for (k in 1:n){
  random_numbers <- randomTest(which(length.signif.bins>0)[1])
  list.chromosomes<-which(length.signif.bins>0)[-1] #list all chromosomes that are in the list, except chr 1
  
  for (j in 1:length(list.chromosomes)){
    tmp_random_numbers <- randomTest(list.chromosomes[j])
    random_numbers <- c(random_numbers, tmp_random_numbers)
    rm(tmp_random_numbers)
  }
  cycle[,k] <- as.numeric(random_numbers)
}

fraction_telomeres <- matrix(0, nrow=n, ncol=1)
colnames(fraction_telomeres) = "fraction_per_iteration"
for (h in 1:n){
  fraction_telomeres[h,] <- length(which(cycle[,h] <= 10000000))
}


#----------------Results Checking and Plotting------------------------------------
observed_value = sum(t.matrix[,"in telomere?"] == 1)
expected_value = round(mean(fraction_telomeres[,1]),digits=1)
p_val = (length(which(fraction_telomeres[,1] >= observed_value))+1)/(n+1)

hist(fraction_telomeres, xlab="Number of Telomeric Bins per Iteration", main=c("MCL loss",paste("p =",signif(p_val, digits=3)),paste("expected = ",expected_value)))
abline(v=c(expected_value, observed_value), col=c("red","forestgreen"), lty=2)



