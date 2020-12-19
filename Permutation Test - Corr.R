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
  length.signif.bins[i] <- length(which(t.matrix[,"chrom"]==i))
}


#Load IF Matrix and exclude chr11,14,X
load("chr14_105_106.Rdata")
matrix[,"start"] <- matrix[,"start"]*1e6
matrix[,"end"] <- matrix[,"end"]*1e6
matrix <- subset(matrix, chrom %in% c(1:10,12:13,15:22))

#Function to random sample for each chromosome
randomTest <- function(chr){
  rows <- which(matrix[,"chrom"]==chr)
  random_rows <- matrix[sample(rows,length.signif.bins[chr],replace=FALSE),]
  chrom.bins <- read.table("chrom.bins.1Mb.txt")
  chr.length <- chrom.bins[chr,2]
  tmp.matrix <- matrix(0, nrow=length.signif.bins[chr], ncol=4)
  colnames(tmp.matrix) <- c("chr","start","end","telomere.dist")
  
  tmp.matrix[,"chr"] <- chr
  
  #Calculate telomere distance
  {if(length.signif.bins[chr]>1){
    for (i in 1:length.signif.bins[chr]){
      bin.start <- random_rows[i,"start"]
      tmp.matrix[i,"start"] <-bin.start
      bin.end <- random_rows[i,"end"]
      tmp.matrix[i,"end"] <-bin.end
      dist1 = (bin.start+bin.end)/2
      dist2 = abs(chr.length*1e06-dist1)
      dist = min(dist1,dist2)
      tmp.matrix[i,"telomere.dist"] <- dist
    }
    telomere.dist <- tmp.matrix[,"telomere.dist"]
  }
    else{
      bin.start = random_rows["start"]
      tmp.matrix[1,"start"] <- as.numeric(bin.start)
      bin.end <- random_rows["end"]
      tmp.matrix[1,"end"] <- as.numeric(bin.end)
      dist1 = (bin.start+bin.end)/2
      dist2 = abs(chr.length*1e06-dist1)
      dist = min(dist1,dist2)
      tmp.matrix[1,"telomere.dist"] <- dist
      telomere.dist <- tmp.matrix[1,"telomere.dist"]
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

fraction_telomeres <- c()
for (h in 1:n){
  fraction_telomeres[h] <- length(which(cycle[,h] <= 10000000))
}


#----------------Results Checking and Plotting------------------------------------
observed_value = sum(t.matrix[,"in_telomere"] == 1)
expected_value = round(mean(fraction_telomeres),digits=1)
p_val = (length(which(fraction_telomeres >= observed_value))+1)/(n+1)

hist(fraction_telomeres, xlab="Number of Telomeric Bins per Iteration", main=c("MCL gain",paste("p =",signif(p_val, digits=3)),paste("expected = ",expected_value)))
abline(v=c(expected_value, observed_value), col=c("red","forestgreen"), lty=2)



