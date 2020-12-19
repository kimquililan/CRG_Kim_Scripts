#_______________________________________________________________________
#
# Telomere Distance Calculation 
#
# Calculates the distance from the telomere of the specific significant bins
# 
#_______________________________________________________________________


#Needed for analysis:
#   signif.bins = matrix of the significant bins containing the columns "bin","chr","start","end"


#Make the matrix
t.matrix <- signif.bins[,colnames(signif.bins) %in% c("bin","chr","start","end", "chrom")]
t.matrix$telomere.dist <- 0
t.matrix$in_telomere <- 0

#Get the chr length, change the format of chr number
chrom.bins.length <- read.table("chrom.bins.1Mb.txt")

#Calculate telomere distance
for(i in 1:dim(t.matrix)[1]){
  row <- match(t.matrix[i,"chrom"],chrom.bins.length[,1])
  chr.length <- chrom.bins.length[row,2]
  bin.start <- t.matrix[i,"start"]
  bin.end <- t.matrix[i,"end"]
  dist1 = (bin.start+bin.end)/2
  dist2 = abs(chr.length*1e06-dist1)
  dist = min(dist1,dist2)
  t.matrix[i,"telomere.dist"] <- dist
  cutoff.max <- chr.length-10
  if(dist/1e6<10 | (dist/1e06>=cutoff.max & dist/1e06<chr.length)){
    t.matrix[i,"in_telomere"] <- TRUE
  }
  else{
    t.matrix[i,"in_telomere"] <- FALSE
  }
  rm(row,chr.length,bin.start,bin.end,dist1,dist2,dist)
}



#-------------------------Results checking-----------------------------

sum(t.matrix[,"in_telomere"])


#For downstream analysis

#For Permutation Test of Telomeric Bins
t.matrix

