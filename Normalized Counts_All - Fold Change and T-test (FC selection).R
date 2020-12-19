#_______________________________________________________________________
#
# Normalized Counts - Fold Change and T-test
# All viewpoints
#
# Calculates the fold change for normalized IF values and filters by selecting bins
# wih >= 1.5-FC then performs T-test on selected bins
# 
#_______________________________________________________________________

Ttest_Norm = function(viewpt_file, cell_group){
  
  setwd("/Users/kimquililan/Google Drive/Data/Normalized2/viewpoints_1Mb_norm")
  load(paste0(viewpt_file,".Rdata"))
  matrix[,"start"] <- matrix[,"start"]*1e6
  matrix[,"end"] <- matrix[,"end"]*1e6
  matrix <- subset(matrix, chrom %in% c(1:10,12:13,15:22))
  
  
  #Choose between MCL or CLL as experimental
  ifelse(cell_group == "MCL", 
         experimental <- c("cMCL1064","cMCL568","nnMCL309","nnMCL817", "nnMCL828"),
         experimental <- c("mCLL3","mCLL110","mCLL1228","mCLL1525","mCLL1532","uCLL12","uCLL182"))
  
  control = c("NBC1","NBC2","NBC3","GCBC1","GCBC2","GCBC3","MBC1","MBC2","MBC3","PBC1","PBC2","PBC3")
  
  #filter rows based on number of NAs
  NAs.exp<-rowSums(is.na(matrix[,experimental]))
  NAs.control<-rowSums(is.na(matrix[,control]))
  rows.to.include<-which(NAs.exp<=(length(experimental)-3) & NAs.control<=(length(control)-3))
  matrix<-matrix[rows.to.include,]
  

  #log transform
  log_matrix <- log(matrix[,c(experimental,control)]+2^(-6),2)
  new_matrix <- data.frame(matrix[,c("bin","chrom","start","end")],log_matrix)

  #Calculate log2(FC)
  new_matrix$fold_change = apply(new_matrix, 1, function(x) mean(x[experimental],na.rm=T) - mean(x[control],na.rm=T))
  
  #select bins based on log2FC as input matrix for wilcoxon test
  new_matrix_input = subset(new_matrix, abs(fold_change) > log2(1.5))

  #T Test
  p.value.matrix <- new_matrix_input[,c("bin","chrom","start","end","fold_change")]
  p.value.matrix$p.value =  apply(new_matrix_input, 1, function(x) t.test(x[experimental], x[control],alternative = "two.sided", na.rm=T)$p.value)
  p.value.matrix$p.adjust = p.adjust(p.value.matrix$p.value, method="BH")

  p.value.matrix <<- p.value.matrix
}

#---------For downstream analysis---------------------------------------

#For Calculating Telomere Distance of significant bins
signif.bins <- p.value.matrix[which(p.value.matrix[,6]<0.1),]

#For Calculating MCL-specific or CLL-specific significant bins
MCL.up= subset(p.value.matrix, p.adjust < 0.1 & fold_change > log2(1.5))
MCL.down = subset(p.value.matrix, p.adjust < 0.1 & fold_change < -log2(1.5))

CLL.up= subset(p.value.matrix, p.adjust < 0.1 & fold_change > log2(1.5))
CLL.down = subset(p.value.matrix, p.adjust < 0.1 & fold_change < -log2(1.5))

similar.up = Reduce(intersect, list(MCL.up$bin,CLL.up$bin))
similar.down = Reduce(intersect, list(MCL.down$bin,CLL.down$bin))

MCL.up.spec = subset(MCL.up, !bin %in% similar.up)
MCL.down.spec = subset(MCL.down, !bin %in% similar.down)

CLL.up.spec = subset(CLL.up, !bin %in% similar.up)
CLL.down.spec = subset(CLL.down, !bin %in% similar.down)


#zscore
load("chr11.68.to.69.Rdata")
experimental <- c("cMCL1064","cMCL568","nnMCL309","nnMCL817", "nnMCL828")
experimental <- c("mCLL3","mCLL110","mCLL1228","mCLL1525","mCLL1532","uCLL12","uCLL182")
control = c("NBC1","NBC2","NBC3","GCBC1","GCBC2","GCBC3","MBC1","MBC2","MBC3","PBC1","PBC2","PBC3")

matrix_input = subset(matrix, bin %in% MCL.up$bin)
zscore_results = t(apply(matrix_input[,c(experimental,control)], 1, function(x) (x-mean(x))/sd(x)))
output = data.frame(matrix_input[,c("bin","chrom","start","end")],zscore_results)

#heatmap
library(gplots)
hmcol = colorRampPalette(c("royalblue4","dodgerblue4","dodgerblue3","steelblue3","lightblue3","lavender","white","papayawhip","wheat2","sandybrown","brown2","brown3","brown4"))(100)
color.breaks=seq(-1.5,1.5,length=101)
sample.colors <- c(rep("orange",5),rep("palegreen", 3),rep("red",3),rep("blue",3),rep("mediumpurple2",3))

heatmap.2(zscore_results, distfun = function(x) {dist(x, method="euclidean")},
          hclustfun = function(x) {hclust(x, method="ward.D2")}, trace="none",density.info="none", col = hmcol, breaks=color.breaks, 
          labRow = FALSE, labCol = FALSE, ColSideColors=sample.colors)

#telomere
#Make the matrix
w=4
t.matrix <- matrices_list[[w]]
t.matrix$telomere.dist <- 0
t.matrix$in_telomere <- 0
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

sum(t.matrix$in_telomere)


#permutation
length.signif.bins <- rep(0,23)

for(i in 1:23){
  length.signif.bins[i] <- length(which(t.matrix[,"chrom"]==i))
}


#Load IF Matrix and exclude chr11,14,X
load(paste0(viewpoints[i],".Rdata"))
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
n=1e3
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

hist(fraction_telomeres, xlab="Number of Telomeric Bins per Iteration", main=c("MCL.down.spec",paste("p =",signif(p_val, digits=3)),paste("expected = ",expected_value)))
abline(v=c(expected_value, observed_value), col=c("red","forestgreen"), lty=2)







#-----
