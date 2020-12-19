#_______________________________________________________________________
#
# Normalized Counts - Fold Change and T-test
#
# Calculates the fold change for normalized IF values and filters by selecting bins
# wih >= 1.5-FC then performs T-test on selected bins
# 
#_______________________________________________________________________


#Wilcoxon Test/T-test of Norm IFs with log2FC filtering
rm(list=ls())
setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/Normalized")
load("Compiled Matrices_0_1Mb_norm.Rdata")

experimental = c("MCL-cMCL1",
                 "MCL-nnMCL2",
                 "MCL-cMCL2",
                 "MCL-nnMCL1",
                 "MCL-nnMCL3")

experimental = c("CLL-uCLL1","CLL-uCLL2","CLL-mCLL1","CLL-mCLL2","CLL-mCLL3","CLL-mCLL4","CLL-mCLL5")

control = c("NBC-NBC1","NBC-NBC2","NBC-NBC3",
            "GCBC-GCBC1","GCBC-GCBC2","GCBC-GCBC3",
            "MBC-MBC1","MBC-MBC2","MBC-MBC3",
            "PBC-PBC1","PBC-PBC2","PBC-PBC3")

start <- matrix[,2]*1e6
end <- (matrix[,2]+1)*1e6
matrix[,2] <- 0:(dim(matrix)[1]-1)
matrix <- cbind(matrix,start,end)
rows.wo.chr11.14.X <- which(!matrix[,34]==11 & !matrix[,34]==14 & !matrix[,34]==23)
matrix <- matrix[rows.wo.chr11.14.X,]

#log transform
matrix2 <- matrix+2^(-4)
log_matrix <- log(matrix2[,3:26],2)
new_matrix <- cbind(matrix[,c(2,34:36)],log_matrix)

#filter rows based on nr. of NAs in MCL and and other
NAs.exp<-rowSums(is.na(new_matrix[,experimental]))
NAs.control<-rowSums(is.na(new_matrix[,control]))
rows.to.include<-which(NAs.exp<=(length(experimental)-3) & NAs.control<=(length(control)-3))
new_matrix<-new_matrix[rows.to.include,]

#Calculate log2(FC)
output <- matrix(0, nrow=dim(new_matrix)[1], ncol=5)
output[,1:4] <- new_matrix[,1:4]
colnames(output) <- c(colnames(new_matrix[,1:4]),"fold change")

for (i in 1:nrow(new_matrix)){
  mean_control <- mean(as.numeric(new_matrix[i,control]),na.rm=TRUE)
  mean_exp <- mean(as.numeric(new_matrix[i,experimental]),na.rm=TRUE)
  
  output[i,"fold change"] <- mean_exp - mean_control
}

#select bins based on log2FC as input matrix for wilcoxon test
rows.with.high.FC<-which(abs(output[,"fold change"])>log2(1.5))
new_matrix_input_wilcoxon<-new_matrix[rows.with.high.FC,]
output2 <- output[rows.with.high.FC,]

#Wilcoxon Test
rows <- nrow(new_matrix_input_wilcoxon)
withTrans <- new_matrix_input_wilcoxon[,experimental]
withoutTrans <- new_matrix_input_wilcoxon[,control]

p.value.matrix <- matrix(0, nrow=rows,ncol=7)
p.value.matrix[,1:4] <- new_matrix_input_wilcoxon[,1:4]
colnames(p.value.matrix)<-c("bin","chr","start","end","p-value","p-adjust","fold change")


for (i in 1:rows){
  bin.w <- as.numeric(withTrans[i,])
  bin.wo <- as.numeric(withoutTrans[i,])
  w.test <- t.test(bin.w,bin.wo,alternative="two.sided",na.rm=TRUE)
  p.value.matrix[i,5] <- w.test$p.value
  p.value.matrix[i,7] <- output2[i,5]
  rm(w.test, bin.w,bin.wo)
}

padjust <- p.adjust(p.value.matrix[,5], method="BH")
p.value.matrix[,6] <- padjust



#--------------------------Results checking-----------------------------
# Expected Results:
# MCL vs normal: 125 significant bins
# CLL vs normal: 227 significant bins

#Significant Bins
length(which(p.value.matrix[,6]<0.1))



#---------For downstream analysis---------------------------------------

#For Calculating Telomere Distance of significant bins
signif.bins <- p.value.matrix[which(p.value.matrix[,6]<0.1),]

#For Calculating MCL-specific or CLL-specific significant bins
results.MCL = p.value.matrix
results.CLL = p.value.matrix
similar_bins = Reduce(intersect, list(results.MCL[which(results.MCL[,"p-adjust"]<0.1),1],results.CLL[which(results.CLL[,"p-adjust"]<0.1),"bin"]))
MCL.specific.bins = results.MCL[which(!results.MCL[,"bin"] %in% similar_bins & results.MCL[,"p-adjust"]<0.1),]
CLL.specific.bins = results.CLL[which(!results.CLL[,"bin"] %in% similar_bins & results.CLL[,"p-adjust"]<0.1),]


#Telomere Distance of Bins with either Gain or Loss 
#change if greater than or equal to 0, and matrix used
#followed by "Telomere Distance Calculation.Rdata" then by Permutataion Test(optional)
signif.bins <- MCL.specific.bins[which(MCL.specific.bins[,"fold change"]<0),1:4] 

#For calculating row zscore then heatmap
signif.rows <- MCL.specific.bins[,1]
