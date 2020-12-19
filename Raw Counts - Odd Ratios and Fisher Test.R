#_______________________________________________________________________
#
# Raw counts - Odds Ratio and Fisher Test
#
# Calculates Relative Risk (RR) for raw counts and performs Fisher Test to calculate
# significance. Additional: Creates volcano plots
# 
#_______________________________________________________________________

rm(list=ls())
setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/Raw")
load("Compiled Matrices_0_1Mb_raw.Rdata")
start <- matrix[,2]*1e6
end <- (matrix[,2]+1)*1e6
matrix[,2] <- 0:(dim(matrix)[1]-1)
matrix <- cbind(matrix,start,end)
rows.wo.chr11.14.X <- which(!matrix[,34]==11 & !matrix[,34]==14 & !matrix[,34]==23)
matrix <- matrix[rows.wo.chr11.14.X,]

#Choose between MCL or CLL as experimental
experimental = c("MCL-cMCL1",
                 "MCL-cMCL2",
                 "MCL-nnMCL1",
                 "MCL-nnMCL2",
                 "MCL-nnMCL3")

experimental = c("CLL-uCLL1","CLL-uCLL2","CLL-mCLL1","CLL-mCLL2","CLL-mCLL3","CLL-mCLL4","CLL-mCLL5")

control = c("NBC-NBC1","NBC-NBC2","NBC-NBC3",
            "GCBC-GCBC1","GCBC-GCBC2","GCBC-GCBC3",
            "MBC-MBC1","MBC-MBC2","MBC-MBC3",
            "PBC-PBC1","PBC-PBC2","PBC-PBC3")

#filter based on CNA
CNA = read.table("CNAs.txt", header=T)

rows.excluded = list()
for (i in 1:nrow(CNA)){
  positions = seq(CNA[i,"start"], CNA[i,"end"]-1, by=1e6)
  chr = CNA[i,"chrom"]
  rows.excluded[[i]] = which(matrix[,"start"] %in% positions & matrix[,"chrom"]==chr)
}

matrix = matrix[-unlist(rows.excluded),]

#filter rows based on number of NAs
NAs.exp<-rowSums(is.na(matrix[,experimental]))
NAs.control<-rowSums(is.na(matrix[,control]))
rows.to.include<-which(NAs.exp<=(length(experimental)-3) & NAs.control<=(length(control)-3))
matrix<-matrix[rows.to.include,]


#Calculation 
output <- matrix(0, nrow=nrow(matrix), ncol=7)
colnames(output) <- c("bin","chr","start","end","Odd Ratio","p-value", "p-adjust")
output[,"bin"] <- matrix[,2]
output[,c("chr","start", "end")] <- matrix[,34:36]

colSum = data.frame(apply(matrix[,c(control,experimental)], 2, sum, na.rm=T))

for (i in 1:nrow(matrix)){ 
  tmp <- matrix(0,nrow=length(c(experimental,control)), ncol=2)
  rownames(tmp) <- c(control, experimental)
  tmp[,1] = colSum[c(control,experimental),1]
  tmp[control,2] = as.numeric(unlist(matrix[i,control]))
  tmp[experimental,2] <- as.numeric(unlist(matrix[i,experimental]))
  tmp <- tmp+1
  a = sum(tmp[experimental,2], na.rm=TRUE)
  b = sum(tmp[experimental,1], na.rm=TRUE)  
  c = sum(tmp[control,2], na.rm=TRUE)
  d = sum(tmp[control,1], na.rm = TRUE)  
  OR =   (a*d)/(b*c)
  output[i,"Odd Ratio"] <- OR
  fishertest.results <- fisher.test(matrix(c(a,b,c,d),nrow=2, ncol=2))
  output[i,"p-value"] <- fishertest.results$p.value
  rm(a,b,c,d,OR,tmp)
}


padjust <- p.adjust(output[,"p-value"], method="BH")
output[,"p-adjust"] <- padjust

#--------Results checking-------
# Expected Results:
# MCL vs normal: 66 significant bins
# CLL vs normal: 164 significant bins

#Significant Bins
length(which(output[,7]<0.1))
length(which(output[,5]>1 & output[,7]<0.1))
length(which(output[,5]<1 & output[,7]<0.1))
length(which(output[,7]>0.1))
output[which(output[,7]<0.1),]

#---------Plotting-------------------
#volcano plot
interchrom <- output
OR <- interchrom[,"Odd Ratio"]
pval <- interchrom[,"p-adjust"]
vplotdata <- cbind(pval,OR)

plot(log2(vplotdata[,2]), -log10(vplotdata[,1]),
     xlim=c(-3,3),ylim=c(0,10),xlab="LOR", ylab="-log10(FDR)",
     pch=20,cex=1, col="grey")
points(log2(vplotdata[which(vplotdata[,"OR"]>1 & vplotdata[,"pval"]<0.1),2]),
       -log10(vplotdata[which(vplotdata[,"OR"]>1 & vplotdata[,"pval"]<0.1),1]), col="red", pch=20, cex=1)
points(log2(vplotdata[which(vplotdata[,"OR"]<1 & vplotdata[,"pval"]<0.1),2]),
       -log10(vplotdata[which(vplotdata[,"OR"]<1 & vplotdata[,"pval"]<0.1),1]), col="blue", pch=20, cex=1)
abline(h=1, v=c(-0.1,0.1),lty=2)

#LOR vs raw cts
plot(output[,"total raw cts"],
     log2(output[,"Odd Ratio"]),
     pch=20, xlab="total raw cts",ylab="LOR",main="2515 significant bins")


#------Getting Data for Downstream Analysis--------
#Distance to telomere
#Get the bin, chr, start and end
inc <- output[which(output[,"Odd Ratio"]>1 & output[,"p-adjust"]<0.1),1:4]
dec <- output[which(output[,"Odd Ratio"]<1 & output[,"p-adjust"]<0.1),1:4]
stable <- output[which(output[,"p-adjust"]>0.1),1:4]


#Chromatin states
#get the bin number
inc <- output[which(output[,"Odd Ratio"]>1 & output[,"p-adjust"]<0.1),"bin"]
dec <- output[which(output[,"Odd Ratio"]<1 & output[,"p-adjust"]<0.1),"bin"]
stable <- output[which(output[,"p-adjust"]>0.1),"bin"]

