#_______________________________________________________________________
#
# Chromosome Length vs Odd Ratio (Scatterplot)
#
# Calculates Odds Ratio for raw counts and performs Fisher Test to calculate
# significance. Also plots the LOR for each sample in each chromosome
# 
#_______________________________________________________________________

rm(list=ls())
setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/Raw")
load("Compiled Matrices_0_1Mb_raw.Rdata")

#choose which experimental or control
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


load("Sum of Interactions per chrom and sample.Rdata")


OR_matrix <- matrix(0, nrow=20, ncol=length(experimental))
colnames(OR_matrix) <- experimental
intxn.sum.chr.sample <- intxn.sum.chr.sample+1

for (i in 1:20){
  for (j in 1:length(experimental)){
    a = intxn.sum.chr.sample[i,experimental[j]]
    b = sum(intxn.sum.chr.sample[,experimental[j]]) - a
    c = sum(intxn.sum.chr.sample[i,control])
    d = sum(intxn.sum.chr.sample[,control]) - c
    OR =   (a*d)/(b*c)
    OR_matrix[i,j] <- OR
    rm(a,b,c,d,OR)
  }
}

list.chr <- c(1:10,12:13,15:22)

chromosomes <- as.numeric(sapply(list.chr, function(x) rep(x,length(experimental))))
cell.grp <- rep("CLL",(20*length(experimental))) #change either "MCL" or "CLL"
Odd_Ratio <- as.numeric(t(OR_matrix))
output <- data.frame(cell.grp, chromosomes, Odd_Ratio)


#Get both output for CLL and MCL
MCL.output <- output
CLL.output <- output
final_output <- rbind(MCL.output, CLL.output)
final_output[,"chromosomes"] <- factor(final_output[,"chromosomes"], levels=list.chr)

#----------------------------Plotting-------------------------------------------

#Colored according to sample
  library(ggplot2)
  ggplot(final_output, aes(x=cell.grp, y=log2(Odd_Ratio), fill=cell.grp)) +
    xlab("Cell") + ylab("Log2 Odd Ratio") + ylim(-0.6,0.6) +
    geom_point(shape=21, size=1.75) +
    theme_classic() + theme(strip.background = element_blank(),strip.text.x = element_blank()) + facet_wrap(~chromosomes, ncol=1, drop=FALSE)+coord_flip() +
    scale_fill_manual(values=c("azure4","orange")) +
    geom_hline(yintercept=0, lty=2)

#Colored according to increased or decreased OR
remark <- sapply(final_output[,"Odd_Ratio"], function(x) ifelse(x>1, "increased", "decreased"))
final_output <- cbind(final_output, remark)
final_output[,"remark"] <- factor(final_output[,"remark"], levels=c("increased", "decreased"))

  #Either MCL or CLL
  library(ggplot2)
  ggplot(final_output[which(cell.grp=="CLL"),], aes(x= chromosomes, y=log2(Odd_Ratio), fill=remark)) +
    xlab("Chromosome") + ylab("Log2 Odd Ratio") + ylim(-0.6,0.6) + 
    geom_point(shape=21, size=1.75) +
    theme_classic() + coord_flip() +
    scale_fill_manual(values=c("red","blue")) +
    ggeom_hline(yintercept=0, lty=2)
  
  #MCL and CLL in one plot
  ggplot(final_output, aes(x=cell.grp, y=log2(Odd_Ratio), fill=remark)) +
    xlab("Cell") + ylab("Log2 Odd Ratio") + ylim(-0.6,0.6) +
    geom_point(shape=21, size=1.75) +
    theme_classic() + theme(strip.background = element_blank(),strip.text.x = element_blank()) + facet_wrap(~chromosomes, ncol=1, drop=FALSE)+coord_flip() +
    scale_fill_manual(values=c("red","blue")) + 
    geom_hline(yintercept=0, lty=2) 