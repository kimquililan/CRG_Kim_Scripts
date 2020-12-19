#_______________________________________________________________________
#
# Chromosome Length vs Raw Counts (Boxplot)
#
# Calculates Fraction of raw counts for each sample for every chromosome
# It also generates a boxplot and calculates the significance through
# Wilcoxon Test
# 
#_______________________________________________________________________

rm(list=ls())
setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/Raw")
load("Compiled Matrices_0_1Mb_raw.Rdata")

list.chr <- c(1:10,12:13,15:22)
load("Sum of Interactions per chrom and sample.Rdata")

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

intxn.sum.chr.sample <- intxn.sum.chr.sample[,c(control,experimental)]
intxn.sum.chr.sample2 <- matrix(0, nrow=nrow(intxn.sum.chr.sample), ncol=ncol(intxn.sum.chr.sample))
colnames(intxn.sum.chr.sample2) <- c(control,experimental)

#Calculate fraction of raw counts
for(i in 1:nrow(intxn.sum.chr.sample)){
  for(j in 1:ncol(intxn.sum.chr.sample)){
    intxn.sum.chr.sample2[i,j]<-intxn.sum.chr.sample[i,j]/sum(intxn.sum.chr.sample[,j])
  }
}

#Statistical Test
wilcox.output <- matrix(0, nrow=20, ncol=3)
wilcox.output[,1] <- list.chr
colnames(wilcox.output) <- c("chr", "p-value","p-adjust")

for (i in 1:20){
  bin.experimental <- as.numeric(intxn.sum.chr.sample2[i,experimental])
  bin.control <- as.numeric(intxn.sum.chr.sample2[i,control])
  
  w.test <- wilcox.test(bin.experimental,bin.control,alternative="two.sided",na.rm=TRUE)
  wilcox.output[i,2] <- w.test$p.value  
}

padjust <- p.adjust(wilcox.output[,"p-value"], method="BH")
wilcox.output[,"p-adjust"] <- padjust


#-----------------------------------Plotting-----------------------------------------

  IF_control <- intxn.sum.chr.sample2[i,control]
  IF_experimental <- intxn.sum.chr.sample2[i,experimental]
  chrom_data <- cbind(rep(1,sum(length(c(control,experimental)))),
                      c(rep("control",length(control)),rep("experimental",length(experimental))),
                      c(IF_control,IF_experimental))
  
  for (i in 2:20){
    IF_control <- intxn.sum.chr.sample2[i,control]
    IF_experimental <- intxn.sum.chr.sample2[i,experimental]
    tmp <- cbind(rep(list.chr[i],sum(length(c(control,experimental)))),
                 c(rep("control",length(control)),rep("experimental",length(experimental))),
                 c(IF_control,IF_experimental))
    chrom_data <- rbind(chrom_data,tmp)
    
    rm(IF_control, IF_experimental,tmp)
  }

  
  chrom_data <- data.frame(chrom_data, row.names = NULL)
  colnames(chrom_data) <- c("chr","cellgrp","ratio_raw")
  chrom_data[,"ratio_raw"] <- as.numeric(as.character(chrom_data[,"ratio_raw"]))
  chrom_data[,"chr"] <- factor(chrom_data[,"chr"], levels=list.chr)
  
  
  library(ggplot2)
  ggplot(chrom_data, aes(x= chr, y=ratio_raw, colour=cellgrp, fill=cellgrp)) +
    xlab("Chromosome") + ylab("Raw Counts (Fraction)") +
    geom_boxplot() +
    theme_classic() +
    scale_color_manual(values=(rep("black",4))) +
    scale_fill_manual(values=c("grey","orange"))
  
  