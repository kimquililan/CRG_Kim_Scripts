#_______________________________________________________________________
#
# Chromosome Length vs Odd Ratio (Barplot)
#
# Calculates Odds Ratio for raw counts and performs Fisher Test to calculate
# significance. Also plots the LOR for each chromosome
# 
#_______________________________________________________________________


rm(list=ls())
setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/Raw")
load("Compiled Matrices_0_1Mb_raw.Rdata")

#choose experimental and the control
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


#Calculate
list.chr <- c(1:10,12:13,15:22)
output <- matrix(0, nrow=20, ncol=5)
colnames(output) <- c("chr","Odd Ratio","remark","p-value", "p-adjust")
output[,"chr"] <- list.chr

load("Sum of Interactions per chrom and sample.Rdata")


for (i in 1:20){
  a = sum(intxn.sum.chr.sample[i,experimental])
  b = sum(intxn.sum.chr.sample[,experimental]) - a
  c = sum(intxn.sum.chr.sample[i,control])
  d = sum(intxn.sum.chr.sample[,control]) - c
  OR =   (a*d)/(b*c)
  output[i,"Odd Ratio"] <- OR
  output[i,"remark"] <- ifelse (OR > 1, "increased","decreased")
  tmp.output <- matrix(c(a,b,c,d),nrow=2, ncol=2)
  fishertest.results <- fisher.test(tmp.output)
  output[i,"p-value"] <- fishertest.results$p.value
  rm(a,b,c,d,OR,tmp.output)
}

padjust <- p.adjust(output[,"p-value"], method="BH")
output[,"p-adjust"] <- padjust
output <- data.frame(output)
output[,"Odd.Ratio"] <- as.numeric(as.character(output[,"Odd.Ratio"]))
output[,"remark"] <- factor(output[,"remark"], levels=c("increased","decreased"))
output[,"chr"] <- factor(output[,"chr"],levels=c(22,21,20,19,18,17,16,15,13,12,10,9,8,7,6,5,4,3,2,1))




#----------------------Plotting------------------------------

#For each sample either MCL or CLL
library(ggplot2)
gb <- ggplot(output, aes(x= chr, y=log2(Odd.Ratio), fill=remark))
gb <- gb + xlab("Chromosome") + ylab("Log2 Odd Ratio") +  ylim(-0.4,0.4)
gb <- gb + geom_bar(stat="identity")
gb <- gb + theme_classic() + coord_flip()
gb <- gb+scale_fill_manual(values=c("red","blue"))
gb

#Both MCL and CLL samples in one plot
MCL.output <- output[,c("chr","Odd.Ratio","remark")]
CLL.output <- output[,c("chr","Odd.Ratio","remark")]
cell.grp <- c(rep("MCL",20),rep("CLL",20))
final_output <- cbind(cell.grp,rbind(MCL.output,CLL.output))
final_output[,"chr"] <- factor(final_output[,"chr"],levels=list.chr)

gb <- ggplot(final_output, aes(x=cell.grp, y=log2(Odd.Ratio), fill=remark))
gb <- gb + xlab("Cell") + ylab("Log2 Odd Ratio") +  ylim(-0.6,0.6)
gb <- gb + geom_bar(stat="identity")
gb <- gb + theme_classic() + theme(strip.background = element_blank(),strip.text.x = element_blank()) + facet_wrap(~chr, ncol=1, drop=FALSE)+coord_flip()
gb <- gb+scale_fill_manual(values=c("red","blue"))
gb

