#_______________________________________________________________________
#
# Length of Chromosome vs % Increased Bins
#
# Generates a scatter plot of %increased bins in each chromosome
# and also a plot comparing CLL and MCL
#_______________________________________________________________________
#
#
# Data needed for analysis: 
#   output: matrix of Odds ratio values of all bins; must have column "chr" and "Odd Ratio"

circos.matrix <- matrix(0,nrow=20,ncol=5)
colnames(circos.matrix) = c("chr", "chr_length", "inc", "dec", "remark")

list.chr = c(1:10,12:13,15:22)
chrom_length_orig <- read.table("chrom.bins.1Mb.txt")

for (i in 1:20){
  current.chrom <- output[which(output[,"chr"]==list.chr[i]),]
  increased = length(which(current.chrom[,"Odd Ratio"]>1))
  length = chrom_length_orig[which(chrom_length_orig[,1]==list.chr[i]),2]
  fraction = (increased/dim(current.chrom)[1])*100
  circos.matrix[i,"chr"] = list.chr[i]
  circos.matrix[i,"chr_length"] = length
  circos.matrix[i,"inc"] = fraction
  circos.matrix[i,"dec"] = 100-fraction
  circos.matrix[i,"remark"] = ifelse (circos.matrix[i,"inc"] > 50, "increased","decreased")
}

circos.matrix <- data.frame(circos.matrix)
circos.matrix[,"chr_length"] = as.numeric(as.character(circos.matrix[,"chr_length"]))
circos.matrix[,"inc"] = as.numeric(as.character(circos.matrix[,"inc"]))
circos.matrix[,"dec"] = as.numeric(as.character(circos.matrix[,"dec"]))
circos.matrix[,"remark"] = factor(circos.matrix[,"remark"], levels=c("increased","decreased"))

#Basic Plot
plot(circos.matrix[,"chr_length"], circos.matrix[,"inc"], ylim=c(0,100), cex=1.25,pch=21, col="black", bg="orange", xlab="chromosome length", ylab="% increased", main="MCL vs normal")

#ggplot
library("ggplot2")
ggplot(circos.matrix, aes(x=chr_length, y=inc, colour=remark)) + 
  geom_point(size=1.5) + scale_colour_manual(values=c("red", "blue")) + ylim(0,100) + 
  theme_classic() + geom_hline(yintercept=50, lty=2) + geom_vline(xintercept=108, lty=2) +
  ylab("Percentage of Increased Bins") + xlab("Length of Chromosome") + theme(legend.position = "none")

#Comparing MCL and CLL
#Note: must already have a circos.matrix for MCL vs normals, and circos.matrix for CLL vs normals
cell.grp = c(rep("MCL",20),rep("CLL",20))
MCL.data = circos.matrix.MCL[,c("chr","chr_length","inc")]
CLL.data = circos.matrix.CLL[,c("chr","chr_length","inc")]
MCL.CLL.data = data.frame(cell.grp,rbind(MCL.data, CLL.data))

library("ggplot2")
large.chrom = c(1:12,21:32)
small.chrom = c(13:20,33:40)
ggplot(MCL.CLL.data[large.chrom,],
       aes(x=cell.grp, y=inc, fill=cell.grp)) + theme_classic() + 
       geom_point(shape=21, size=2.75,) + geom_line(aes(group=chr)) +
       scale_fill_manual(values=c("yellow2","orange")) +
       ylab("Percentage of Increased Bins") + theme(legend.position = "none")
