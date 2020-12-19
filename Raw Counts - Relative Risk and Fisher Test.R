#_______________________________________________________________________
#
# Raw counts - Relative Risk and Fisher Test
#
# Calculates Relative Risk (RR) for raw counts and performs Fisher Test to calculate
# significance. Additional: Creates volcano plots
# Requires new_matrix from "Create_Corr_Matrix" function
#_______________________________________________________________________

#Relative_Risk("chr14_105_106", "MCL")

Relative_Risk = function(viewpt_file, cell_group){

  setwd("~/Google Drive/Data/CoverageCorrected")
  load(paste0(viewpt_file,".Rdata"))
  matrix[,"start"] <- matrix[,"start"]*1e6
  matrix[,"end"] <- matrix[,"end"]*1e6
  matrix <- subset(matrix, chrom %in% c(1:11,12:14,15:22))
  
  
  #Choose between MCL or CLL as experimental
  ifelse(cell_group == "MCL", 
         experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
         experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
  
  control = c("NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3","PBC_1","PBC_2","PBC_3")
  
  #filter rows based on number of NAs
  NAs.exp<-rowSums(is.na(matrix[,experimental]))
  NAs.control<-rowSums(is.na(matrix[,control]))
  rows.to.include<-which(NAs.exp<=(length(experimental)-3) & NAs.control<=(length(control)-3))
  matrix<-matrix[rows.to.include,]
  
  
  #Calculation 
  output <- matrix[,c("bin","chrom","start", "end")]
  output$RR = 0 
  output$p.value = 0
  
  colSum = data.frame(apply(matrix[,c(experimental,control)], 2, sum, na.rm=T))
  
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
    RR =   (a*d)/(b*c)
    output[i,"RR"] <- RR
    output[i,"p.value"] <- fisher.test(matrix(c(a,b,c,d),nrow=2, ncol=2))$p.value
    rm(a,b,c,d,RR,tmp)
  }
  
  
  padjust <- p.adjust(output[,"p.value"], method="BH")
  output$p.adjust <- padjust

  output <<- output
}

#-----------------------Results checking------------------------
# Expected Results:
# MCL vs normal: 119 significant bins
# CLL vs normal: 205 significant bins

#Significant Bins
length(which(output[,"p.adjust"]<0.1))
length(which(output[,"RR"]>1 & output[,"p.adjust"]<0.1))
length(which(output[,"RR"]<1 & output[,"p.adjust"]<0.1))


#-------Chromosome Length vs Percentage Increased Bins----------
#Percentage_Inc_Bins(output, "MCL")
#output dependent on Relative_Risk function results

Percentage_Inc_Bins = function(output, cell_group){
  circos.matrix <- matrix(0,nrow=20,ncol=5)
  colnames(circos.matrix) = c("chr", "chr_length", "inc", "dec", "remark")
  
  list.chr = c(1:10,12:13,15:22)
  chrom_length_orig <- read.table("chrom.bins.1Mb.txt")

  for (i in 1:20){
    current.chrom <- output[which(output[,"chrom"]==list.chr[i]),]
    increased = length(which(current.chrom[,"RR"]>1))
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

 
  #ggplot
  library("ggplot2")
  dotplot = ggplot(circos.matrix, aes(x=chr_length, y=inc, colour=remark)) + 
    geom_point(size=1.5) + scale_colour_manual(values=c("red", "blue")) + ylim(0,100) + 
    theme_classic() + geom_hline(yintercept=50, lty=2) + geom_vline(xintercept=108, lty=2) + ggtitle(paste0(cell_group," vs normal B cells"))+
    ylab("Percentage of Increased Bins") + xlab("Length of Chromosome") + theme(legend.position = "none")
  
  print(dotplot)
  ifelse(cell_group == "MCL", circos.matrix.MCL <<- circos.matrix, circos.matrix.CLL <<- circos.matrix)
}

#---------Percentage_Inc_Bins_11(output, "MCL")
Percentage_Inc_Bins_11 = function(output, cell_group){
  circos.matrix <- matrix(0,nrow=20,ncol=5)
  colnames(circos.matrix) = c("chr", "chr_length", "inc", "dec", "remark")
  
  list.chr = c(1:10,12:13,15:22)
  chrom_length_orig <- read.table("chrom.bins.1Mb.txt")
  
  for (i in 1:20){
    current.chrom <- output[which(output[,"chrom"]==list.chr[i]),]
    increased = length(which(current.chrom[,"RR"]>1))
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
  
  
  #ggplot
  library("ggplot2")
  dotplot = ggplot(circos.matrix, aes(x=chr_length, y=inc, colour=remark)) + 
    geom_point(size=1.5) + scale_colour_manual(values=c("red", "blue")) + ylim(0,100) + 
    theme_classic() + geom_hline(yintercept=50, lty=2) + geom_vline(xintercept=136, lty=2) + ggtitle(paste0(cell_group," vs normal B cells"))+
    ylab("Percentage of Increased Bins") + xlab("Length of Chromosome") + theme(legend.position = "none")
  
  print(dotplot)
  ifelse(cell_group == "MCL", circos.matrix.MCL <<- circos.matrix, circos.matrix.CLL <<- circos.matrix)
}

#Comparing MCL and CLL
#Inc_Bins_Comparison(circos.matrix.MCL, circos.matrix.CLL, "large")
Inc_Bins_Comparison = function(circos.matrix.MCL, circos.matrix.CLL, size){
  cell.grp = c(rep("MCL",20),rep("CLL",20))
  MCL.data = circos.matrix.MCL[,c("chr","chr_length","inc")]
  CLL.data = circos.matrix.CLL[,c("chr","chr_length","inc")]
  MCL.CLL.data = data.frame(cell.grp,rbind(MCL.data, CLL.data))
  kim = data.frame(CLL.data, MCL.data[,3])
  remark = c()
for(i in 1:20){
remark[i] = ifelse(kim[i,4]>kim[i,3], "increased", "decreased")
}
  kim = data.frame(kim,remark)
  colnames(kim) = c("chr","chr_length","CLL","MCL","remark")
  print(kim)
  library("ggplot2")
  large.chrom = c(1:12,21:32)
  small.chrom = c(13:20,33:40)
  ifelse(size == "large", chrom.size <- large.chrom, chrom.size <- small.chrom)
  ggplot(MCL.CLL.data[chrom.size,],
         aes(x=cell.grp, y=inc, fill=cell.grp)) + theme_classic() + 
          geom_point(shape=21, size=2.75,) + geom_line(aes(group=chr)) +
          scale_fill_manual(values=c("yellow2","orange")) + ggtitle(paste0(size))+
          ylab("Percentage of Increased Bins") + theme(legend.position = "none")
}

#Inc_Bins_Comparison_11(circos.matrix.MCL, circos.matrix.CLL, "large")
Inc_Bins_Comparison_11 = function(circos.matrix.MCL, circos.matrix.CLL, size){
  cell.grp = c(rep("MCL",20),rep("CLL",20))
  MCL.data = circos.matrix.MCL[,c("chr","chr_length","inc")]
  CLL.data = circos.matrix.CLL[,c("chr","chr_length","inc")]
  MCL.CLL.data = data.frame(cell.grp,rbind(MCL.data, CLL.data))
  kim = data.frame(CLL.data, MCL.data[,3])
  remark = c()
  for(i in 1:20){
    remark[i] = ifelse(kim[i,4]>kim[i,3], "increased", "decreased")
  }
  kim = data.frame(kim,remark)
  colnames(kim) = c("chr","chr_length","CLL","MCL","remark")
  print(kim)
  library("ggplot2")
  large.chrom = c(1:10,21:30)
  small.chrom = c(11:20,31:40)
  ifelse(size == "large", chrom.size <- large.chrom, chrom.size <- small.chrom)
  ggplot(MCL.CLL.data[chrom.size,],
         aes(x=cell.grp, y=inc, fill=cell.grp)) + theme_classic() + 
    geom_point(shape=21, size=2.75,) + geom_line(aes(group=chr)) +
    scale_fill_manual(values=c("yellow2","orange")) + ggtitle(paste0(size))+
    ylab("Percentage of Increased Bins") + theme(legend.position = "none")
}

#---------ColSum------------------------------
#ColSum_Interaction("chr11_68_69")
ColSum_Interaction <- function(viewpt_file){
  setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/CoverageCorrected")
  load(paste0(viewpt_file,".Rdata"))
  
  list.chr = c(1:10,12:13,15:22)
  intxn.sum.chr.sample= list() 
  for (i in 1:20){
    current.chr = subset(matrix, matrix[,"chrom"]==list.chr[i])
    intxn.sum.chr.sample[[i]] = apply(current.chr[,-c(1:4)],2,sum, na.rm=T)
  }
  
  intxn.sum.chr.sample = data.frame(intxn.sum.chr.sample) 
  colnames(intxn.sum.chr.sample) = list.chr
  intxn.sum.chr.sample = t(intxn.sum.chr.sample)
  intxn.sum.chr.sample <<- intxn.sum.chr.sample
}

#---------Barplot-----------------------------
#Chromosome_vs_RR_Barplot(intxn.sum.chr.sample, "MCL")
Chromosome_vs_RR_Barplot = function(intxn.sum.chr.sample, cell_group){
  ifelse(cell_group == "MCL", 
         experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
         experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
  
  control = c("NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3","PBC_1","PBC_2","PBC_3")
  
  list.chr <- c(1:10,12:13,15:22)
  output <- matrix(0, nrow=20, ncol=5)
  colnames(output) <- c("chr","RR","remark","p-value", "p-adjust")
  output[,"chr"] <- list.chr
  
  for (i in 1:20){
    a = sum(intxn.sum.chr.sample[i,experimental])
    b = sum(intxn.sum.chr.sample[,experimental]) 
    c = sum(intxn.sum.chr.sample[i,control])
    d = sum(intxn.sum.chr.sample[,control]) 
    RR =   (a*d)/(b*c)
    output[i,"RR"] <- RR
    output[i,"remark"] <- ifelse (RR > 1, "increased","decreased")
    tmp.output <- matrix(c(a,b,c,d),nrow=2, ncol=2)
    fishertest.results <- fisher.test(tmp.output)
    output[i,"p-value"] <- fishertest.results$p.value
    rm(a,b,c,d,OR,tmp.output)
  }
  
  padjust <- p.adjust(output[,"p-value"], method="BH")
  output[,"p-adjust"] <- padjust
  output <- data.frame(output)
  output[,"RR"] <- as.numeric(as.character(output[,"RR"]))
  output[,"remark"] <- factor(output[,"remark"], levels=c("increased","decreased"))
  output[,"chr"] <- factor(output[,"chr"],levels=c(22,21,20,19,18,17,16,15,13,12,10,9,8,7,6,5,4,3,2,1))
  
  
  
  library(ggplot2)
  ggplot(output, aes(x= chr, y=log2(RR), fill=remark)) + xlab("Chromosome") + ylab("Log2 Relative Risk") + ylim(-0.5,0.5) + 
    geom_bar(stat="identity") + theme_classic() + coord_flip() + scale_fill_manual(values=c("red","blue")) + ggtitle(paste0(cell_group, " vs normal B cells"))
  
}

#--------Scatterplot------------------------
#Chromosome_vs_RR_Scatterplot(intxn.sum.chr.sample, "MCL")
Chromosome_vs_RR_Scatterplot = function(intxn.sum.chr.sample, cell_group){
  ifelse(cell_group == "MCL", 
         experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
         experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
  
  control = c("NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3","PBC_1","PBC_2","PBC_3")
  
  RR_matrix <- matrix(0, nrow=20, ncol=length(experimental))
  colnames(RR_matrix) <- experimental
  intxn.sum.chr.sample <- intxn.sum.chr.sample+1
  
  for (i in 1:20){
    for (j in 1:length(experimental)){
      a = intxn.sum.chr.sample[i,experimental[j]]
      b = sum(intxn.sum.chr.sample[,experimental[j]]) 
      c = sum(intxn.sum.chr.sample[i,control])
      d = sum(intxn.sum.chr.sample[,control]) 
      RR =   (a*d)/(b*c)
      RR_matrix[i,j] <- RR
      rm(a,b,c,d,RR)
    }
  }
  
  list.chr <- c(1:10,12:13,15:22)
  
  chromosomes <- as.numeric(sapply(list.chr, function(x) rep(x,length(experimental))))
  cell.grp <- rep(cell_group,(20*length(experimental)))
  RR_matrix <- as.numeric(t(RR_matrix))
  output <- data.frame(cell.grp, chromosomes, RR_matrix)
  output[,"chromosomes"] <- factor(output[,"chromosomes"], levels=sort(list.chr, decreasing=T))
  remark <- sapply(output[,"RR_matrix"], function(x) ifelse(x>1, "increased", "decreased"))
  output <- cbind(output, remark)
  output[,"remark"] <- factor(output[,"remark"], levels=c("increased", "decreased"))
  
  library(ggplot2)
  ggplot(output, aes(x= chromosomes, y=log2(RR_matrix), fill=remark)) +
    xlab("Chromosome") + ylab("Log2 Relative Risk") + ylim(-1.5,1.5) + 
    geom_point(shape=21, size=1.75) +
    theme_classic() + coord_flip() +
    scale_fill_manual(values=c("red","blue")) +
    geom_hline(yintercept=0, lty=2) + ggtitle(paste0(cell_group, " vs normal B cells"))
}

#-----------Boxplot-----------------------------







#---------Plotting-------------------
#volcano plot
interchrom <- output
RR <- interchrom[,"RR"]
pval <- interchrom[,"p.adjust"]
vplotdata <- cbind(pval,RR)

plot(log2(vplotdata[,2]), -log10(vplotdata[,1]),
     xlim=c(-3,3),ylim=c(0,10),xlab="Log2 Relative Risk", ylab="-log10(FDR)",
     pch=20,cex=1, col="grey")
points(log2(vplotdata[which(vplotdata[,"RR"]>1 & vplotdata[,"pval"]<0.1),2]),
       -log10(vplotdata[which(vplotdata[,"RR"]>1 & vplotdata[,"pval"]<0.1),1]), col="red", pch=20, cex=1)
points(log2(vplotdata[which(vplotdata[,"RR"]<1 & vplotdata[,"pval"]<0.1),2]),
       -log10(vplotdata[which(vplotdata[,"RR"]<1 & vplotdata[,"pval"]<0.1),1]), col="blue", pch=20, cex=1)
abline(h=1, v=c(-0.1,0.1),lty=2)

