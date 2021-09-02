# ALLUVIAL
#
# uses the packages: ggalluviall, ggplot2, data.table
#
# updated as of: 2021/09/01
#
#
#---------------------Create the matrix_100kb-----------------------------------------
setwd("~/Google Drive/Data/CompartmentScores")
list.samples = c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828",
                 "mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182",
                 "NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3",
                 "PBC_1","PBC_2","PBC_3")
list.chr = c(1:22)

tmp = read.table(paste0("chr.",list.chr[1],".compartment.data.scaled.txt"), header=T)
matrix_100kb = data.frame(chrom = rep(list.chr[1],nrow(tmp)),tmp)

for(i in 2:22){
  tmp = read.table(paste0("chr.",list.chr[i],".compartment.data.scaled.txt"), header=T)
  matrix_100kb_tmp = data.frame(chrom= rep(list.chr[i],nrow(tmp)),tmp)
  matrix_100kb = rbind(matrix_100kb, matrix_100kb_tmp, stringsAsFactors = FALSE)
}

#Matrix of CS for every sample in every 100kb bin; lacks bin number! 
matrix_100kb = matrix_100kb[,c("chrom",list.samples)]

#incorporate bin number
current = subset(matrix_100kb, chrom == 1)
tmp = data.frame(bin = 0:(nrow(current)-1))
tmp$start = tmp$bin*1e5
tmp$end = tmp$start + 1e5

for(i in 2:22){
  current = subset(matrix_100kb, chrom == i)
  tmp2 = data.frame(bin = 0:(nrow(current)-1))
  tmp2$start = tmp2$bin*1e5
  tmp2$end = tmp2$start + 1e5
  tmp = rbind(tmp, tmp2, stringsAsFactors = FALSE)
}

#Matrix of CS for every sample in 100kb windows
matrix_100kb = data.frame(tmp, matrix_100kb)

rm(tmp2, tmp, matrix_100kb_tmp, current)

#---------------------Subset Matrix for 1Mb--------------------------------

chrom.lengths <- ceiling(table(matrix_100kb$chrom)/10)

SubsetMatrix = function(chr){
  SubsetMatrix_100Kb = subset(matrix_100kb, chrom == chr)
  chromosome = 0
  SubsetMatrix_1Mb = c()
  ListOfTen= subset(SubsetMatrix_100Kb, start %in% seq(chromosome*1e6, (chromosome+1)*1e6-1e5, by=1e5))
  AveListOfTen = apply(ListOfTen[,-(1:4)],2, median, na.rm=T)
  SubsetMatrix_1Mb = data.frame(bin = 0, start = chromosome*1e6, end = (chromosome+1)*1e6, chrom = chr, t(AveListOfTen))
  
  for(chromosome in 1:(chrom.lengths[chr]-1)){
    ListOfTen= subset(SubsetMatrix_100Kb, start %in% seq(chromosome*1e6, (chromosome+1)*1e6-1e5, by=1e5))
    AveListOfTen = apply(ListOfTen[,-(1:4)],2, median, na.rm=T)
    AveListOfTen = data.frame(bin = chromosome, start = chromosome*1e6, end = (chromosome+1)*1e6, chrom = chr, t(AveListOfTen))
    SubsetMatrix_1Mb = rbind(SubsetMatrix_1Mb, AveListOfTen, stringsAsFactors =FALSE)
  }
  SubsetMatrix_1Mb <<- SubsetMatrix_1Mb
}

curr.matrix = list()
for(i in 1:22){
  curr.matrix[[i]] <- SubsetMatrix(i)
}

#Matrix 100kb CS converted to 1Mb
matrix_1Mb = do.call("rbind", curr.matrix)

rm(SubsetMatrix_1Mb, curr.matrix)

#---------------------Data Imputation, CS analysis, and Integrate IF data------------------------
#Data Imputation = setting which is experimental between MCL and CLL
#Analysis = labeling of Compartment State
#Integrate = IF data to data frame

#Step 1: Extract the data of IF
setwd("~/Google Drive/Data/Interaction Scores (Genome Wide)")
Matrix_Int_CS = function(cell_group){
  IntScores = list()
  for(i in 1:22){
    IntScores[[i]] = read.table(paste0(cell_group,".chr.",i,".int.diff.short.minus.long.txt"))
  }
  IntScores <<- IntScores
}
IntScores_MCL = do.call(rbind,Matrix_Int_CS("MCL"))
IntScores_CLL = do.call(rbind,Matrix_Int_CS("CLL"))

#Step 2: Combine Matrix 1Mb with 3) IF      
tmp = data.frame(matrix_1Mb, Int_MCL = unlist(IntScores_MCL), Int_CLL = unlist(IntScores_CLL))

#Step 3: 1) data imputation 2) labelling of compartment state
DataImputation = function(cell_group, matrix){
  #cell_group ...... either MCL or CLL
  #matrix...........either matrix_100kb or matrix_1Mb
  
  ifelse(cell_group == "MCL", 
         experimental <<- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
         experimental <<- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
  
  control <<- c("NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3","PBC_1","PBC_2","PBC_3")
  matrix <- matrix
}
AnalysisMatrix = function(cell_group, SelectedMatrix){
  matrix = DataImputation(cell_group,SelectedMatrix)
  
  #Difference of PC1 values
  ave_exp = apply(matrix, 1, function(x) median(x[experimental]))
  ave_control = apply(matrix, 1, function(x) median(x[control]))
  new_matrix = data.frame(matrix,ave_exp,ave_control)
  
  new_matrix$diff = new_matrix$ave_exp - new_matrix$ave_control
  
  #Compartment states label
  new_matrix$cs_status = with(new_matrix, ifelse(ave_control > 0 & ave_exp > 0, "A-A",
                                                 ifelse(ave_control>0 & ave_exp < 0, "A-B",
                                                        ifelse(ave_control < 0 & ave_exp < 0,"B-B",
                                                               ifelse(ave_control <  0 & ave_exp > 0, "B-A", NA)))))
  
  #Define the cell group which was compared to Normals
  new_matrix$cell.grp = cell_group
  new_matrix <<- new_matrix
}
matrix_1Mb_MCL = AnalysisMatrix("MCL",tmp)[,c("bin","start","end","chrom", "Int_MCL","ave_exp","ave_control","diff", "cs_status")]
matrix_1Mb_CLL = AnalysisMatrix("CLL",tmp)[,c("bin","start","end","chrom","Int_CLL","ave_exp","ave_control","diff","cs_status")]

rm(IntScores, IntScores_CLL, IntScores_MCL, new_matrix, tmp) 




#---------------------Merge CS and IF data-----------------------------------
#
#.....generates a data.frame of only the significantly increased/decreased bins
#
#------------EDIT STARTING HERE
#.....sample - values are any of the following "MCL", "CLL", "MCL.not.CLL"
#.....select.matrix - whether it is "matrix_1Mb_MCL" or "matrix_1Mb_CLL"

sample = "MCL.not.CLL" 
select.matrix = matrix_1Mb_MCL

#--------------DO NOT EDIT
setwd("/Users/kimquililan/Google Drive/Data/Alluvial/Files")
listData = list.files()
listData
sample_list = listData[grep(paste0("^",sample,".[ud]"), listData, ignore.case = TRUE)]
IF_data = rbind.data.frame(read.table(sample_list[1],header=TRUE), read.table(sample_list[2],header=TRUE))
CS_data= apply(IF_data, 1, function(x) select.matrix[which(select.matrix$start == x["start"]*1e6 & select.matrix$chrom == x["chr"]),c("ave_exp","ave_control","diff", "cs_status")])
library(data.table)
CS_data = rbindlist(CS_data)
colnames(CS_data) = c("CS.Normal", paste0("CS.",sample), "CS.change", "CS.status")
alluvial_data = data.frame(IF_data, CS_data, sample = rep(sample,nrow(IF_data)))
alluvial_data$IF =  with(alluvial_data, ifelse(log2FC > 0, "inc", "dec"))
alluvial_data[, paste0("CS.status.",sample)] = with(alluvial_data, ifelse(paste0("CS.",sample) > 0, "A", "B")) 
alluvial_data$chr.2 = with(alluvial_data, ifelse(chr.1 == 11, 
                                                 ifelse(start.1 > 0 & start.1 < 69, "11.1", "11.2"), chr.1))
alluvial_data$chr.2 = factor(alluvial_data$chr.2, levels = c(1:10,11.1,11.2,12:22))



#---------------------Create Alluvial Plots-------------------------------------
#
#....generates a data.frame and plot of interactions from each bin of each chromosome to every other chromosome
#....For example, in a 1Mb window, chromosome 1 has 249 bins, so the data.frame for chr1 is 249 x 23 with
#....a maximum of 2 lines (inc/dec) arising from it. ("stable" white, "inc" red3, "dec" mediumblue).
#....It also colors every bin according to its compartment score change 
#....("A-A" indianred, "A-B" turquoise1, "B-A" hotpink, "B-B" steelblue)
#

library(ggplot2)
library(ggalluvial)

df = select.matrix[,c("bin","chrom","start","end", "cs_status")]
df$bin = 0:(nrow(df)-1)

df = df[rep(seq_len(nrow(df)), each = 23), ]
df$chr.1 = rep(c(1:10,11.1,11.2,12:22),2887)
df$IF = rep("stable",66401)

for (i in 1:nrow(alluvial_data)){
  index = which(df$bin == alluvial_data$bin[i] & df$chr.1 == alluvial_data$chr.2[i])
  df$IF[index] = alluvial_data$IF[i]   
}
table(df$IF)
table(alluvial_data$IF)


setwd("/Users/kimquililan/Google Drive/Data/Alluvial")
chrom_list <- read.table("chrom.bins.1Mb.txt")

setwd(paste0("/Users/kimquililan/Google Drive/Data/Alluvial/Plots/",sample))

#for all other chromosomes
for(i in c(1:10,12:22)){
  chr = i
  start = 0e+06
  end = chrom_list[i,2]*1e06
  
  current = df[which(df$chrom == chr & (df$start >= start & df$start < end)),]
  current$chr.2 = format(round(current$chr.1,1),nsmall = 1)
  
  pdf(paste0(sample,"_chr_",chr,"_cs.pdf"))
  figure_cs <- ggplot(data = current,
                      aes(axis1 = bin, axis2 = chr.2, fill=cs_status)) +
    scale_x_discrete(limits = c("Bin", "To Chromosome"), expand = c(.2, .05)) +
    geom_stratum(width = 1/8, size=0.1) + 
    scale_fill_discrete(type=c("indianred", "turquoise1","hotpink","steelblue"), na.value="gray80") + 
    theme_void() +
    ggtitle(paste0("from chr",chr,": ", start," - ",end))
  print(figure_cs)
  dev.off()
  
  current$IF = factor(current$IF, levels=c("stable","dec","inc"))
  pdf(paste0(sample,"_chr_",chr,"_inc.pdf"))
  figure_inc <- ggplot(data = current,
                       aes(axis1 = bin, axis2 = chr.2)) +
    scale_x_discrete(limits = c("Bin", "To Chromosome"), expand = c(.2, .05)) +
    geom_alluvium(aes(fill = IF), width=1/12,alpha=1) +
    geom_stratum(width = 1/8, size=0.1) + 
    scale_fill_discrete(type=c("white","white","red3")) + 
    theme_void() +
    ggtitle(paste0("from chr",chr,": ", start," - ",end))
  print(figure_inc)
  dev.off()
  
  current$IF = factor(current$IF, levels=c("stable","inc","dec"))
  pdf(paste0(sample,"_chr_",chr,"_dec.pdf"))
  figure_dec <- ggplot(data = current,
                       aes(axis1 = bin, axis2 = chr.2)) +
    scale_x_discrete(limits = c("Bin", "To Chromosome"), expand = c(.2, .05)) +
    geom_alluvium(aes(fill = IF), width=1/12, alpha=1) +
    geom_stratum(width = 1/8, size=0.1) + 
    scale_fill_discrete(type=c("white","white","mediumblue")) + 
    theme_void() +
    ggtitle(paste0("from chr",chr,": ", start," - ",end))
  print(figure_dec)
  dev.off()
}




#for chromosome 11.1 and 11.2
chr = 11
start = 69e+06 
end = 136e+06

current = df[which(df$chrom == chr & (df$start >= start & df$start < end)),]
current$chr.2 = format(round(current$chr.1,1),nsmall = 1)

pdf(paste0(sample,"_chr_11.2_cs.pdf"))
figure_cs <- ggplot(data = current,
                    aes(axis1 = bin, axis2 = chr.2, fill=cs_status)) +
  scale_x_discrete(limits = c("Bin", "To Chromosome"), expand = c(.2, .05)) +
  geom_stratum(width = 1/8, size=0.1) + 
  scale_fill_discrete(type=c("indianred","turquoise1","hotpink","steelblue"), na.value="gray80") + 
  theme_void() +
  ggtitle(paste0("from chr",chr,": ", start," - ",end))
print(figure_cs)
dev.off()

current$IF = factor(current$IF, levels=c("stable","dec","inc"))
pdf(paste0(sample,"_chr_11.2_inc.pdf"))
figure_inc <- ggplot(data = current,
                     aes(axis1 = bin, axis2 = chr.2)) +
  scale_x_discrete(limits = c("Bin", "To Chromosome"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = IF), width=1/12, alpha=1) +
  geom_stratum(width = 1/8, size=0.1) + 
  scale_fill_discrete(type=c("white","white","red3")) + 
  theme_void() +
  ggtitle(paste0("from chr",chr,": ", start," - ",end))
print(figure_inc)
dev.off()

current$IF = factor(current$IF, levels=c("stable","inc","dec"))
pdf(paste0(sample,"_chr_11.2_dec.pdf"))
figure_dec <- ggplot(data = current,
                     aes(axis1 = bin, axis2 = chr.2)) +
  scale_x_discrete(limits = c("Bin", "To Chromosome"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = IF), width=1/12, alpha=1) +
  geom_stratum(width = 1/8, size=0.1) + 
  scale_fill_discrete(type=c("white","white","mediumblue")) + 
  theme_void() +
  ggtitle(paste0("from chr",chr,": ", start," - ",end))
print(figure_dec)
dev.off()

