#_______________________________________________________________________
#
# t-SNE plot
#
# Generates tSNE plot from matrix of IF values
#_______________________________________________________________________

tsneplot <- function(filename){
  set.seed(13)
 # load("chr11_69_70.Rdata")
  load(paste0(filename,".Rdata"))
  matrix[,"start"] <- matrix[,"start"]*1e6
  matrix[,"end"] <- matrix[,"end"]*1e6
  matrix <- subset(matrix, chrom %in% c(1:10,12:13,15:22) & apply(matrix[,-(1:4)], 1, sd)>0)
  matrix <- log2(matrix[,-(1:4)]+1)
  matrix <- t(na.omit(matrix))
  #mycolours <- c("palegreen","palegreen","palegreen", "red", "red", "red", "blue","blue","blue", "mediumpurple2","mediumpurple2","mediumpurple2", "yellow","yellow","yellow","yellow","yellow","yellow","yellow","orange","orange","orange","orange","orange")
  mycolours <- c(rep("orange",5),rep("yellow",7),rep("palegreen", 3),rep("red",3),rep("blue",3),rep("mediumpurple2",3))
  library("Rtsne")
  tsne_results <- Rtsne(matrix, dims=2,perplexity=5, pca_center=TRUE, pca_scale=TRUE, max_iter = 10000)
  plot(tsne_results$Y, frame=FALSE, bg = mycolours, col="black", lwd=1, pch = 21, cex = 1.5, main=c(filename,"scale=T"))
}

files = c("chr11_68_69", "chr11_69_70", "chr11_70_71",
          "chr14_104_105", "chr14_105_106")

for(i in 1:5){
  tsneplot(files[i])
}


set.seed(13)
load("chr14_106_107.Rdata")
matrix[,"start"] <- matrix[,"start"]*1e6
matrix[,"end"] <- matrix[,"end"]*1e6
matrix <- subset(matrix, chrom %in% c(1:10,12:13,15:22))
matrix <- log2(matrix[,-c(1:4,16)]+1)
idx <- which(rowSums(matrix == log2(1)) == 23)
matrix <- matrix[-idx,]
matrix <- t(na.omit(matrix))
mycolours <- c(rep("orange",5),rep("yellow",6),rep("palegreen", 3),rep("red",3),rep("blue",3),rep("mediumpurple2",3))
library("Rtsne")
tsne_results <- Rtsne(matrix, dims=2,perplexity=5, pca_center=TRUE, pca_scale=TRUE, max_iter = 10000)
plot(tsne_results$Y, frame=FALSE, bg = mycolours, col="black", lwd=1, pch = 21, cex = 1.5, main=c("chr14_106_107","scale=T"))

