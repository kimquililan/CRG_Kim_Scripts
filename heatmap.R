#_______________________________________________________________________
#
# Heatmap
#
# Generates heatmap using a Euclidean method of distance calculation
# and a clustering method of ward.D2
# 
#_______________________________________________________________________

#Needed for analysis:
#   output = any matrix with samples as columns and bins as rows, and samples 
#            are ordered by control then experimental


library(gplots)
hmcol = colorRampPalette(c("royalblue4","dodgerblue4","dodgerblue3","steelblue3","lightblue3","lavender","white","papayawhip","wheat2","sandybrown","brown2","brown3","brown4"))(100)
color.breaks=seq(-1.5,1.5,length=101)
sample.colors<-c(rep("palegreen",3), 
                 rep("red",3), 
                 rep("blue",3),
                 rep("mediumpurple2",3),
                 #rep("yellow2",7),
                 rep("orange",5))
                 
heatmap.2(output, distfun = function(x) {dist(x, method="euclidean")},
          hclustfun = function(x) {hclust(x, method="ward.D2")}, trace="none",density.info="none", col = hmcol, breaks=color.breaks, 
           labRow = FALSE, labCol = FALSE, ColSideColors=sample.colors)



