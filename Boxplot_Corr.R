experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828")
experimental2 <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182")

control = c("NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3","PBC_1","PBC_2","PBC_3")


BoxplotMatrix = function(filename){
  load(paste0(filename,".Rdata"))
  list.chr = c(1:10,12:13,15:22)
  matrix.sum = aggregate(matrix[,c(control,experimental,experimental2)], by=list(matrix$chrom), sum, na.rm=TRUE)
  matrix.sum = matrix.sum[which(matrix.sum$Group.1 %in% list.chr),]
  matrix.sum[,1] = as.numeric(as.character(matrix.sum[,1]))
  matrix.sum = matrix.sum[order(matrix.sum$Group.1),]
  
  matrix.fraction = matrix(0,nrow=20,ncol=24)
  rownames(matrix.fraction)= list.chr
  for(j in 1:24){
    for(i in 1:20){
      matrix.fraction[i,j] = matrix.sum[i,j+1]/sum(matrix.sum[,j+1])
    }
  }
  
  colnames(matrix.fraction) = c(rep("Normals",12), rep("MCL",5), rep("CLL",7))
  matrix.fraction= t(matrix.fraction)
  
  output = data.frame(chrom=rep("1",24), cell.grp=rownames(matrix.fraction), viewpoint=rep(filename,24),value=matrix.fraction[,1])
  
  for(i in 2:20){
    tmp = data.frame(chrom=rep(list.chr[i],24), cell.grp=rownames(matrix.fraction), viewpoint=rep(filename,24),value=matrix.fraction[,i])
    output = rbind.data.frame(output,tmp)
    rm(tmp)
  }
  
  return(assign(paste0("output_",filename), output, env=.GlobalEnv))
  
}

BoxplotMatrix("chr11_68_69")
BoxplotMatrix("chr11_69_70")
BoxplotMatrix("chr11_70_71")
BoxplotMatrix("chr14_104_105")
BoxplotMatrix("chr14_105_106")
BoxplotMatrix("chr14_106_107")

output <- rbind.data.frame(output_chr11_68_69,output_chr11_69_70,output_chr11_70_71,output_chr14_104_105,output_chr14_105_106,output_chr14_106_107)


library(ggplot2)
library(gridExtra)

Boxplot_Plot= function(chr){
p1= ggplot(output[which(output$chrom == chr & output$cell.grp == "Normals"),], aes(x=viewpoint,y=value,fill=cell.grp)) + geom_boxplot() + 
  scale_fill_manual(values="grey") + xlab("Normals") + ylim(0,0.15) + theme(panel.background = element_rect(fill="white"),
                                           axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
                                            panel.border = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(), legend.position = "none",
                                           axis.text.x = element_text(angle = 90))

p2 = ggplot(output[which(output$chrom == chr & output$cell.grp == "MCL"),], aes(x=viewpoint,y=value,fill=cell.grp)) + geom_boxplot() + 
  scale_fill_manual(values="orange") + ylim(0,0.15) + xlab("MCL") + theme(panel.background = element_rect(fill="white"),
                                           axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
                                           panel.border = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(), legend.position = "none",
                                           axis.title.y=element_blank(),
                                           axis.text.y=element_blank(),
                                           axis.ticks.y=element_blank(),
                                           axis.line.y = element_blank(),axis.text.x = element_text(angle = 90))

p3 = ggplot(output[which(output$chrom == chr & output$cell.grp == "CLL"),], aes(x=viewpoint,y=value,fill=cell.grp)) + geom_boxplot() + 
  scale_fill_manual(values="yellow2") + ylim(0,0.15) + xlab("CLL") + theme(panel.background = element_rect(fill="white"),
                                           axis.line = element_line(size = 0.5, linetype = "solid",colour = "black"),
                                           panel.border = element_blank(),
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(), legend.position = "none",
                                           axis.title.y=element_blank(),
                                           axis.text.y=element_blank(),
                                           axis.ticks.y=element_blank(),
                                           axis.line.y = element_blank(),axis.text.x = element_text(angle = 90))

png(paste0("chr",chr,"_boxplot.png"))
grid.arrange(p1,p2,p3,nrow=1,top=textGrob(paste0("chr",chr),gp=gpar(fontsize=15,font=3)) )
dev.off()
}
