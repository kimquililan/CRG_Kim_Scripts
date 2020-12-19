# ------------------------------------------------------------------------------------
# 
# COMPARTMENT SCORES
#
# Compares compartment scores between MCL and normals / CLL and normals
#
#-------------------------------------------------------------------------------------

#Create the matrix
list.samples = c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828",
                 "mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182",
                 "NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3",
                 "PBC_1","PBC_2","PBC_3")
list.chr = c(1:22)

tmp = read.table(paste0("chr.",list.chr[1],".compartment.data.scaled.txt"), header=T)
matrix = data.frame(chrom = rep(list.chr[1],nrow(tmp)),tmp)

for(i in 2:22){
  tmp = read.table(paste0("chr.",list.chr[i],".compartment.data.scaled.txt"), header=T)
  matrix2 = data.frame(chrom= rep(list.chr[i],nrow(tmp)),tmp)
  matrix = rbind(matrix, matrix2, stringsAsFactors = FALSE)
}

matrix = matrix[,c("chrom",list.samples)]

current = subset(matrix, chrom == 1)
tmp = data.frame(bin = 0:(nrow(current)-1))
tmp$start = tmp$bin*1e5
tmp$end = tmp$start + 1e5

for(i in 2:22){
  current = subset(matrix, chrom == i)
  tmp2 = data.frame(bin = 0:(nrow(current)-1))
  tmp2$start = tmp2$bin*1e5
  tmp2$end = tmp2$start + 1e5
  tmp = rbind(tmp, tmp2, stringsAsFactors = FALSE)
}

matrix = data.frame(tmp, matrix)

#---------------------------Data Imputation------------------------#

cell_group = "CLL"
cell_group = "MCL"
ifelse(cell_group == "MCL", 
       experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
       experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))

control = c("NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3","PBC_1","PBC_2","PBC_3")

NAs.exp<-rowSums(is.na(matrix[,experimental]))
NAs.control<-rowSums(is.na(matrix[,control]))
rows.to.include<-which(NAs.exp<=(length(experimental)-3) & NAs.control<=(length(control)-3))
matrix<-matrix[rows.to.include,]
dim(matrix)

#----------------------------Analysis-----------------------------
#Difference of PC1 values
mean_exp = apply(matrix, 1, function(x) mean(x[experimental]))
mean_control = apply(matrix, 1, function(x) mean(x[control]))
new_matrix = data.frame(matrix,mean_exp,mean_control)

new_matrix$diff = new_matrix$mean_exp - new_matrix$mean_control

#Compartment states label
new_matrix$cs_status = with(new_matrix, ifelse(mean_control > 0 & mean_exp > 0, "A-A",
                        ifelse(mean_control>0 & mean_exp < 0, "A-B",
                               ifelse(mean_control < 0 & mean_exp < 0,"B-B",
                                      ifelse(mean_control <  0 & mean_exp > 0, "B-A", NA)))))

#Define the cell group which was compared to Normals
new_matrix$cell.grp = cell_group

                     
#---------------------------Visualization------------------------

#Genome wide              
  #Barplots
    par(mfrow=c(2,1))
    for(i in 1:22){
      barplot(new_matrix[which(new_matrix$chrom == i),"mean_control"],main=paste0("chr",i),ylim=c(-1,1),ylab="Normal",col=ifelse(new_matrix[which(new_matrix$chrom == i),"mean_control"]>0,"red","blue"),border=NA)
      barplot(new_matrix[which(new_matrix$chrom == i),"mean_exp"],ylim=c(-1,1),ylab="MCL",col=ifelse(new_matrix[which(new_matrix$chrom == i),"mean_exp"]>0,"red","blue"),border=NA)
    }
    for(i in 1:22){
      barplot(new_matrix_MCL[which(new_matrix_MCL$chrom == i),"diff"], main=paste0("chr",i), ylim=c(-1,1), ylab="Mean PC1 values (MCL-normal)", col=ifelse(new_matrix_MCL[which(new_matrix_MCL$chrom == i),"diff"]>0,"red","blue"),border=NA)
    }

  #Boxplots
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.border = element_rect(colour = "black", fill=NA, size=1))
    data = rbind(new_matrix_CLL,new_matrix_MCL)
    data$chrom = factor(data$chrom)
    ggplot(data, aes(x=chrom, y=diff, fill=cell.grp)) + geom_boxplot() + facet_wrap(~chrom, scale="free") + scale_fill_manual(values=c("yellow","orange")) + geom_hline(yintercept=0,lty=2) + ggtitle("All bins")
    ggplot(data, aes(x=chrom, y=diff, fill=cell.grp)) + geom_boxplot() + scale_fill_manual(values=c("yellow","orange")) + geom_hline(yintercept=0,lty=2) + ggtitle("All bins") + theme_minimal()
    ggplot(new_matrix_CLL, aes(y=diff, x=chrom)) + geom_boxplot(fill="yellow") + geom_hline(yintercept=0,lty=2) + ggtitle("All bins - CLL vs normals") 
    ggplot(new_matrix_MCL, aes(y=diff, x=chrom)) + geom_boxplot(fill="orange") + geom_hline(yintercept=0,lty=2) + ggtitle("All bins - MCL vs normals")

    boxplot(diff~chrom, data=new_matrix_CLL,pch=19,ylim=c(-1,1),col="yellow",main="CLL vs normals")
    boxplot(diff~chrom, data=new_matrix_MCL,pch=19,ylim=c(-1,1),col="orange",main="MCL vs normals")

  #Histogram
    #Difference only
      par(mfrow=c(3,3))
      i=0              
      for(i in 1:22){
        d.2 = density(subset(new_matrix_MCL, chrom == i)$diff, na.rm=T)
        plot(d.2,type="n", main=paste0("chrom",i))
        lines(d.2,col="black")
        abline(v=0, lty=2, col="grey")
        mean_density = mean(subset(new_matrix_MCL, chrom == i)$diff,na.rm=T)
        abline(v=mean_density, lty=2, col="grey4")
      }

    #5 MCls and 12 normals
      par(mfrow=c(3,3))
      i=0              
      for(i in 1:22){

        exp_hist = apply(subset(new_matrix_CLL, chrom == i)[experimental], 2, function(x) density(x, na.rm=T))
        ctr_hist = apply(subset(new_matrix_CLL, chrom == i)[control], 2, function(x) density(x, na.rm=T))
        plot(exp_hist[[1]], type="n", xlim=c(-2,2), main=paste0("chr",i))
        for(j in 1:length(experimental)){ lines(exp_hist[[j]], col="yellow4")}
        for(k in 1:length(control)){ lines(ctr_hist[[k]], col="grey")}
      }
  
    #Mean MCL and Mean Normal 
      i=0                  
      for(i in 1:22){
        exp = density(subset(new_matrix_CLL, chrom == i)$mean_exp, na.rm=T)
        ctr = density(subset(new_matrix_CLL, chrom == i)$mean_control, na.rm=T)
        plot(exp, xlim=c(-2,2), main=paste0("chr",i), type="n")
        lines(exp, col="yellow4")
        lines(ctr, col="grey")
        }


                   
#Grouped by Compartment State
  #Pie Chart
  pie_chart = data.frame(table(new_matrix[,"cs_status"]))
  pie_chart$percentage = pie_chart[,2]/sum(pie_chart[,2])*100
  pie(pie_chart$percentage, labels = NA, col=c("grey","blue","red","black"))

  perc_chr15_22 = c()
  perc_chr15_22[1] = length(which(new_matrix$cs_status == "A-A" & new_matrix$chrom %in% 15:22))
  perc_chr15_22[2] = length(which(new_matrix$cs_status == "A-B" & new_matrix$chrom %in% 15:22))
  perc_chr15_22[3] = length(which(new_matrix$cs_status == "B-A" & new_matrix$chrom %in% 15:22))
  perc_chr15_22[4] = length(which(new_matrix$cs_status == "B-B" & new_matrix$chrom %in% 15:22))
  pie_chart$perc_chr15_22 = perc_chr15_22/pie_chart$Freq*100
  pie_chart$perc_rest = 100-pie_chart$perc_chr15_22
  counts = t(pie_chart[,c("perc_chr15_22","perc_rest")])
  barplot(counts, col=c("blue","grey"), legend=rownames(counts),legend.position = "topright",ylab="Percentage of Compartments")

  #Boxplot
  boxplot(diff~chrom, data=new_matrix[which(new_matrix$cs_status == "A-B"),],pch=19)
  boxplot(diff~chrom, data=new_matrix[which(new_matrix$cs_status == "B-A"),],pch=19)
  boxplot(diff~chrom, data=new_matrix[which(new_matrix$cs_status == "B-B"),],pch=19)
  boxplot(diff~chrom, data=new_matrix[which(new_matrix$cs_status == "A-B"),],pch=19,ylim=c(-1,1))
  boxplot(diff~chrom, data=new_matrix[which(new_matrix$cs_status == "B-A"),],pch=19,ylim=c(-1,1))
  boxplot(diff~chrom, data=new_matrix[which(new_matrix$cs_status %in% c("A-B","B-A")),],pch=19,ylim=c(-1,1),col="yellow",main="CLL vs normals")
  abline(h=0,lty=2)


