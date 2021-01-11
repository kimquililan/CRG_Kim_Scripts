# ------------------------------------------------------------------------------------
# 
# COMPARTMENT SCORES
#
# Compares compartment scores between MCL and normals / CLL and normals
#
# updated as of: 2021/01/11
#-------------------------------------------------------------------------------------

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

matrix_100kb = matrix_100kb[,c("chrom",list.samples)]

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

matrix_100kb = data.frame(tmp, matrix_100kb)


#--------------------------Subset Matrix for 1Mb--------------------------------

chrom.lengths <- ceiling(table(matrix_100kb$chrom)/10)

SubsetMatrix = function(chr){
  SubsetMatrix_100Kb = subset(matrix_100kb, chrom == chr)
  chromosome = 0
  SubsetMatrix_1Mb = c()
  ListOfTen= subset(SubsetMatrix_100Kb, start %in% seq(chromosome*1e6, (chromosome+1)*1e6-1e5, by=1e5))
  AveListOfTen = apply(ListOfTen[,-(1:4)],2, mean, na.rm=T)
  SubsetMatrix_1Mb = data.frame(bin = 0, start = chromosome*1e6, end = (chromosome+1)*1e6, chrom = chr, t(AveListOfTen))
    
  for(chromosome in 1:(chrom.lengths[chr]-1)){
    ListOfTen= subset(SubsetMatrix_100Kb, start %in% seq(chromosome*1e6, (chromosome+1)*1e6-1e5, by=1e5))
    AveListOfTen = apply(ListOfTen[,-(1:4)],2, mean, na.rm=T)
    AveListOfTen = data.frame(bin = chromosome, start = chromosome*1e6, end = (chromosome+1)*1e6, chrom = chr, t(AveListOfTen))
    SubsetMatrix_1Mb = rbind(SubsetMatrix_1Mb, AveListOfTen, stringsAsFactors =FALSE)
  }
  SubsetMatrix_1Mb <<- SubsetMatrix_1Mb
}

curr.matrix = list()
for(i in 1:22){
  curr.matrix[[i]] <- SubsetMatrix(i)
}

matrix_1Mb = do.call("rbind", curr.matrix)

#-------------------------Data Imputation and Analysis of matrix 100kb------------------------
#Data Imputation = cleaning the data 
#Analysis = labeling of Compartment State

DataImputation = function(cell_group, matrix){
  #cell_group ...... either MCL or CLL
  #matrix...........either matrix_100kb or matrix_1Mb
  
  ifelse(cell_group == "MCL", 
         experimental <<- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
         experimental <<- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
  
  control <<- c("NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3","PBC_1","PBC_2","PBC_3")
  
  NAs.exp<-rowSums(is.na(matrix[,experimental]))
  NAs.control<-rowSums(is.na(matrix[,control]))
  rows.to.include<-which(NAs.exp<=(length(experimental)-3) & NAs.control<=(length(control)-3))
  matrix <- matrix[rows.to.include,]
}

AnalysisMatrix = function(cell_group, SelectedMatrix){
  matrix = DataImputation(cell_group,SelectedMatrix)
  
  #Difference of PC1 values
  ave_exp = apply(matrix, 1, function(x) mean(x[experimental]))
  ave_control = apply(matrix, 1, function(x) mean(x[control]))
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

new_matrix_MCL = AnalysisMatrix("MCL",matrix_100kb)
new_matrix_CLL = AnalysisMatrix("CLL",matrix_100kb)
                     
#---------------------------Visualization------------------------

#Genome wide              
  #Barplots
    par(mfrow=c(2,1))
    new_matrix = new_matrix_MCL
    for(i in 1:22){
      barplot(new_matrix[which(new_matrix$chrom == i),"ave_control"],main=paste0("chr",i),ylim=c(-1,1),ylab="Normal",col=ifelse(new_matrix[which(new_matrix$chrom == i),"ave_control"]>0,"red","blue"),border=NA)
      barplot(new_matrix[which(new_matrix$chrom == i),"ave_exp"],ylim=c(-1,1),ylab="MCL",col=ifelse(new_matrix[which(new_matrix$chrom == i),"ave_exp"]>0,"red","blue"),border=NA)
    }
    for(i in 1:22){
      barplot(new_matrix_MCL[which(new_matrix_MCL$chrom == i),"diff"], main=paste0("chr",i), ylim=c(-1,1), ylab="Mean PC1 values (MCL-normal)", col=ifelse(new_matrix_MCL[which(new_matrix_MCL$chrom == i),"diff"]>0,"red","blue"),border=NA)
    }

  #Boxplots
    data = rbind(new_matrix_CLL,new_matrix_MCL)
    data$chrom = factor(data$chrom)
    library(ggplot2)
    ggplot(data, aes(x=chrom, y=diff, fill=cell.grp)) + geom_boxplot() + facet_wrap(~chrom, scale="free") + scale_fill_manual(values=c("yellow","orange")) + geom_hline(yintercept=0,lty=2) + ggtitle("All bins")
    ggplot(data, aes(x=chrom, y=diff, fill=cell.grp)) + geom_boxplot() + scale_fill_manual(values=c("yellow","orange")) + geom_hline(yintercept=0,lty=2) + ggtitle("All bins") + theme_minimal()
    ggplot(subset(data, cell.grp == "CLL"), aes(y=diff, x=chrom)) + geom_boxplot(fill="yellow") + geom_hline(yintercept=0,lty=2) + ggtitle("All bins - CLL vs normals") 
    ggplot(subset(data, cell.grp == "MCL"), aes(y=diff, x=chrom)) + geom_boxplot(fill="orange") + geom_hline(yintercept=0,lty=2) + ggtitle("All bins - MCL vs normals")

    boxplot(diff~chrom, data=new_matrix_CLL,pch=19,ylim=c(-0.06,0.06),col="yellow",main="CLL vs normals"); abline(h=0,lty=2)
    boxplot(diff~chrom, data=new_matrix_MCL,pch=19,ylim=c(-0.06,0.06),col="orange",main="MCL vs normals"); abline(h=0,lty=2)

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
        plot(exp_hist[[1]], type="n", xlim=c(-0.1,0.1), main=paste0("chr",i))
        for(j in 1:length(experimental)){ lines(exp_hist[[j]], col="yellow4")}
        for(k in 1:length(control)){ lines(ctr_hist[[k]], col="grey")}
      }
  
    #Mean MCL and Mean Normal 
      i=0                  
      for(i in 1:22){
        exp = density(subset(new_matrix_CLL, chrom == i)$ave_exp, na.rm=T)
        ctr = density(subset(new_matrix_CLL, chrom == i)$ave_control, na.rm=T)
        plot(exp, xlim=c(-0.1,0.1), main=paste0("chr",i), type="n")
        lines(exp, col="yellow4")
        lines(ctr, col="grey")
        }


                   
#Grouped by Compartment State
  #Pie Chart
  new_matrix = new_matrix_MCL
  pie_chart = data.frame(table(new_matrix[,"cs_status"]))
  pie_chart$percentage = pie_chart[,2]/sum(pie_chart[,2])*100
  pie(pie_chart$percentage, labels = NA, col=c("grey","blue","red","black"))
  pie_chart
  
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
  boxplot(diff~chrom, data=new_matrix[which(new_matrix$cs_status %in% c("A-B","B-A")),],pch=19,col="yellow",main="CLL vs normals")
  abline(h=0,lty=2)

#-----------------------------Interaction Score vs CS-------------------------

  setwd("~/Google Drive/Data/Interaction Scores (Genome Wide)")
  
  #Merged matrix
        Matrix_Int_CS = function(cell_group){
          IntScores = list()
          for(i in 1:22){
            IntScores[[i]] = read.table(paste0(cell_group,".chr.",i,".int.diff.short.minus.long.txt"))
          }
          IntScores <<- IntScores
        }
        IntScores_MCL = do.call(rbind,Matrix_Int_CS("MCL"))
        IntScores_CLL = do.call(rbind,Matrix_Int_CS("CLL"))
        
        tmp = data.frame(matrix_1Mb, Int_MCL = unlist(IntScores_MCL), Int_CLL = unlist(IntScores_CLL))
        matrix_1Mb_MCL = AnalysisMatrix("MCL",tmp)[,c("bin","start","end","chrom", "Int_MCL","ave_exp","ave_control","diff")]
        matrix_1Mb_CLL = AnalysisMatrix("CLL",tmp)[,c("bin","start","end","chrom","Int_CLL","ave_exp","ave_control","diff")]
        
  #Individual samples matrix
  list.samples = c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828",
                   "mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182")
  
  persample = function(chr){
    IntPerChrom= list()
    for(j in 1:length(list.samples)){
      IntPerChrom[[list.samples[j]]] = read.table(paste0(list.samples[j],".chr.",i,".int.diff.short.minus.long.txt"))
    }
    IntPerChrom <<- IntPerChrom
  }
  
  IntAllChrom = list()
  for(i in 1:22){
    IntAllChrom[[i]] = do.call(cbind,persample(i))
  }
  IntAllChrom = do.call(rbind, IntAllChrom)
  colnames(IntAllChrom) = paste0(list.samples,"_Int")  #int scores ALL 2887
  CS = matrix_1Mb[,-c(1:4)]  
  colnames(CS) = paste0(names(CS),"_CS")
  MeanControl = apply(matrix_1Mb, 1, function(x) median(as.numeric(x[control])))  
  CSDiff = t(apply(matrix_1Mb, 1, function(x)  x[list.samples]-mean(x[control])))
  colnames(CSDiff) =  paste0(list.samples, "_CSDiff")
  
  tmp = data.frame(matrix_1Mb[,1:4],MeanControl, CS, CSDiff, IntAllChrom)
  DataImputation2 = function(cell_group,selected_matrix){
    
    ifelse(cell_group == "MCL", 
           experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
           experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
    
    control <- c("NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3","PBC_1","PBC_2","PBC_3")
    
    NAs.exp<-rowSums(is.na(selected_matrix[,paste0(experimental,"_CS")]))
    NAs.control<-rowSums(is.na(selected_matrix[,paste0(control,"_CS")]))
    rows.to.include<-which(NAs.exp<=(length(experimental)-3) & NAs.control<=(length(control)-3))
    selected_matrix <<- selected_matrix[rows.to.include,!(names(selected_matrix) %in% paste0(control,"_CS"))]
  }
  matrix_1Mb_MCL5 = DataImputation2("MCL",tmp)
  matrix_1Mb_MCL5 = matrix_1Mb_MCL5[,!grepl("*CLL",names(matrix_1Mb_MCL5))]
  matrix_1Mb_CLL7 = DataImputation2("CLL",tmp)
  matrix_1Mb_CLL7 = matrix_1Mb_CLL7[,!grepl("*MCL",names(matrix_1Mb_CLL7))]
  
  
  #-------Interaction Counts vs MCL/CLL/Normal CS
  
  Plottin = function(matrix, x, y, chromosome, cell_group){
    plot(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],pch=19,cex=0.5, xlab="Interaction Scores", ylab=c("Compartment Score",cell_group),main=paste("chromosome",chromosome))
  }
  #Plottin(matrix_1Mb_MCL, "Int_MCL", "ave_exp", 11, "MCL" )
  #Plottin(matrix_1Mb_MCL, "Int_MCL", "diff", 14, "MCL)
  
  Plottin2 = function(matrix, x, y, chromosome, cell_group){
    plot(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],pch=19,cex=0.7, xlab=c("Compartment Score","Normal"), ylab=c("Compartment Score",cell_group),main=paste("chromosome",chromosome))
  }
  #Plottin2(matrix_1Mb_MCL, "ave_control", "diff", 11, "MCL-normal" )
  
  
  #Triple
  par(mfrow=c(3,1))
  
  Plottin3_orig = function(matrix, x, y, chromosome, ylabel, cell_group,ylimdown=NULL,ylimup=NULL){
    plot(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],pch=19,cex=0.5, ylim=c(ylimdown,ylimup),xlab="Genomic Location", ylab=c(ylabel,cell_group),main=paste("chromosome",chromosome),col=colour)
    lines(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],lwd=2,col=ifelse(matrix[which(matrix$chrom == chromosome),y]>0,"red","blue"));abline(h=0,lty=2)
  }
  
  Plottin3 = function(matrix, x, y, chromosome, ylabel, cell_group,ylimdown=NULL,ylimup=NULL){
    plot(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],pch=19,cex=0.5,ylim=c(ylimdown,ylimup),xlab="Genomic Location", ylab=c(ylabel,cell_group),main=paste("chromosome",chromosome),col=ifelse(matrix[which(matrix$chrom == chromosome),y]>0,"red","blue"))
    abline(h=0, lty=2)
    segments(x0=matrix[which(matrix$chrom == chromosome),x][-length(which(matrix$chrom == chromosome))],
              y0=matrix[which(matrix$chrom == chromosome),y][-length(which(matrix$chrom == chromosome))],
              x1=matrix[which(matrix$chrom == chromosome),x][-1],
              y1=matrix[which(matrix$chrom == chromosome),y][-1], col=ifelse(matrix[which(matrix$chrom == chromosome),y][-1]>0,"red","blue"))

    }
  Plottin3(matrix_1Mb_MCL, "bin", "Int_MCL", 11, "Interaction Scores", "MCL",-5,5)  
  Plottin3(matrix_1Mb_MCL, "bin", "diff", 11, "Compartment Scores", "MCL-normal")  
  Plottin3(matrix_1Mb_MCL, "bin", "ave_exp", 11, "Compartment Scores", "MCL")
  
  Plottin3(matrix_1Mb_MCL, "bin", "Int_MCL", 14, "Interaction Scores", "MCL",-5,5)  
  Plottin3(matrix_1Mb_MCL, "bin", "diff", 14, "Compartment Scores", "MCL-normal")  
  Plottin3(matrix_1Mb_MCL, "bin", "ave_exp", 14, "Compartment Scores", "MCL")
  
  Plottin3(matrix_1Mb_CLL, "bin", "Int_CLL", 11, "Interaction Scores", "CLL",-5,5)  
  Plottin3(matrix_1Mb_CLL, "bin", "diff", 11, "Compartment Scores", "CLL-normal")  
  Plottin3(matrix_1Mb_CLL, "bin", "ave_exp", 11, "Compartment Scores", "CLL")
  
  Plottin3(matrix_1Mb_CLL, "bin", "Int_CLL", 14, "Interaction Scores", "CLL",-5,5)  
  Plottin3(matrix_1Mb_CLL, "bin", "diff", 14, "Compartment Scores", "CLL-normal")  
  Plottin3(matrix_1Mb_CLL, "bin", "ave_exp", 14, "Compartment Scores", "CLL")
  
#------Per sample
  PlotFig1a = function(matrix, x, y, chromosome, cell_group){
    plot(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],pch=19,cex=0.5, xlab=c("Interaction Scores",cell_group), ylab=c("Compartment Score",cell_group),main=paste("chromosome",chromosome))
  }
  
  Fig1a = function(cell_group, matrix, chromosome){
  
  ifelse(cell_group == "MCL", 
         experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
         experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
  
  for(i in 1:length(experimental)){
    PlotFig1a(matrix, paste0(experimental[i],"_Int"), paste0(experimental[i],"_CS"), chromosome, experimental[i])
  }
  }

  
  
  PlotFig2 = function(matrix, x, y, chromosome, cell_group){
    plot(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],pch=19,cex=0.5, xlab=c("Interaction Scores",cell_group), ylab=c("Compartment Score",paste0(cell_group,"-normal")),main=paste("chromosome",chromosome))
  }
  
  Fig2 = function(cell_group, matrix, chromosome){
    
    ifelse(cell_group == "MCL", 
           experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
           experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
    
    for(i in 1:length(experimental)){
      PlotFig2(matrix, paste0(experimental[i],"_Int"), paste0(experimental[i],"_CSDiff"), chromosome, experimental[i])
    }
  }
  
  PlotFig3 = function(matrix, x, y, chromosome, cell_group){
    plot(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],pch=19,cex=0.5, xlab=c("Compartment Score","Normal"), ylab=c("Compartment Score",paste0(cell_group,"-normal")),main=paste("chromosome",chromosome))
  }
  Fig3 = function(cell_group, matrix, chromosome){
    
    ifelse(cell_group == "MCL", 
           experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
           experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
    
    for(i in 1:length(experimental)){
      PlotFig3(matrix, "MeanControl", paste0(experimental[i],"_CSDiff"), chromosome, experimental[i])
    }
  }
  
  
  
  
  
  
  
  
  par(mfrow=c(3,1))
  
  
  PlotFig4 = function(matrix, x, y, chromosome, ylabel, cell_group,ylimdown=NULL,ylimup=NULL){
    
    
    
    plot(matrix[which(matrix$chrom == chromosome),x], matrix[which(matrix$chrom == chromosome),y],pch=19,cex=0.25,ylim=c(ylimdown,ylimup),xlab="Genomic Location", ylab=c(ylabel,cell_group),main=paste("chromosome",chromosome),col=ifelse(matrix[which(matrix$chrom == chromosome),y]>0,"red","blue"))
    abline(h=0, lty=2)
    segments(x0=matrix[which(matrix$chrom == chromosome),x][-length(which(matrix$chrom == chromosome))],
             y0=matrix[which(matrix$chrom == chromosome),y][-length(which(matrix$chrom == chromosome))],
             x1=matrix[which(matrix$chrom == chromosome),x][-1],
             y1=matrix[which(matrix$chrom == chromosome),y][-1], col=ifelse(matrix[which(matrix$chrom == chromosome),y][-1]>0,"red","blue"),lwd=0.5)
    
  }
  Fig4 = function(cell_group, matrix, y, chromosome,ylabel,ylimdown=NULL,ylimup=NULL){
    
    ifelse(cell_group == "MCL", 
           experimental <- c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828"),
           experimental <- c("mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182"))
    
    for(i in 1:length(experimental)){
      PlotFig4(matrix, "bin", paste0(experimental[i],"_",y), chromosome, ylabel, experimental[i], ylimdown, ylimup)
    }
  }
  
  Fig4("MCL", matrix_1Mb_MCL5, "Int", 11, "Interaction Score")
  
  
  