#--------------------------------------------------------------------
#
# CREATE MATRIX 
#
# Creates the matrix for coverage corrected files for every viewpoint
# file selected
#
#--------------------------------------------------------------------

Create_Corr_Matrix = function(viewpt){
  setwd("C:/Users/admin/Google Drive (kimberlyquililan@gmail.com)/Data/CoverageCorrected")
  list.samples = c("cMCL_1064","cMCL_568","nnMCL_309","nnMCL_817", "nnMCL_828",
                   "mCLL_3","mCLL_110","mCLL_1228","mCLL_1525","mCLL_1532","uCLL_12","uCLL_182",
                   "NBC_1","NBC_2","NBC_3","GCBC_1","GCBC_2","GCBC_3","MBC_1","MBC_2","MBC_3",
                   "PBC_1","PBC_2","PBC_3")
  
  tmp = read.table(paste0(list.samples[1],".data.chr11.chr14.breakpoint.bins.coverage.corrected.txt"), header=T)
  new_matrix = tmp[,c("bin","chrom","start","end")]
  new_matrix$viewpt = tmp[,viewpt]
  names(new_matrix)[length(names(new_matrix))]<- list.samples[1]
  rm(tmp)
  
  for(i in 2:24){
    tmp = read.table(paste0(list.samples[i],".data.chr11.chr14.breakpoint.bins.coverage.corrected.txt"), header=T)
    new_matrix$viewpt = tmp[,viewpt]
    names(new_matrix)[length(names(new_matrix))]<- list.samples[i]
    rm(tmp)
  }
  
  new_matrix <<- new_matrix
}

Create_Norm_Matrix = function(viewpt){
  setwd("/Users/kimquililan/Google Drive/Data/Normalized2/viewpoints_1Mb_norm")
  list.samples = c("cMCL1064","cMCL568","nnMCL309","nnMCL817", "nnMCL828",
                   "mCLL3","mCLL110","mCLL1228","mCLL1525","mCLL1532","uCLL12","uCLL182",
                   "NBC1","NBC2","NBC3","GCBC1","GCBC2","GCBC3","MBC1","MBC2","MBC3",
                   "PBC1","PBC2","PBC3")
  
  tmp = read.table(paste0(list.samples[1],".norm.counts.viewpoints.chr11.14.txt"), header=T)
  new_matrix = tmp[,c("bin","chr","start","end")]
  colnames(new_matrix) = c("bin","chrom","start","end")
  new_matrix$viewpt = tmp[,viewpt]
  names(new_matrix)[length(names(new_matrix))]<- list.samples[1]
  rm(tmp)
  
  for(i in 2:24){
    tmp = read.table(paste0(list.samples[i],".norm.counts.viewpoints.chr11.14.txt"), header=T)
    new_matrix$viewpt = tmp[,viewpt]
    names(new_matrix)[length(names(new_matrix))]<- list.samples[i]
    rm(tmp)
  }
  
  new_matrix <<- new_matrix
}

viewpoints = c("chr11.68.to.69",	"chr11.69.to.70",	"chr11.70.to.71",	"chr14.104.to.105",	"chr14.105.to.106","chr14.106.to.107")
for(i in 1:6){
  Create_Norm_Matrix(viewpoints[i])
  matrix = new_matrix
  save(matrix, file=paste0(viewpoints[i],".Rdata"))
}
