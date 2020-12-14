# Module match
library(readr)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(gridExtra)

cancerss = list(c('BRCA',4), c('LUAD',3), c('LSCC',3), c('OV',2), c('CO',2), c('HN',3), c('UCEC',4), c('GBM',4), c('CCRCC',3))
protcol = colorRampPalette(brewer.pal(4,"Dark2"))
RNAcol = colorRampPalette(brewer.pal(4,"Set2"))

for (can in cancerss){
  prot = read_csv(paste("Results/", can[1], "_prot_groups.csv", sep=''))
  trans = read_csv(paste("Results/", can[1], "_RNA_groups.csv", sep=''))
  prot.short = prot[,c(1:2)]
  trans.short = trans[, c(1:2)]
  colnames(prot.short) = c('Gene', 'cluster')
  colnames(trans.short) = c('Gene', 'cluster')
  
  mgg = merge(prot.short, trans.short, by="Gene", suffixes = c(".prot",".RNA"))
  prot = prot[prot$X1 %in% mgg$Gene, c(c('X1', 'clusters'), mgg$Gene)]
  trans = trans[trans$X1 %in% mgg$Gene, c(c('X1', 'clusters'), mgg$Gene)]
  row.names(prot) = prot$X1
  row.names(trans) = trans$X1
  
  prot[is.na(prot)] = 0
  trans[is.na(trans)] = 0
  
  print(can[1])
  for (clustp in 1:can[2]){
    for (clustt in 1:can[2]){
      print(paste("Prot:", clustp, " RNA:", clustt, sep=""))
      print(nrow(mgg[mgg$cluster.prot==clustp & mgg$cluster.RNA==clustt,])/nrow(mgg))
    }
  }
  write.csv(mgg, paste("Results/", can[1], "_groups.csv", sep=''), row.names = FALSE)
  
  
  protcollist = protcol(can[2])
  RNAcollist =RNAcol(can[2])
  names(protcollist) = c(1:can[2])
  names(RNAcollist) = c(1:can[2])
  
  anno = HeatmapAnnotation(prot.cluster = mgg$cluster.prot,
                           RNA.cluster = mgg$cluster.RNA,
                           col=list(prot.cluster=protcollist,
                                    RNA.cluster=RNAcollist))
  protplot = Heatmap(as.matrix(prot[,-c(1,2)]), top_annotation = anno, column_title = paste(can[1], ' proteomics', sep=''), show_row_names=FALSE, show_column_names=FALSE)
  RNAplot = Heatmap(as.matrix(trans[,-c(1,2)]), top_annotation = anno, column_title = paste(can[1], ' transcriptomics', sep=''), show_row_names=FALSE, show_column_names=FALSE)

  pdf(file=paste('Results/', can[1], "-prot-corr-groups.pdf", sep=''),
      width=8,height=8)
  draw(protplot)
  dev.off()
  
  pdf(file=paste('Results/', can[1], "-RNA-corr-groups.pdf", sep=''),
      width=8,height=8)
  draw(RNAplot)
  dev.off()
  
}






