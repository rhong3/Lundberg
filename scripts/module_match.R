# Module match
library(readr)
library(ComplexHeatmap)

cancerss = list(c('BRCA',4), c('LUAD',3), c('LSCC',3), c('OV',2), c('CO',2), c('HN',3), c('UCEC',4), c('GBM',4), c('CCRCC',3))
for (can in cancerss){
  prot = read_csv(paste("Results/", can[1], "_prot_groups.csv", sep=''))
  trans = read_csv(paste("Results/", can[1], "_RNA_groups.csv", sep=''))
  prot = prot[,c(1:2)]
  trans = trans[, c(1:2)]
  colnames(prot) = c('Gene', 'cluster')
  colnames(trans) = c('Gene', 'cluster')
  mgg = merge(prot, trans, by="Gene", suffixes = c(".prot",".RNA"))
  print(can[1])
  for (clustp in 1:can[2]){
    for (clustt in 1:can[2]){
      print(paste("Prot:", clustp, " RNA:", clustt, sep=""))
      print(nrow(mgg[mgg$cluster.prot==clustp & mgg$cluster.RNA==clustt,])/nrow(mgg))
    }
  }
  write.csv(mgg, paste("Results/", can[1], "_groups.csv", sep=''), row.names = FALSE)
}


