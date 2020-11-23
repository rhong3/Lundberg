library(readr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pheatmap)

cancers = c('BRCA', 'LUAD', 'LSCC', 'OV', 'CO', 'HN', 'UCEC', 'GBM', 'CCRCC')
ensg_ccdprotein <- read_csv("ensg_ccdprotein.csv")
ensg_ccdprotein_transcript_regulated <- read_csv("ensg_ccdprotein_transcript_regulated.csv")
ensg_ccdprotein_nontranscript_regulated <- read_csv("ensg_ccdprotein_nontranscript_regulated.csv")

# Prot-RNA corr
all = data.frame(matrix(ncol = 3, nrow = 0))
colnames(all) <- c('correlation coefficient', 'transcript regulated', 'type')
summ = data.frame(matrix(ncol = 7, nrow = 0))
colnames(summ) <- c('type', "proteins intersect", "RNA intersect", "all sample intersect", "tumor sample intersect", "tumor gene intersect", "corr genes")

for (can in cancers){
  prot = read_csv(paste(can, "_proteomics.csv", sep=''))
  trans = read_csv(paste(can, "_transcriptomics.csv", sep=''))
  prot[, 2:ncol(prot)] = prot[, order(colnames(prot[, 2:ncol(prot)]))]
  trans[, 2:ncol(trans)] = trans[, order(colnames(trans[, 2:ncol(trans)]))]
  prot[, 2:ncol(prot)] <- sapply(prot[, 2:ncol(prot)] , as.numeric)
  trans[, 2:ncol(trans)] <- sapply(trans[, 2:ncol(trans)] , as.numeric)
  
  prot = prot[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(prot)))]
  trans = trans[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(trans)))]
  
  print(can)
  aa = ncol(prot)-1
  bb = ncol(trans)-1
  print(paste('proteins intersect: ', aa, sep=''))
  print(paste('RNA intersect: ', bb, sep=''))
  
  pr = prot$Patient_ID
  tr = trans$Patient_ID
  
  row.names(prot) = pr
  prot = prot[intersect(pr, tr), ]
  row.names(trans) = tr
  trans = trans[intersect(pr, tr), ]
  
  prot = prot[order(prot$Patient_ID),]
  trans = trans[order(trans$Patient_ID),]
  
  cc = nrow(prot)
  print(paste('all sample intersect: ', cc, sep=''))
  
  prot.tumor = prot[!grepl("\\.N", prot$Patient_ID), ]
  trans.tumor = trans[!grepl("\\.N", trans$Patient_ID), ]
  prot.tumor = prot[!grepl("\\.C", prot$Patient_ID), ]
  trans.tumor = trans[!grepl("\\.C", trans$Patient_ID), ]
  
  dd = nrow(prot.tumor)
  print(paste('tumor sample intersect: ', dd, sep=''))
  
  prot.tumor.x = prot.tumor[,intersect(colnames(prot.tumor), colnames(trans.tumor))]
  trans.tumor.x = trans.tumor[,intersect(colnames(prot.tumor), colnames(trans.tumor))]
  
  ee = ncol(prot.tumor.x)-1
  print(paste('tumor gene intersect: ', ee, sep=''))
  
  inter_cor = diag(cor(as.matrix(prot.tumor.x[, 2:ncol(prot.tumor.x)]), as.matrix(trans.tumor.x[, 2:ncol(trans.tumor.x)]),  method = "spearman", use="pairwise.complete.obs"))
  inter_cor = data.frame(inter_cor)
  inter_cor$regulated = row.names(inter_cor) %in% ensg_ccdprotein_transcript_regulated$SYMBOL
  inter_cor = na.omit(inter_cor)
  colnames(inter_cor) = c('correlation coefficient', 'transcript regulated')
  inter_cor$type = can
  
  ff = nrow(inter_cor)
  print(paste('corr genes: ', ff, sep=''))
  
  all = rbind(all, inter_cor)
  summ[nrow(summ)+1, ] = c(can, aa, bb, cc, dd, ee, ff)
  
  p = ggplot(inter_cor, aes(x=`correlation coefficient`, fill=`transcript regulated`, color=`transcript regulated`)) +
    geom_histogram(position="identity", alpha=0.5, binwidth=0.05)+scale_color_brewer(palette="Dark2")+
    scale_fill_brewer(palette="Dark2")+
    labs(title=paste(can, "protein and mRNA Spearman correlation", sep=' '))+ 
    theme_classic()+
    theme(legend.position="top", plot.title = element_text(hjust = 0.5)) 
  pdf(file=paste('Results/', can, "_protein-mRNA-corr.pdf", sep=''),
      width=10,height=8)
  grid.arrange(p,nrow=1, ncol=1)
  dev.off()
}

p = ggplot(all, aes(x=`correlation coefficient`, fill=`transcript regulated`, color=`transcript regulated`)) +
  geom_histogram(position="identity", alpha=0.5, binwidth=0.05)+scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  labs(title=paste('All types', "protein and mRNA Spearman correlation", sep=' '))+ 
  theme_classic()+
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) 
pdf(file=paste('Results/', 'All types', "_protein-mRNA-corr.pdf", sep=''),
    width=10,height=8)
grid.arrange(p,nrow=1, ncol=1)
dev.off()

write.csv(all, 'Results/protein-mRNA-corr.csv', row.names = TRUE)
write.csv(summ, 'Results/summary_protein-mRNA-corr.csv', row.names = FALSE)


# # Corr within group
# for (can in cancers){
#   prot = read_csv(paste(can, "_proteomics.csv", sep=''))
#   trans = read_csv(paste(can, "_transcriptomics.csv", sep=''))
#   prot[, 2:ncol(prot)] = prot[, order(colnames(prot[, 2:ncol(prot)]))]
#   trans[, 2:ncol(trans)] = trans[, order(colnames(trans[, 2:ncol(trans)]))]
#   prot[, 2:ncol(prot)] <- sapply(prot[, 2:ncol(prot)] , as.numeric)
#   trans[, 2:ncol(trans)] <- sapply(trans[, 2:ncol(trans)] , as.numeric)
#   
#   prot.tumor = prot[!grepl("\\.N", prot$Patient_ID), ]
#   trans.tumor = trans[!grepl("\\.N", trans$Patient_ID), ]
#   prot.tumor = prot[!grepl("\\.C", prot$Patient_ID), ]
#   trans.tumor = trans[!grepl("\\.C", trans$Patient_ID), ]
#   prot.tumor = prot.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(prot.tumor)))]
#   trans.tumor = trans.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(trans.tumor)))]
#   matt.prot = cor(prot.tumor[, 2:ncol(prot.tumor)],  method = "spearman", use="pairwise.complete.obs")
#   matt.trans = cor(trans.tumor[, 2:ncol(trans.tumor)],  method = "spearman", use="pairwise.complete.obs")
#   write.csv(matt.prot, paste('Results/', can, '_protein-corr-ALL.csv', sep=''), row.names = TRUE)
#   write.csv(matt.trans, paste('Results/', can, '_RNA-corr-ALL.csv', sep=''), row.names = TRUE)
#   matt.prot[is.na(matt.prot)] <- 0
#   matt.trans[is.na(matt.trans)] <- 0
#   pp = pheatmap(matt.prot, show_rownames=FALSE, show_colnames=FALSE, main=paste(can, ' proteomics Spearman correlation ALL', sep=''))
#   pt = pheatmap(matt.trans, show_rownames=FALSE, show_colnames=FALSE, main=paste(can, ' transcriptomics Spearman correlation ALL', sep=''))
#   pdf(file=paste('Results/', can, "_ALL-internal-corr.pdf", sep=''),
#       width=20,height=11)
#   grid.arrange(pp$gtable, pt$gtable, nrow=1, ncol=2)
#   dev.off()
#   
#   
#   prot.tumor = prot[!grepl("\\.N", prot$Patient_ID), ]
#   trans.tumor = trans[!grepl("\\.N", trans$Patient_ID), ]
#   prot.tumor = prot[!grepl("\\.C", prot$Patient_ID), ]
#   trans.tumor = trans[!grepl("\\.C", trans$Patient_ID), ]
#   prot.tumor = prot.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein_nontranscript_regulated$SYMBOL, colnames(prot.tumor)))]
#   trans.tumor = trans.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein_nontranscript_regulated$SYMBOL, colnames(trans.tumor)))]
#   matt.prot = cor(prot.tumor[, 2:ncol(prot.tumor)],  method = "spearman", use="pairwise.complete.obs")
#   matt.trans = cor(trans.tumor[, 2:ncol(trans.tumor)],  method = "spearman", use="pairwise.complete.obs")
#   write.csv(matt.prot, paste('Results/', can, '_protein-corr-nontranscript.csv', sep=''), row.names = TRUE)
#   write.csv(matt.trans, paste('Results/', can, '_RNA-corr-nontranscript.csv', sep=''), row.names = TRUE)
#   matt.prot[is.na(matt.prot)] <- 0
#   matt.trans[is.na(matt.trans)] <- 0
#   pp = pheatmap(matt.prot, show_rownames=FALSE, show_colnames=FALSE, main=paste(can, ' proteomics Spearman correlation nontranscript', sep=''))
#   pt = pheatmap(matt.trans, show_rownames=FALSE, show_colnames=FALSE, main=paste(can, ' transcriptomics Spearman correlation nontranscript', sep=''))
#   pdf(file=paste('Results/', can, "_nontranscript-internal-corr.pdf", sep=''),
#       width=20,height=11)
#   grid.arrange(pp$gtable, pt$gtable, nrow=1, ncol=2)
#   dev.off()
#   
#   prot.tumor = prot[!grepl("\\.N", prot$Patient_ID), ]
#   trans.tumor = trans[!grepl("\\.N", trans$Patient_ID), ]
#   prot.tumor = prot[!grepl("\\.C", prot$Patient_ID), ]
#   trans.tumor = trans[!grepl("\\.C", trans$Patient_ID), ]
#   prot.tumor = prot.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein_transcript_regulated$SYMBOL, colnames(prot.tumor)))]
#   trans.tumor = trans.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein_transcript_regulated$SYMBOL, colnames(trans.tumor)))]
#   matt.prot = cor(prot.tumor[, 2:ncol(prot.tumor)],  method = "spearman", use="pairwise.complete.obs")
#   matt.trans = cor(trans.tumor[, 2:ncol(trans.tumor)],  method = "spearman", use="pairwise.complete.obs")
#   write.csv(matt.prot, paste('Results/', can, '_protein-corr-transcript.csv', sep=''), row.names = TRUE)
#   write.csv(matt.trans, paste('Results/', can, '_RNA-corr-transcript.csv', sep=''), row.names = TRUE)
#   matt.prot[is.na(matt.prot)] <- 0
#   matt.trans[is.na(matt.trans)] <- 0
#   pp = pheatmap(matt.prot, show_rownames=FALSE, show_colnames=FALSE, main=paste(can, ' proteomics Spearman correlation transcript', sep=''))
#   pt = pheatmap(matt.trans, show_rownames=FALSE, show_colnames=FALSE, main=paste(can, ' transcriptomics Spearman correlation transcript', sep=''))
#   pdf(file=paste('Results/', can, "_transcript-internal-corr.pdf", sep=''),
#       width=20,height=11)
#   grid.arrange(pp$gtable, pt$gtable, nrow=1, ncol=2)
#   dev.off()
# }
# 
