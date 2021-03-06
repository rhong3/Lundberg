library(readr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)

cancers = c('BRCA', 'LUAD', 'LSCC', 'OV', 'CO', 'HN', 'UCEC', 'GBM', 'CCRCC')
ensg_ccdprotein <- read_csv("ensg_ccdprotein.csv")
ensg_ccdprotein_transcript_regulated <- read_csv("ensg_ccdprotein_transcript_regulated.csv")
ensg_ccdprotein_nontranscript_regulated <- read_csv("ensg_ccdprotein_nontranscript_regulated.csv")

# # Prot-RNA corr
# all = data.frame(matrix(ncol = 3, nrow = 0))
# colnames(all) <- c('correlation coefficient', 'transcript regulated', 'type')
# summ = data.frame(matrix(ncol = 7, nrow = 0))
# colnames(summ) <- c('type', "proteins intersect", "RNA intersect", "all sample intersect", "tumor sample intersect", "tumor gene intersect", "corr genes")
# 
# for (can in cancers){
#   prot = read_csv(paste(can, "_proteomics.csv", sep=''))
#   trans = read_csv(paste(can, "_transcriptomics.csv", sep=''))
#   prot[, 2:ncol(prot)] = prot[, order(colnames(prot[, 2:ncol(prot)]))]
#   trans[, 2:ncol(trans)] = trans[, order(colnames(trans[, 2:ncol(trans)]))]
#   prot[, 2:ncol(prot)] <- sapply(prot[, 2:ncol(prot)] , as.numeric)
#   trans[, 2:ncol(trans)] <- sapply(trans[, 2:ncol(trans)] , as.numeric)
# 
#   prot = prot[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(prot)))]
#   trans = trans[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(trans)))]
# 
#   print(can)
#   aa = ncol(prot)-1
#   bb = ncol(trans)-1
#   print(paste('proteins intersect: ', aa, sep=''))
#   print(paste('RNA intersect: ', bb, sep=''))
# 
#   pr = prot$Patient_ID
#   tr = trans$Patient_ID
# 
#   row.names(prot) = pr
#   prot = prot[intersect(pr, tr), ]
#   row.names(trans) = tr
#   trans = trans[intersect(pr, tr), ]
# 
#   prot = prot[order(prot$Patient_ID),]
#   trans = trans[order(trans$Patient_ID),]
# 
#   cc = nrow(prot)
#   print(paste('all sample intersect: ', cc, sep=''))
# 
#   prot.tumor = prot[!grepl("\\.N", prot$Patient_ID), ]
#   trans.tumor = trans[!grepl("\\.N", trans$Patient_ID), ]
#   prot.tumor = prot[!grepl("\\.C", prot$Patient_ID), ]
#   trans.tumor = trans[!grepl("\\.C", trans$Patient_ID), ]
# 
#   dd = nrow(prot.tumor)
#   print(paste('tumor sample intersect: ', dd, sep=''))
# 
#   prot.tumor.x = prot.tumor[,intersect(colnames(prot.tumor), colnames(trans.tumor))]
#   trans.tumor.x = trans.tumor[,intersect(colnames(prot.tumor), colnames(trans.tumor))]
# 
#   ee = ncol(prot.tumor.x)-1
#   print(paste('tumor gene intersect: ', ee, sep=''))
# 
#   inter_cor = diag(cor(as.matrix(prot.tumor.x[, 2:ncol(prot.tumor.x)]), as.matrix(trans.tumor.x[, 2:ncol(trans.tumor.x)]),  method = "spearman", use="pairwise.complete.obs"))
#   inter_cor = data.frame(inter_cor)
#   inter_cor$regulated = row.names(inter_cor) %in% ensg_ccdprotein_transcript_regulated$SYMBOL
#   inter_cor = na.omit(inter_cor)
#   colnames(inter_cor) = c('correlation coefficient', 'transcript regulated')
#   inter_cor$type = can
# 
#   ff = nrow(inter_cor)
#   print(paste('corr genes: ', ff, sep=''))
# 
#   all = rbind(all, inter_cor)
#   summ[nrow(summ)+1, ] = c(can, aa, bb, cc, dd, ee, ff)
# 
#   p = ggplot(inter_cor, aes(x=`correlation coefficient`, fill=`transcript regulated`, color=`transcript regulated`)) +
#     geom_histogram(position="identity", alpha=0.5, binwidth=0.05)+scale_color_brewer(palette="Dark2")+
#     scale_fill_brewer(palette="Dark2")+
#     labs(title=paste(can, "protein and mRNA Spearman correlation", sep=' '))+
#     theme_classic()+
#     theme(legend.position="top", plot.title = element_text(hjust = 0.5))
#   pdf(file=paste('Results/', can, "_protein-mRNA-corr.pdf", sep=''),
#       width=10,height=8)
#   grid.arrange(p,nrow=1, ncol=1)
#   dev.off()
# }
# 
# p = ggplot(all, aes(x=`correlation coefficient`, fill=`transcript regulated`, color=`transcript regulated`)) +
#   geom_histogram(position="identity", alpha=0.5, binwidth=0.05)+scale_color_brewer(palette="Dark2")+
#   scale_fill_brewer(palette="Dark2")+
#   labs(title=paste('All types', "protein and mRNA Spearman correlation", sep=' '))+
#   theme_classic()+
#   theme(legend.position="top", plot.title = element_text(hjust = 0.5))
# pdf(file=paste('Results/', 'All types', "_protein-mRNA-corr.pdf", sep=''),
#     width=10,height=8)
# grid.arrange(p,nrow=1, ncol=1)
# dev.off()
# 
# write.csv(all, 'Results/protein-mRNA-corr.csv', row.names = TRUE)
# write.csv(summ, 'Results/summary_protein-mRNA-corr.csv', row.names = FALSE)






cancerss = list(c('BRCA',4,4), c('LUAD',3,3), c('LSCC',3,3), c('OV',2,2), c('CO',2,2), c('HN',3,3), c('UCEC',4,4), c('GBM',4,4), c('CCRCC',3,3))


# Corr within group
for (can in cancerss){
  prot = read_csv(paste(can[1], "_proteomics.csv", sep=''))
  trans = read_csv(paste(can[1], "_transcriptomics.csv", sep=''))
  prot[, 2:ncol(prot)] = prot[, order(colnames(prot[, 2:ncol(prot)]))]
  trans[, 2:ncol(trans)] = trans[, order(colnames(trans[, 2:ncol(trans)]))]
  prot[, 2:ncol(prot)] <- sapply(prot[, 2:ncol(prot)] , as.numeric)
  trans[, 2:ncol(trans)] <- sapply(trans[, 2:ncol(trans)] , as.numeric)

  prot.tumor = prot[!grepl("\\.N", prot$Patient_ID), ]
  trans.tumor = trans[!grepl("\\.N", trans$Patient_ID), ]
  prot.tumor = prot[!grepl("\\.C", prot$Patient_ID), ]
  trans.tumor = trans[!grepl("\\.C", trans$Patient_ID), ]
  prot.tumor = prot.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(prot.tumor)))]
  trans.tumor = trans.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(trans.tumor)))]
  matt.prot = cor(prot.tumor[, 2:ncol(prot.tumor)],  method = "spearman", use="pairwise.complete.obs")
  matt.trans = cor(trans.tumor[, 2:ncol(trans.tumor)],  method = "spearman", use="pairwise.complete.obs")
  write.csv(matt.prot, paste('Results/', can[1], '_protein-corr-ALL.csv', sep=''), row.names = TRUE)
  write.csv(matt.trans, paste('Results/', can[1], '_RNA-corr-ALL.csv', sep=''), row.names = TRUE)
  matt.prot[is.na(matt.prot)] <- 0
  matt.trans[is.na(matt.trans)] <- 0
  
  annodf.prot = ensg_ccdprotein[ensg_ccdprotein$SYMBOL %in% row.names(matt.prot), c(3:6)]
  annodf.prot[annodf.prot=='N/A'] = NA
  annodf.prot <- sapply(annodf.prot, as.factor)
  rownames(annodf.prot) = ensg_ccdprotein[ensg_ccdprotein$SYMBOL %in% row.names(matt.prot),]$SYMBOL
  annodf.prot = data.frame(annodf.prot)
  
  Reason_color = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(annodf.prot$Reason)))
  names(Reason_color) <- levels(annodf.prot$Reason)
  Prot_Pseudotime_color = colorRampPalette(brewer.pal(9,"Purples"))(length(levels(annodf.prot$Prot_Pseudotime)))
  names(Prot_Pseudotime_color) <- levels(annodf.prot$Prot_Pseudotime)[order(levels(annodf.prot$Prot_Pseudotime))]
  Location_color = colorRampPalette(brewer.pal(8,"Set2"))(length(levels(annodf.prot$Location)))
  names(Location_color) <- levels(annodf.prot$Location)
  RNA_Pseudotime_color = colorRampPalette(brewer.pal(9,"Greens"))(length(levels(annodf.prot$RNA_Pseudotime)))
  names(RNA_Pseudotime_color) <- levels(annodf.prot$RNA_Pseudotime)[order(levels(annodf.prot$RNA_Pseudotime))]
  
  ann_colors_prot = list(
    Reason = Reason_color,
    Prot_Pseudotime = Prot_Pseudotime_color,
    Location = Location_color,
    RNA_Pseudotime = RNA_Pseudotime_color
  )
  
  annodf.trans = ensg_ccdprotein[ensg_ccdprotein$SYMBOL %in% row.names(matt.trans), c(3:6)]
  annodf.trans[annodf.trans=='N/A'] = NA
  annodf.trans <- sapply(annodf.trans, as.factor)
  rownames(annodf.trans) = ensg_ccdprotein[ensg_ccdprotein$SYMBOL %in% row.names(matt.trans),]$SYMBOL
  annodf.trans = data.frame(annodf.trans)
  
  Reason_color = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(annodf.trans$Reason)))
  names(Reason_color) <- levels(annodf.trans$Reason)
  Prot_Pseudotime_color = colorRampPalette(brewer.pal(9,"Purples"))(length(levels(annodf.trans$Prot_Pseudotime)))
  names(Prot_Pseudotime_color) <- levels(annodf.trans$Prot_Pseudotime)[order(levels(annodf.trans$Prot_Pseudotime))]
  Location_color = colorRampPalette(brewer.pal(8,"Set2"))(length(levels(annodf.trans$Location)))
  names(Location_color) <- levels(annodf.trans$Location)
  RNA_Pseudotime_color = colorRampPalette(brewer.pal(9,"Greens"))(length(levels(annodf.trans$RNA_Pseudotime)))
  names(RNA_Pseudotime_color) <- levels(annodf.trans$RNA_Pseudotime)[order(levels(annodf.trans$RNA_Pseudotime))]
  
  ann_colors_trans = list(
    Reason = Reason_color,
    Prot_Pseudotime = Prot_Pseudotime_color,
    Location = Location_color,
    RNA_Pseudotime = RNA_Pseudotime_color
  )
  
  pp = pheatmap(matt.prot, annotation_row = annodf.prot, annotation_col = annodf.prot, annotation_colors = ann_colors_prot, annotation_legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, main=paste(can[1], ' proteomics Spearman correlation ALL', sep=''))
  pt = pheatmap(matt.trans, annotation_row = annodf.trans, annotation_col =annodf.trans, annotation_colors = ann_colors_trans, annotation_legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, main=paste(can[1], ' transcriptomics Spearman correlation ALL', sep=''))
  pdf(file=paste('Results/', can[1], "_ALL-internal-corr.pdf", sep=''),
      width=20,height=11)
  grid.arrange(pp$gtable, pt$gtable, nrow=1, ncol=2)
  dev.off()
  
  matt.prot = read_csv(paste('Results/', can[1], '_protein-corr-ALL.csv', sep=''))
  ppg = cbind(matt.prot, clusters = cutree(pp$tree_row, k=can[2]))
  ppg = ppg[,c(ncol(ppg), 2:(ncol(ppg)-1))]
  write.csv(ppg, paste('Results/', can[1], "_prot_groups.csv", sep=''))
  matt.trans = read_csv(paste('Results/', can[1], '_RNA-corr-ALL.csv', sep=''))
  ptg = cbind(matt.trans, clusters = cutree(pt$tree_row, k=can[3]))
  ptg = ptg[,c(ncol(ptg), 2:(ncol(ptg)-1))]
  write.csv(ptg, paste('Results/', can[1], "_RNA_groups.csv", sep=''))
  

  prot.tumor = prot[!grepl("\\.N", prot$Patient_ID), ]
  trans.tumor = trans[!grepl("\\.N", trans$Patient_ID), ]
  prot.tumor = prot[!grepl("\\.C", prot$Patient_ID), ]
  trans.tumor = trans[!grepl("\\.C", trans$Patient_ID), ]
  prot.tumor = prot.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein_nontranscript_regulated$SYMBOL, colnames(prot.tumor)))]
  trans.tumor = trans.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein_nontranscript_regulated$SYMBOL, colnames(trans.tumor)))]
  matt.prot = cor(prot.tumor[, 2:ncol(prot.tumor)],  method = "spearman", use="pairwise.complete.obs")
  matt.trans = cor(trans.tumor[, 2:ncol(trans.tumor)],  method = "spearman", use="pairwise.complete.obs")
  write.csv(matt.prot, paste('Results/', can[1], '_protein-corr-nontranscript.csv', sep=''), row.names = TRUE)
  write.csv(matt.trans, paste('Results/', can[1], '_RNA-corr-nontranscript.csv', sep=''), row.names = TRUE)
  matt.prot[is.na(matt.prot)] <- 0
  matt.trans[is.na(matt.trans)] <- 0
  
  annodf.prot = ensg_ccdprotein_nontranscript_regulated[ensg_ccdprotein_nontranscript_regulated$SYMBOL %in% row.names(matt.prot), c(3:6)]
  annodf.prot[annodf.prot=='N/A'] = NA
  annodf.prot <- sapply(annodf.prot, as.factor)
  rownames(annodf.prot) = ensg_ccdprotein_nontranscript_regulated[ensg_ccdprotein_nontranscript_regulated$SYMBOL %in% row.names(matt.prot),]$SYMBOL
  annodf.prot = data.frame(annodf.prot)
  
  
  Reason_color = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(annodf.prot$Reason)))
  names(Reason_color) <- levels(annodf.prot$Reason)
  Prot_Pseudotime_color = colorRampPalette(brewer.pal(9,"Purples"))(length(levels(annodf.prot$Prot_Pseudotime)))
  names(Prot_Pseudotime_color) <- levels(annodf.prot$Prot_Pseudotime)[order(levels(annodf.prot$Prot_Pseudotime))]
  Location_color = colorRampPalette(brewer.pal(8,"Set2"))(length(levels(annodf.prot$Location)))
  names(Location_color) <- levels(annodf.prot$Location)
  RNA_Pseudotime_color = colorRampPalette(brewer.pal(9,"Greens"))(length(levels(annodf.prot$RNA_Pseudotime)))
  names(RNA_Pseudotime_color) <- levels(annodf.prot$RNA_Pseudotime)[order(levels(annodf.prot$RNA_Pseudotime))]
  
  ann_colors_prot = list(
    Reason = Reason_color,
    Prot_Pseudotime = Prot_Pseudotime_color,
    Location = Location_color,
    RNA_Pseudotime = RNA_Pseudotime_color
  )
  
  annodf.trans = ensg_ccdprotein_nontranscript_regulated[ensg_ccdprotein_nontranscript_regulated$SYMBOL %in% row.names(matt.trans), c(3:6)]
  annodf.trans[annodf.trans=='N/A'] = NA
  annodf.trans <- sapply(annodf.trans, as.factor)
  rownames(annodf.trans) = ensg_ccdprotein_nontranscript_regulated[ensg_ccdprotein_nontranscript_regulated$SYMBOL %in% row.names(matt.trans),]$SYMBOL
  annodf.trans = data.frame(annodf.trans)
  
  Reason_color = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(annodf.trans$Reason)))
  names(Reason_color) <- levels(annodf.trans$Reason)
  Prot_Pseudotime_color = colorRampPalette(brewer.pal(9,"Purples"))(length(levels(annodf.trans$Prot_Pseudotime)))
  names(Prot_Pseudotime_color) <- levels(annodf.trans$Prot_Pseudotime)[order(levels(annodf.trans$Prot_Pseudotime))]
  Location_color = colorRampPalette(brewer.pal(8,"Set2"))(length(levels(annodf.trans$Location)))
  names(Location_color) <- levels(annodf.trans$Location)
  RNA_Pseudotime_color = colorRampPalette(brewer.pal(9,"Greens"))(length(levels(annodf.trans$RNA_Pseudotime)))
  names(RNA_Pseudotime_color) <- levels(annodf.trans$RNA_Pseudotime)[order(levels(annodf.trans$RNA_Pseudotime))]
  
  ann_colors_trans = list(
    Reason = Reason_color,
    Prot_Pseudotime = Prot_Pseudotime_color,
    Location = Location_color,
    RNA_Pseudotime = RNA_Pseudotime_color
  )
  
  pp = pheatmap(matt.prot, annotation_row = annodf.prot, annotation_col = annodf.prot, annotation_colors = ann_colors_prot, annotation_legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, main=paste(can[1], ' proteomics Spearman correlation nontranscript', sep=''))
  pt = pheatmap(matt.trans, annotation_row = annodf.trans, annotation_col = annodf.trans, annotation_colors = ann_colors_trans, annotation_legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, main=paste(can[1], ' transcriptomics Spearman correlation nontranscript', sep=''))
  pdf(file=paste('Results/', can[1], "_nontranscript-internal-corr.pdf", sep=''),
      width=20,height=11)
  grid.arrange(pp$gtable, pt$gtable, nrow=1, ncol=2)
  dev.off()

  
  prot.tumor = prot[!grepl("\\.N", prot$Patient_ID), ]
  trans.tumor = trans[!grepl("\\.N", trans$Patient_ID), ]
  prot.tumor = prot[!grepl("\\.C", prot$Patient_ID), ]
  trans.tumor = trans[!grepl("\\.C", trans$Patient_ID), ]
  prot.tumor = prot.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein_transcript_regulated$SYMBOL, colnames(prot.tumor)))]
  trans.tumor = trans.tumor[ ,c("Patient_ID", intersect(ensg_ccdprotein_transcript_regulated$SYMBOL, colnames(trans.tumor)))]
  matt.prot = cor(prot.tumor[, 2:ncol(prot.tumor)],  method = "spearman", use="pairwise.complete.obs")
  matt.trans = cor(trans.tumor[, 2:ncol(trans.tumor)],  method = "spearman", use="pairwise.complete.obs")
  write.csv(matt.prot, paste('Results/', can[1], '_protein-corr-transcript.csv', sep=''), row.names = TRUE)
  write.csv(matt.trans, paste('Results/', can[1], '_RNA-corr-transcript.csv', sep=''), row.names = TRUE)
  matt.prot[is.na(matt.prot)] <- 0
  matt.trans[is.na(matt.trans)] <- 0
  
  annodf.prot = ensg_ccdprotein_transcript_regulated[ensg_ccdprotein_transcript_regulated$SYMBOL %in% row.names(matt.prot), c(3:6)]
  annodf.prot[annodf.prot=='N/A'] = NA
  annodf.prot <- sapply(annodf.prot, as.factor)
  rownames(annodf.prot) = ensg_ccdprotein_transcript_regulated[ensg_ccdprotein_transcript_regulated$SYMBOL %in% row.names(matt.prot),]$SYMBOL
  annodf.prot = data.frame(annodf.prot)
  
  Reason_color = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(annodf.prot$Reason)))
  names(Reason_color) <- levels(annodf.prot$Reason)
  Prot_Pseudotime_color = colorRampPalette(brewer.pal(9,"Purples"))(length(levels(annodf.prot$Prot_Pseudotime)))
  names(Prot_Pseudotime_color) <- levels(annodf.prot$Prot_Pseudotime)[order(levels(annodf.prot$Prot_Pseudotime))]
  Location_color = colorRampPalette(brewer.pal(8,"Set2"))(length(levels(annodf.prot$Location)))
  names(Location_color) <- levels(annodf.prot$Location)
  RNA_Pseudotime_color = colorRampPalette(brewer.pal(9,"Greens"))(length(levels(annodf.prot$RNA_Pseudotime)))
  names(RNA_Pseudotime_color) <- levels(annodf.prot$RNA_Pseudotime)[order(levels(annodf.prot$RNA_Pseudotime))]
  
  ann_colors_prot = list(
    Reason = Reason_color,
    Prot_Pseudotime = Prot_Pseudotime_color,
    Location = Location_color,
    RNA_Pseudotime = RNA_Pseudotime_color
  )
  
  annodf.trans = ensg_ccdprotein_transcript_regulated[ensg_ccdprotein_transcript_regulated$SYMBOL %in% row.names(matt.trans), c(3:6)]
  annodf.trans[annodf.trans=='N/A'] = NA
  annodf.trans <- sapply(annodf.trans, as.factor)
  rownames(annodf.trans) = ensg_ccdprotein_transcript_regulated[ensg_ccdprotein_transcript_regulated$SYMBOL %in% row.names(matt.trans),]$SYMBOL
  annodf.trans = data.frame(annodf.trans)
  
  Reason_color = colorRampPalette(brewer.pal(8,"Dark2"))(length(levels(annodf.trans$Reason)))
  names(Reason_color) <- levels(annodf.trans$Reason)
  Prot_Pseudotime_color = colorRampPalette(brewer.pal(9,"Purples"))(length(levels(annodf.trans$Prot_Pseudotime)))
  names(Prot_Pseudotime_color) <- levels(annodf.trans$Prot_Pseudotime)[order(levels(annodf.trans$Prot_Pseudotime))]
  Location_color = colorRampPalette(brewer.pal(8,"Set2"))(length(levels(annodf.trans$Location)))
  names(Location_color) <- levels(annodf.trans$Location)
  RNA_Pseudotime_color = colorRampPalette(brewer.pal(9,"Greens"))(length(levels(annodf.trans$RNA_Pseudotime)))
  names(RNA_Pseudotime_color) <- levels(annodf.trans$RNA_Pseudotime)[order(levels(annodf.trans$RNA_Pseudotime))]
  
  ann_colors_trans = list(
    Reason = Reason_color,
    Prot_Pseudotime = Prot_Pseudotime_color,
    Location = Location_color,
    RNA_Pseudotime = RNA_Pseudotime_color
  )
  

  pp = pheatmap(matt.prot, annotation_row = annodf.prot, annotation_col = annodf.prot, annotation_colors = ann_colors_prot, annotation_legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, main=paste(can[1], ' proteomics Spearman correlation transcript', sep=''))
  pt = pheatmap(matt.trans, annotation_row = annodf.trans, annotation_col = annodf.trans, annotation_colors = ann_colors_trans, annotation_legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, main=paste(can[1], ' transcriptomics Spearman correlation transcript', sep=''))
  pdf(file=paste('Results/', can[1], "_transcript-internal-corr.pdf", sep=''),
      width=20,height=11)
  grid.arrange(pp$gtable, pt$gtable, nrow=1, ncol=2)
  dev.off()
}

