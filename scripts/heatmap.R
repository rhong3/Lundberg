# Heatmap of samples
library(readr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)

cancers = c('BRCA', 'LUAD', 'LSCC', 'OV', 'CO', 'HN', 'UCEC', 'GBM', 'CCRCC')
ensg_ccdprotein <- read_csv("ensg_ccdprotein.csv")
ensg_ccdprotein_transcript_regulated <- read_csv("ensg_ccdprotein_transcript_regulated.csv")
ensg_ccdprotein_nontranscript_regulated <- read_csv("ensg_ccdprotein_nontranscript_regulated.csv")

for (can in cancers){
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
  
  for(i in 1:ncol(prot.tumor)){
    prot.tumor[is.na(prot.tumor[,i]), i] <- mean(unlist(prot.tumor[,i]), na.rm = TRUE)
  }
  
  for(i in 1:ncol(trans.tumor)){
    trans.tumor[is.na(trans.tumor[,i]), i] <- mean(unlist(trans.tumor[,i]), na.rm = TRUE)
  }
  
  ppid = prot.tumor$Patient_ID
  tpid = trans.tumor$Patient_ID
  prot.tumor = prot.tumor[,-1]
  trans.tumor = trans.tumor[,-1]
  row.names(prot.tumor) = ppid
  row.names(trans.tumor) = tpid
  
  annodf.prot = ensg_ccdprotein[ensg_ccdprotein$SYMBOL %in% colnames(prot.tumor), c(3:6)]
  annodf.prot[annodf.prot=='N/A'] = NA
  annodf.prot <- sapply(annodf.prot, as.factor)
  rownames(annodf.prot) = ensg_ccdprotein[ensg_ccdprotein$SYMBOL %in% colnames(prot.tumor),]$SYMBOL
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
  
  annodf.trans = ensg_ccdprotein[ensg_ccdprotein$SYMBOL %in% colnames(trans.tumor), c(3:6)]
  annodf.trans[annodf.trans=='N/A'] = NA
  annodf.trans <- sapply(annodf.trans, as.factor)
  rownames(annodf.trans) = ensg_ccdprotein[ensg_ccdprotein$SYMBOL %in% colnames(trans.tumor),]$SYMBOL
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
  
  pp = pheatmap(prot.tumor,  annotation_col = annodf.prot, annotation_colors = ann_colors_prot, annotation_legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, main=paste(can[1], ' proteomics ALL', sep=''))
  pt = pheatmap(trans.tumor, annotation_col = annodf.trans, annotation_colors = ann_colors_trans, annotation_legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, main=paste(can[1], ' transcriptomics ALL', sep=''))
  pdf(file=paste('Results/', can[1], "_ALL.pdf", sep=''),
      width=30,height=20)
  grid.arrange(pp$gtable, pt$gtable, nrow=2, ncol=1)
  dev.off()
  
}


