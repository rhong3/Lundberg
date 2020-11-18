# Correlation matrices within proteomics and protein-RNA correlations
library(readr)
endometrial_proteomics <- read_csv("endometrial_proteomics.csv")
endometrial_transcriptomics <- read_csv("endometrial_transcriptomics.csv")
endometrial_proteomics[, 2:ncol(endometrial_proteomics)] = endometrial_proteomics[, order(colnames(endometrial_proteomics[, 2:ncol(endometrial_proteomics)]))]
endometrial_transcriptomics[, 2:ncol(endometrial_transcriptomics)] = endometrial_transcriptomics[, order(colnames(endometrial_transcriptomics[, 2:ncol(endometrial_transcriptomics)]))]

endometrial_proteomics[, 2:ncol(endometrial_proteomics)] <- sapply(endometrial_proteomics[, 2:ncol(endometrial_proteomics)] , as.numeric)
endometrial_transcriptomics[, 2:ncol(endometrial_transcriptomics)] <- sapply(endometrial_transcriptomics[, 2:ncol(endometrial_transcriptomics)] , as.numeric)

ensg_ccdprotein <- read_csv("ensg_ccdprotein.csv")
ensg_ccdprotein_transcript_regulated <- read_csv("ensg_ccdprotein_transcript_regulated.csv")
ensg_ccdprotein_nontranscript_regulated <- read_csv("ensg_ccdprotein_nontranscript_regulated.csv")

# 381
endometrial_proteomics = endometrial_proteomics[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(endometrial_proteomics)))]
pr = endometrial_proteomics$Patient_ID
# 524
endometrial_transcriptomics = endometrial_transcriptomics[ ,c("Patient_ID", intersect(ensg_ccdprotein$SYMBOL, colnames(endometrial_transcriptomics)))]
tr = endometrial_transcriptomics$Patient_ID

# 109 samples
endometrial_proteomics = endometrial_proteomics[intersect(rownames(endometrial_proteomics), rownames(endometrial_transcriptomics)), ]
row.names(endometrial_proteomics) = pr
endometrial_transcriptomics = endometrial_transcriptomics[intersect(rownames(endometrial_proteomics), rownames(endometrial_transcriptomics)), ]
row.names(endometrial_transcriptomics) = tr
endometrial_proteomics = endometrial_proteomics[order(endometrial_proteomics$Patient_ID),]
endometrial_transcriptomics = endometrial_transcriptomics[order(endometrial_transcriptomics$Patient_ID),]

# Tumor 95 samples
endometrial_proteomics.tumor = endometrial_proteomics[!grepl("\\.N", endometrial_proteomics$Patient_ID), ]
endometrial_transcriptomics.tumor = endometrial_transcriptomics[!grepl("\\.N", endometrial_transcriptomics$Patient_ID), ]

# prot-trans (95x382)
endometrial_proteomics.tumor.x = endometrial_proteomics.tumor[,intersect(colnames(endometrial_proteomics.tumor), colnames(endometrial_transcriptomics.tumor))]
endometrial_transcriptomics.tumor.x = endometrial_transcriptomics.tumor[,intersect(colnames(endometrial_proteomics.tumor), colnames(endometrial_transcriptomics.tumor))]

inter_cor = diag(cor(as.matrix(endometrial_proteomics.tumor.x[, 2:ncol(endometrial_proteomics.tumor.x)]), as.matrix(endometrial_transcriptomics.tumor.x[, 2:ncol(endometrial_transcriptomics.tumor.x)]),  method = "spearman", use="pairwise.complete.obs"))
inter_cor = data.frame(inter_cor)
inter_cor$regulated = row.names(inter_cor) %in% ensg_ccdprotein_transcript_regulated$SYMBOL
inter_cor = na.omit(inter_cor) #378
colnames(inter_cor) = c('correlation coefficient', 'transcript regulated')

library(ggplot2)
library(ggpubr)
library(gridExtra)
# Change histogram plot fill colors by groups
p = ggplot(inter_cor, aes(x=`correlation coefficient`, fill=`transcript regulated`, color=`transcript regulated`)) +
  geom_histogram(position="identity", alpha=0.5, binwidth=0.05)+scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  labs(title="UCEC protein and mRNA Spearman correlation")+ 
  theme_classic()+
  theme(legend.position="top", plot.title = element_text(hjust = 0.5)) 
pdf(file=paste("UCEC_protein-mRNA-corr.pdf", sep=''),
    width=10,height=8)
grid.arrange(p,nrow=1, ncol=1)
dev.off()

