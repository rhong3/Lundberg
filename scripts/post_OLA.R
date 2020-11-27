# Post OLA overlap
# outlier summary
library(readr)

summ = data.frame(gene=character(), group=character())
for (i in c('Basal', 'non-Basal', 'LumA', 'non-LumA', 'LumB', 'Her2', 'Normal')){
  ola = read_csv(paste("Results/OLA/BRCA_subtype/", i, ".outlier_sites_in_", i, ".txt", sep=''), col_names = FALSE)
  ola$group = i
  colnames(ola) = c('gene', 'group')
  summ = rbind(summ, ola)
}

summ$CCD = 'No'
transcript = read_csv('ensg_ccdprotein_transcript_regulated.csv')
non_transcript = read_csv('ensg_ccdprotein_nontranscript_regulated.csv')

summ[which(summ$gene %in% transcript$SYMBOL), ]$CCD = 'Yes-regulated'
summ[which(summ$gene %in% non_transcript$SYMBOL), ]$CCD = 'Yes-non-regulated'

write_csv(summ, 'Results/OLA/BRCA_subtype/summary.csv')


