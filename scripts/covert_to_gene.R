BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)
library(readr)

# ensg_ccdprotein <- read_csv("ensg_ccdprotein.csv", col_names = FALSE)
# ensg_ccdprotein_transcript_regulated <- read_csv("ensg_ccdprotein_transcript_regulated.csv", col_names = FALSE)
# ensg_ccdprotein_nontranscript_regulated <- read_csv("ensg_ccdprotein_nontranscript_regulated.csv", col_names = FALSE)

# colnames(ensg_ccdprotein) = c('GENEID')
# colnames(ensg_ccdprotein_transcript_regulated) = c('GENEID')
# colnames(ensg_ccdprotein_nontranscript_regulated) = c('GENEID')
# 
# ensembl.genes = ensg_ccdprotein$GENEID
# geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
# output = merge(ensg_ccdprotein, geneIDs1, by='GENEID')
# write.csv(output, 'ensg_ccdprotein.csv', row.names=FALSE)
# 
# ensembl.genes = ensg_ccdprotein_transcript_regulated$GENEID
# geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
# output = merge(ensg_ccdprotein_transcript_regulated, geneIDs1, by='GENEID')
# write.csv(output, 'ensg_ccdprotein_transcript_regulated.csv', row.names=FALSE)
# 
# ensembl.genes = ensg_ccdprotein_nontranscript_regulated$GENEID
# geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
# output = merge(ensg_ccdprotein_nontranscript_regulated, geneIDs1, by='GENEID')
# write.csv(output, 'ensg_ccdprotein_nontranscript_regulated.csv', row.names=FALSE)


library(readr)
ensg_ccdprotein <- read_csv("ensg_ccdprotein.csv")
ensg_ccdprotein_transcript_regulated <- read_csv("ensg_ccdprotein_transcript_regulated.csv")
ensg_ccdprotein_nontranscript_regulated <- read_csv("ensg_ccdprotein_nontranscript_regulated.csv")

CC = c('CCND1', 'CCNE1', 'CCNA2', 'CCNB1')
intersect(CC, ensg_ccdprotein$SYMBOL)
intersect(CC, ensg_ccdprotein_transcript_regulated$SYMBOL)
intersect(CC, ensg_ccdprotein_nontranscript_regulated$SYMBOL)


