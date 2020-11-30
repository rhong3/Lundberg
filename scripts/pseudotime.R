library(readr)
library(EnsDb.Hsapiens.v79)

ensg_ccdprotein <- read_csv("ensg_ccdprotein.csv")
ensg_ccdprotein_transcript_regulated <- read_csv("ensg_ccdprotein_transcript_regulated.csv")
ensg_ccdprotein_nontranscript_regulated <- read_csv("ensg_ccdprotein_nontranscript_regulated.csv")
# pseudo time
CCDProteins <- read_csv("CCDProteins.csv")
CCDRNA <- read_csv("CCDTranscripts.csv")
CCDRNA = CCDRNA[, c(1,3)]

colnames(CCDProteins) = c('GENEID', 'Reason', 'Prot_Pseudotime', 'Location')
colnames(CCDRNA) = c('GENEID', 'RNA_Pseudotime')

metabolic = read_csv('20201126_MetabolicEnzymes_FenyoLab.csv')
metabolic = metabolic[,-2]
colnames(metabolic)[1] = 'GENEID' 

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CCDProteins$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
CCDProteins = merge(CCDProteins, geneIDs1, by='GENEID')
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= CCDRNA$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
CCDRNA = merge(CCDRNA, geneIDs1, by='GENEID')
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= metabolic$GENEID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
metabolic = merge(metabolic, geneIDs1, by='GENEID')

CCDProteins = CCDProteins[,-1]
CCDRNA = CCDRNA[,-1]
metabolic = metabolic[,-1]

mggCCD= merge(CCDProteins, CCDRNA, by='SYMBOL', all=TRUE)

mggCCD = merge(mggCCD, metabolic, by='SYMBOL', all=TRUE)

ensg_ccdprotein = merge(ensg_ccdprotein, mggCCD, by='SYMBOL')
ensg_ccdprotein_transcript_regulated = merge(ensg_ccdprotein_transcript_regulated, mggCCD, by='SYMBOL')
ensg_ccdprotein_nontranscript_regulated = merge(ensg_ccdprotein_nontranscript_regulated, mggCCD, by='SYMBOL')

write.csv(ensg_ccdprotein, 'ensg_ccdprotein.csv', row.names=FALSE)
write.csv(ensg_ccdprotein_transcript_regulated, 'ensg_ccdprotein_transcript_regulated.csv', row.names=FALSE)
write.csv(ensg_ccdprotein_nontranscript_regulated, 'ensg_ccdprotein_nontranscript_regulated.csv', row.names=FALSE)
write.csv(mggCCD, 'CCD_merged.csv', row.names=FALSE)


