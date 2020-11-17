## Outlier Analysis prep
#' 
#' Created on 1/15/2020
#' 
#' @author: RH

# 5-yr survival
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$survival_5yr),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/5-yr-survival/5yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/5-yr-survival/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/5-yr-survival/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$survival_5yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/5-yr-survival/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$survival_5yr==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/5-yr-survival/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$survival_5yr),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/5-yr-survival/5yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/5-yr-survival/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/5-yr-survival/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_5yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/5-yr-survival/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_5yr==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/5-yr-survival/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$survival_5yr),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/5-yr-survival/5yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/5-yr-survival/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/5-yr-survival/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$survival_5yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/5-yr-survival/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$survival_5yr==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/5-yr-survival/ola_data.csv")


# Gender
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$gender),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/gender/gender-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/gender/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/gender/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$gender==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/gender/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$gender==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/gender/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$gender),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/gender/gender-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/gender/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/gender/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$gender==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/gender/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$gender==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/gender/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$gender),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/gender/gender-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/gender/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/gender/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$gender==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/gender/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$gender==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/gender/ola_data.csv")


# Survival 1 year
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$survival_1yr),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/1-yr-survival/1yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/1-yr-survival/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/1-yr-survival/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$survival_1yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/1-yr-survival/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$survival_1yr==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/1-yr-survival/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$survival_1yr),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/1-yr-survival/1yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/1-yr-survival/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/1-yr-survival/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_1yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/1-yr-survival/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_1yr==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/1-yr-survival/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$survival_1yr),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/1-yr-survival/1yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/1-yr-survival/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/1-yr-survival/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$survival_1yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/1-yr-survival/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$survival_1yr==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/1-yr-survival/ola_data.csv")

# Survival 6 month
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$survival_6mo),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/6-month-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$survival_6mo==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$survival_6mo==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/6-month-survival/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$survival_6mo),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/6-month-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_6mo==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_6mo==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/6-month-survival/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$survival_6mo),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/6-month-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$survival_6mo==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$survival_6mo==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/6-month-survival/ola_data.csv")

# BRAF
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$BRAF),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/BRAF/BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/BRAF/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/BRAF/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/BRAF/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$BRAF==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/BRAF/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$BRAF),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/BRAF/BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/BRAF/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/BRAF/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/BRAF/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$BRAF==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/BRAF/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$BRAF),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/BRAF/BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/BRAF/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/BRAF/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/BRAF/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$BRAF==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/BRAF/ola_data.csv")

# NRAS
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$NRAS),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/NRAS/NRAS-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/NRAS/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/NRAS/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$NRAS==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/NRAS/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$NRAS==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/NRAS/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$NRAS),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/NRAS/NRAS-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/NRAS/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/NRAS/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$NRAS==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/NRAS/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$NRAS==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/NRAS/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$NRAS),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/NRAS/NRAS-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/NRAS/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/NRAS/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$NRAS==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/NRAS/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$NRAS==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/NRAS/ola_data.csv")


## stage
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[prot$dis.stage>2,]
prot.ola = prot[!is.na(prot$dis.stage),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/stage/stage-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/stage/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/stage/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$dis.stage==3),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/stage/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$dis.stage==4),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/stage/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[trans$dis.stage>2,]
trans.ola = trans[!is.na(trans$dis.stage),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/stage/stage-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/stage/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/stage/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$dis.stage==3),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/stage/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$dis.stage==4),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/stage/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[phospho$dis.stage>2,]
phospho.ola = phospho[!is.na(phospho$dis.stage),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/stage/stage-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/stage/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/stage/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$dis.stage==3),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/stage/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$dis.stage==4),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/stage/ola_data.csv")


## primary to metastasis 400
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$dmfs.days),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM400/PtoM400-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM400/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM400/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$dmfs.days>400),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM400/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$dmfs.days<=400),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM400/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$dmfs.days),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM400/PtoM400-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM400/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM400/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$dmfs.days>400),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM400/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$dmfs.days<=400),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM400/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$dmfs.days),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM400/PtoM400-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM400/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM400/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$dmfs.days>400),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM400/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$dmfs.days<=400),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM400/ola_data.csv")


## primary to metastasis 600
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$dmfs.days),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/PtoM600-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$dmfs.days>600),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$dmfs.days<=600),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/PtoM600/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$dmfs.days),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM600/PtoM600-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM600/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM600/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$dmfs.days>600),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM600/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$dmfs.days<=600),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/PtoM600/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$dmfs.days),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM600/PtoM600-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM600/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM600/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$dmfs.days>600),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM600/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$dmfs.days<=600),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/PtoM600/ola_data.csv")

# Survival 3 year
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$survival_3yr),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-survival/3yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-survival/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-survival/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$survival_3yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-survival/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$survival_3yr==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-survival/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$survival_3yr),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-survival/3yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-survival/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-survival/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_3yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-survival/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$survival_3yr==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-survival/ola_data.csv")


# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$survival_3yr),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-survival/3yr-survival-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-survival/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-survival/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$survival_3yr==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-survival/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$survival_3yr==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-survival/ola_data.csv")


# Survival 3 year alive BRAF
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot = prot[prot$os.days >= 1095,]
prot = prot[!is.na(prot$os.days),]
prot = prot[!is.na(prot$os.events),]
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$BRAF),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-alive-BRAF/3-yr-alive-BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-alive-BRAF/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-alive-BRAF/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-alive-BRAF/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$BRAF==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-alive-BRAF/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans = trans[trans$os.days >= 1095,]
trans = trans[!is.na(trans$os.days),]
trans = trans[!is.na(trans$os.events),]
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$BRAF),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-alive-BRAF/3-yr-alive-BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-alive-BRAF/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-alive-BRAF/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-alive-BRAF/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$BRAF==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-alive-BRAF/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho = phospho[phospho$os.days >= 1095,]
phospho = phospho[!is.na(phospho$os.days),]
phospho = phospho[!is.na(phospho$os.events),]
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$BRAF),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-alive-BRAF/3-yr-alive-BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-alive-BRAF/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-alive-BRAF/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-alive-BRAF/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$BRAF==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-alive-BRAF/ola_data.csv")


# Survival 3 year dead BRAF
# proteomics
prot = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics_clinical.csv", row.names=1)
prot = prot[prot$os.events == 1 & prot$os.days < 1095,]
prot = prot[!is.na(prot$os.days),]
prot = prot[!is.na(prot$os.events),]
prot.oridata = read.csv("~/documents/Segundo_Melanoma/Data/proteomics/proteomics.csv")
prot.ola = prot[!is.na(prot$BRAF),]
write.csv(prot.ola, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-dead-BRAF/3-yr-dead-BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-dead-BRAF/samples.txt")
writeLines(rownames(prot.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-dead-BRAF/G1.txt")
writeLines(rownames(prot.ola[which(prot.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-dead-BRAF/G2.txt")
writeLines(rownames(prot.ola[which(prot.ola$BRAF==0),]), fileConn)
close(fileConn)
prot.oladata = prot.oridata[,match(rownames(prot.ola), colnames(prot.oridata))]
prot.oladata = cbind(prot.oridata[111:113], prot.oladata)
write.csv(prot.oladata, "~/documents/Segundo_Melanoma/Results/proteomics/OLA/3-yr-dead-BRAF/ola_data.csv")

# transcriptomics
trans = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics_clinical.csv", row.names=1)
trans = trans[trans$os.events == 1 & trans$os.days < 1095,]
trans = trans[!is.na(trans$os.days),]
trans = trans[!is.na(trans$os.events),]
trans.oridata = read.csv("~/documents/Segundo_Melanoma/Data/transcriptomics/transcriptomics.csv")
trans.ola = trans[!is.na(trans$BRAF),]
write.csv(trans.ola, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-dead-BRAF/3-yr-dead-BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-dead-BRAF/samples.txt")
writeLines(rownames(trans.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-dead-BRAF/G1.txt")
writeLines(rownames(trans.ola[which(trans.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-dead-BRAF/G2.txt")
writeLines(rownames(trans.ola[which(trans.ola$BRAF==0),]), fileConn)
close(fileConn)
trans.oladata = trans.oridata[,match(rownames(trans.ola), colnames(trans.oridata))]
trans.oladata = cbind(trans.oridata[1:2], trans.oladata)
write.csv(trans.oladata, "~/documents/Segundo_Melanoma/Results/transcriptomics/OLA/3-yr-dead-BRAF/ola_data.csv")

# phospho
phospho = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho_clinical.csv", row.names=1)
phospho = phospho[phospho$os.events == 1 & phospho$os.days < 1095,]
phospho = phospho[!is.na(phospho$os.days),]
phospho = phospho[!is.na(phospho$os.events),]
phospho.oridata = read.csv("~/documents/Segundo_Melanoma/Data/phospho/phospho.csv")
phospho.ola = phospho[!is.na(phospho$BRAF),]
write.csv(phospho.ola, "~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-dead-BRAF/3-yr-dead-BRAF-clinical.csv")
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-dead-BRAF/samples.txt")
writeLines(rownames(phospho.ola), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-dead-BRAF/G1.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$BRAF==1),]), fileConn)
close(fileConn)
fileConn<-file("~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-dead-BRAF/G2.txt")
writeLines(rownames(phospho.ola[which(phospho.ola$BRAF==0),]), fileConn)
close(fileConn)
phospho.oladata = phospho.oridata[,match(rownames(phospho.ola), colnames(phospho.oridata))]
phospho.oladata = cbind(phospho.oridata[1:4], phospho.oladata)
write.csv(phospho.oladata, "~/documents/Segundo_Melanoma/Results/phospho/OLA/3-yr-dead-BRAF/ola_data.csv")
