## Outlier Analysis prep
#' 
#' Created on 11/27/2020
#' 
#' @author: RH

# BRCA subtypes
# proteomics
prot = read.csv("BRCA_proteomics.csv", row.names=1)
prot = prot[rowSums(is.na(prot)) < ncol(prot)*0.3, ]
for(i in 1:ncol(prot)){
  prot[is.na(prot[,i]), i] <- mean(unlist(prot[,i]), na.rm = TRUE)
}
clinical = read.csv('ICA_BRCA_clinical.csv', row.names=1)
clinical = clinical[, 1:5]
write.csv(clinical, "Results/OLA/BRCA_subtype/BRCA_subtype.csv")

fileConn<-file("Results/OLA/BRCA_subtype/samples.txt")
writeLines(rownames(clinical), fileConn)
close(fileConn)

fileConn<-file("Results/OLA/BRCA_subtype/Basal.txt")
writeLines(rownames(clinical[which(clinical$PAM50_Basal==1),]), fileConn)
close(fileConn)
fileConn<-file("Results/OLA/BRCA_subtype/non_Basal.txt")
writeLines(rownames(clinical[which(clinical$PAM50_Basal==0),]), fileConn)
close(fileConn)

fileConn<-file("Results/OLA/BRCA_subtype/Her2.txt")
writeLines(rownames(clinical[which(clinical$PAM50_Her2==1),]), fileConn)
close(fileConn)
fileConn<-file("Results/OLA/BRCA_subtype/non_Her2.txt")
writeLines(rownames(clinical[which(clinical$PAM50_Her2==0),]), fileConn)
close(fileConn)

fileConn<-file("Results/OLA/BRCA_subtype/LumA.txt")
writeLines(rownames(clinical[which(clinical$PAM50_LumA==1),]), fileConn)
close(fileConn)
fileConn<-file("Results/OLA/BRCA_subtype/non_LumA.txt")
writeLines(rownames(clinical[which(clinical$PAM50_LumA==0),]), fileConn)
close(fileConn)

fileConn<-file("Results/OLA/BRCA_subtype/LumB.txt")
writeLines(rownames(clinical[which(clinical$PAM50_LumB==1),]), fileConn)
close(fileConn)
fileConn<-file("Results/OLA/BRCA_subtype/non_LumB.txt")
writeLines(rownames(clinical[which(clinical$PAM50_LumB==0),]), fileConn)
close(fileConn)

fileConn<-file("Results/OLA/BRCA_subtype/Normal.txt")
writeLines(rownames(clinical[which(clinical$PAM50_Normal==1),]), fileConn)
close(fileConn)
fileConn<-file("Results/OLA/BRCA_subtype/non_Normal.txt")
writeLines(rownames(clinical[which(clinical$PAM50_Normal==0),]), fileConn)
close(fileConn)

prot = data.frame(t(prot))
prot$gene = rownames(prot)

prot.oladata = prot[,match(rownames(clinical), colnames(prot))]
prot.oladata = cbind(prot[123], prot.oladata)
write.csv(prot.oladata, "Results/OLA/BRCA_subtype/ola_data.csv")

