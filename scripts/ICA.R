# ICA preparation
library(readr)
library(fastDummies)
library(dplyr)

exclude = c('Sample_ID', 'Patient_ID', 'Sample_Tumor_Normal', 'Sample.IDs', 'Human.Readable.Label', 'Replicate_Measurement_IDs', 'specimen_id', 'slide_id', 'Experiment', 'Channel')

cancers = c('BRCA', 'LUAD', 'LSCC', 'OV', 'CO', 'HN', 'UCEC', 'GBM', 'CCRCC')
for (can in cancers){
  clinical = read_csv(paste(can, "_clinical.csv", sep=''))
  prot = read_csv(paste(can, "_proteomics.csv", sep=''))
  prot = prot[rowSums(is.na(prot)) < ncol(prot)*0.3, ]
  clinical = clinical[clinical$Patient_ID %in% prot$Patient_ID,]
  
  clinical = clinical[clinical$Sample_Tumor_Normal== 'Tumor',]
  clinical = clinical[ , colSums(is.na(clinical)) == 0]
  rn = clinical$Patient_ID
  clinical = clinical[, !colnames(clinical) %in% exclude]
  indx <- sapply(clinical, is.factor)
  clinical[indx] <- lapply(clinical[indx], function(x) as.numeric(as.character(x)))
  
  clinical = dummy_cols(clinical)
  nums <- unlist(lapply(clinical, is.numeric))
  clinical = clinical[, nums]
  
  row.names(clinical) = rn

  prot = prot[prot$Patient_ID %in% rownames(clinical), ]
  
  for(i in 2:ncol(prot)){
    prot[is.na(prot[,i]), i] <- mean(unlist(prot[,i]), na.rm = TRUE)
  }
  
  write.csv(clinical, paste('ICA_', can, '_clinical.csv', sep=''), row.names=TRUE)
  write.csv(prot, paste('ICA_', can, '_proteomics.csv', sep=''), row.names=FALSE)
  
}


