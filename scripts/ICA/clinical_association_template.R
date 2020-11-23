#!/usr/bin/env Rscript

# Clinical association template
# Plot average mixing score matrix with clinical annotations
# Find association between clinical features and signatures extracted by ICA
# Clinical data should be pre-processed such that ordinals are converted in to numeric values,
# and binary variables are in 0s and 1s.

# Command line arguments: exp_prefix, clinical_data, ica_rdata, out_dir p_value


args = commandArgs(trailingOnly=TRUE)
message(Sys.time())

exp_prefix=args[1]
message(paste('exp_prefix:',exp_prefix))

clinical_data=args[2]
message(paste('clinical_data:',clinical_data))

ica_rdata=args[3]
message(paste('ica_rdata:',ica_rdata))

out_dir=args[4]
message(paste('out_dir:',out_dir))

p_value=args[5]
message(paste('p_value_threshold:',p_value))

# load required functions

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(stringr)

source('./ICA_Clusters_Functions.R')
# source('~/R_Functions/heatmap.3.R')

clinical=read.table(file=clinical_data,
                    header=T,sep=',',row.names = 1)
#dim(clinical)

#ica_rdata='/media/lwk/lab/OV/retro/data/ov_pr_prosp_ICA.Rdata'
temp_space=new.env()
ica_vars_loaded=load(ica_rdata,temp_space)
#exp_prefix='ov_pr_prosp'


ica_vars=sub(".*_", "", ica_vars_loaded)

for (i in 1:length(ica_vars)){
  if(gsub(ica_vars[i],ica_vars_loaded[i],replacement = '')!=paste(exp_prefix,'_',sep='')){
    stop('Experiment name and ICA results data do not match. Double check the input!')
  }
  assign(ica_vars[i],get(ica_vars_loaded[i],envir=temp_space))
}

# if (!all(rownames(clinical)==Clean_proteomics_samples)){
#   stop('Sample names in clinical data do not match ICA results. Double check the input!')
# }


#dim(ica$A)


## association test (linear regression)

icP=findCor2(ica$A,clustering=iccls$clustering,clinical.data=clinical)
icPave=apply(icP,2,function(x)tapply(x,iccls$clustering,function(x)sum(x<as.numeric(p_value),na.rm=T)))

write.table(data.frame(cluster = iccls$clustering,icP),
            file=paste(out_dir,'/',exp_prefix,'_IC_Clinical_Correlation_P_Value_raw_data.tsv',sep=''),
            sep='\t',row.names = T, col.names = NA)

write.table(icPave[apply(icPave,1,sum)>0,apply(icPave,2,sum)>0],
            file=paste(out_dir,'/',exp_prefix,'_IC_Clinical_Correlation_P_Value.tsv',sep=''),
            sep='\t',row.names = T, col.names = NA)

write.table(icPave,
            file=paste(out_dir,'/',exp_prefix,'_IC_Clinical_Correlation_P_Value_all.tsv',sep=''),
            sep='\t',row.names = T, col.names = NA)

# heatmap of significal clinical association counts
icPave_plot=melt(icPave)
icPave_plot=subset(icPave_plot,
                      X1%in%names(which(apply(icPave,1,sum)>0))&X2%in%colnames(icPave)[which(apply(icPave,2,sum)>0)])
colnames(icPave_plot)[1:2]=c('IC','clinical_feature')
icPave_plot[,1]=factor(icPave_plot[,1])
icPave_plot$IC=paste('IC',str_pad(icPave_plot$IC,2,pad='0'),sep='_')
png(filename=paste(out_dir,paste(exp_prefix,'IC_cluster_clinical_association.png',sep='_'),sep='/'),
    width = 1200, height = 480, units = "px", 
    pointsize = 12,
    bg = 'transparent')

ggplot(icPave_plot,aes(x=IC,y=clinical_feature))+
  theme_classic(base_size=12)+
  theme(axis.line=element_blank())+
  theme(axis.text.x = element_text(angle=90,vjust=0.5))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient2(low='white',high='purple')+
  geom_text(aes(label=round(value,digits=2)),size=2)
graphics.off()



