
cancers = c('BRCA', 'LUAD', 'LSCC', 'OV', 'CO', 'HN', 'UCEC', 'GBM', 'CCRCC')
for (can in cancers){
  exp_prefix=can
  
  out_dir='Results'
  
  p_value=0.0005
  
  clinical_data=paste('ICA_', can, '_clinical.csv', sep='')
  
  ica_rdata=paste('Results/', can, '_ICA.Rdata', sep='')
  
  # load required functions
  
  library(gplots)
  library(RColorBrewer)
  library(ggplot2)
  library(gridExtra)
  library(stringr)
  
  source('scripts/ICA/ICA_Clusters_Functions.R')
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

  
  p = ggplot(icPave_plot,aes(x=IC,y=clinical_feature))+
    theme_classic(base_size=30)+
    theme(axis.line=element_blank())+
    theme(axis.text.x = element_text(angle=90,vjust=0.5))+
    geom_tile(aes(fill=value))+
    scale_fill_gradient2(low='white',high='purple')+
    geom_text(aes(label=round(value,digits=10)),size=10)
  
  pdf(file=paste(out_dir,paste(exp_prefix,'IC_cluster_clinical_association.pdf',sep='_'),sep='/'),
      width = 30, height = 25, 
      bg = 'transparent')
  grid.arrange(p,nrow=1, ncol=1)
  dev.off()
}


