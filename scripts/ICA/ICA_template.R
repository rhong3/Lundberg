cancers = c('BRCA', 'LUAD', 'LSCC', 'OV', 'CO', 'HN', 'UCEC', 'GBM', 'CCRCC')

for (can in cancers){
  exp_prefix=can
  
  out_dir='Results'
  
  n_runs=10
  
  n_components=30
  
  run_seed=NA
  
  save_Rdata='y'
  
  library(stringr)
  
  # read pre-processed data: row names are genes and column names are samples
  
  data_file=paste('ICA_', can, '_proteomics.csv', sep='')
  
  # all variables : 'dat','samples','genes','ica','dis','iccls',
  #                 'centroid','mixing',repeats','tsne'
  
  
  dat=read.table(file=data_file,sep=',',header=T,row.names = 1)
  dat = t(dat)
  samples=colnames(dat)
  ns=dim(dat)[2]
  genes=rownames(dat)
  
  if(n_components%in%c(NA,'',NaN)){
    n_components=ns
  }
  
  source('scripts/ICA/ICA_Clusters_Functions.R')
  
  message('Start extracting independent components...')
  if(is.na(run_seed)){
    run_seed = 123456
  }
  
  ica=icaRepeats(dat,seed=run_seed,nrun=n_runs,nc=n_components)
  message('clustering of extracted components...')
  dis=dist_cal(ics=ica$S)
  
  iccls=icPam(dis,nc=n_components)
  centroid=t(apply(ica$S,1,function(x)tapply(x,iccls$clustering,mean)))
  
  centroid_file=paste(out_dir,'/',exp_prefix,'_IC_centroid.txt',sep='')
  write.table(centroid,
              centroid_file,
              sep='\t',col.names = NA)
  message(paste('signatures saved as:',centroid_file))
  
  mixing = apply(ica$A,2,function(x)tapply(x,as.factor(iccls$clustering),mean))
  rownames(mixing)=paste('IC',str_pad(rownames(mixing),3,pad='0'),sep='_')
  mixing_file=paste(out_dir,'/',exp_prefix,'_IC_mean_mixing_score.txt',sep='')
  write.table(mixing,
              mixing_file,
              sep='\t',col.names = NA)
  message(paste('mean mixing scores saved as:',mixing_file))
  
  
  repeats=tapply(ica$run.ind,iccls$clustering,function(x)length(unique(x)))
  
  png(filename=paste(out_dir,paste(exp_prefix,'cluster_consistency.png',sep='_'),sep='/'),
      width = 480, height = 480, units = "px", 
      pointsize = 12,
      bg = 'transparent')
  
  plot(repeats,iccls$silinfo$clus.avg.widths,
       pch=16,col='skyblue',cex=2,
       xlab='present in different runs',ylab='cluster silhouette', main='cluster_consistency_assessment')
  text(1:ns,x=repeats,y=iccls$silinfo$clus.avg.widths)
  
  graphics.off()
  
  message('Cluster consistency plot generated')
  
  library(tsne)
  message('Calculating tsne representation of all ICs')
  tsne=tsne(dis)
  tsne=data.frame(tsne)
  tsne=cbind(tsne,X3=iccls$clustering)
  
  colnames(tsne)=c('X1','X2','X3')
  
  png(filename=paste(out_dir,paste(exp_prefix,'signature_tsne.png',sep='_'),sep='/'),
      width = 480, height = 480, units = "px", 
      pointsize = 12,
      bg = 'transparent')
  cls2d(tsne,title=paste(exp_prefix,'tsne',sep='_'))
  graphics.off()
  message('TSNE plot generated')
  
  
  if (save_Rdata=='y'){
    var_list=c('dat','samples','genes','ica','dis','iccls','centroid','mixing','repeats','tsne')
    var_list_save=paste(exp_prefix,var_list,sep='_')
    n=length(var_list)
    for (i in 1:n){
      assign(var_list_save[i],get(var_list[i]))
    }
    outfile=paste(out_dir,'/',exp_prefix,'_','ICA.Rdata',sep='')
    save(list=var_list_save,
         file=outfile)
    message(paste('Save variables as:',outfile))
  }
}

