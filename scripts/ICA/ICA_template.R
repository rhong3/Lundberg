#!/usr/bin/env Rscript

# command line arguments: a config file contains all arguments to ICA_template, tab-delimited 
#                         columns with headers: argument, value, description
# Arguments in ICA_template:
## exp_prefix
## data_dir 
## out_dir
## n_runs
## n_components
## run_seed
## save_Rdata

# import all arguments from config file
args = commandArgs(trailingOnly=TRUE)
message(Sys.time())
config_file=args[1]
config = read.table(config_file, sep='\t',header=T, as.is = 1:3)
rownames(config)=config$argument

exp_prefix=config['exp_prefix','value']
message(paste('exp_prefix:',exp_prefix))

data_dir=config['data_dir','value']
message(paste('data_dir:',data_dir))

out_dir=config['out_dir','value']
message(paste('out_dir:',out_dir))

n_runs=as.numeric(config['n_runs','value'])
message(paste('n_runs:',n_runs))

n_components=as.numeric(config['n_components','value'])
message(paste('n_components:',n_components))

run_seed=as.numeric(config['run_seed','value'])
message(paste('run_seed:',run_seed))

save_Rdata=config['save_Rdata','value']
message(paste('save Rdata?',save_Rdata))


library(stringr)

# read pre-processed data: row names are genes and column names are samples

data_file=paste(data_dir,'/',exp_prefix,'.csv',sep='')
message(paste('data file:',data_file))

# all variables : 'dat','samples','genes','ica','dis','iccls',
#                 'centroid','mixing',repeats','tsne'


dat=read.table(file=data_file,sep=',',header=T,row.names = 1)
samples=colnames(dat)
ns=dim(dat)[2]
genes=rownames(dat)

if(n_components%in%c(NA,'',NaN)){
  n_components=ns
}

source('./ICA_Clusters_Functions.R')

message('Start extracting independent components...')
if(is.na(run_seed)){
  run_seed = 123456
}

ica=icaRepeats(dat,seed=run_seed,nrun=n_runs,nc=n_components)
message('PAM clustering of extracted components...')
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

