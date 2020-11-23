library(fastICA)
library(ggplot2)
library(tsne)
library(cluster)
library(reshape)

##ICAClusters functions
icaRepeats=function(data,seed,nrun,nc=NULL){
    
    require(fastICA)
    message(paste('Random seed:', seed))
    message(paste('Number of runs:', nrun))
    message(paste('Dimension of component:',dim(data)[1]))
    
    if (is.null(nc)){nc=dim(data)[2]}
    message(paste('Number of components extracted:', nc))
    
    set.seed(seed)
    run.seed=runif(nrun,min=0,max=10e6)
    
    #SList=list()
    #AList=list()
    allS=data.frame(row.names=rownames(data))
    allA=data.frame()
    
    for (i in 1:nrun){
        
        initM=matrix(nrow=nc,ncol=nc)
        set.seed(run.seed[i])
        init.seed=runif(nc,min=0,max=10e6)
        for (j in 1:nc){
            set.seed(init.seed[j])
            initM[,j]=rnorm(nc)
        }
        #print(dim(initM))
        
        icaRes=fastICA(data,n.comp=nc,fun='logcosh',method='C',w.init=initM)
        Sn=icaRes$S
        An=icaRes$A
        
        #print(Sn)
        for (k in 1:nc){
            if(sum(Sn[,k]>3)<sum(Sn[,k]<(-3))){
                Sn[,k]=-Sn[,k]
                An[k,]=-An[k,]
            }
        }
        #SList[[i]]=Sn
        #AList[[i]]=An
        allS=cbind(allS,Sn)
        allA=rbind(allA,An)
        message(paste('Run',i, 'completed'))
    }
    
    #ic=do.call(cbind,SList)
    #mix=do.call(rbind,AList)
    run.ind=rep(1:nrun,each=nc)
    ic.ind=rep(1:nc,nrun)
    
    colnames(allA)=colnames(data)
    
    
    return(list(run.ind=run.ind,ic.ind=ic.ind,S=allS,A=allA))
    
}

dist_cal=function(ics,dist_method='spearman'){
  if (dist_method=='spearman'){
    ## cor function computes correlation between columns
    d=(1-cor(ics,method='spearman'))/2
    d=as.dist(d)
  }
  if(dist_method=='euclidean'){
    ## dist function computes distance between rows
    d=dist(t(ics),method='euclidean')
  }
  
  return(d)
}

icPam=function(d,nc=10){
    require(cluster)
    cls=pam(d,k=nc)
    return(cls)
}


icKmeans=function(d,nc=10,nstart=5){
  cls=kmeans(d,k=nc,nstart=nstart)
  return(cls)
}

icDbscan=function(d,eps,minPts){
  require(dbscan)
  cls=dbscan(d,eps=eps,minPts = minPts)
  return(cls)
}

cls2d=function(data,title=NULL){
    plot(data$X1,data$X2,xlab='X1',ylab='X2',
         #xlim=c(-10,10),ylim=c(-10,10),
         pch=16,
         col=alpha(data$X3,0.2))
    text(unique(data$X3),
         x=tapply(data$X1,data$X3,mean),
         y=tapply(data$X2,data$X3,mean),col='black')
    if(!is.null(title)){
      title(title)
    }
}

#cls3d=function(data){
#    library(scatterplot3d)
#    library(ggplot2)
    
#    s3d=with(data,scatterplot3d(X1,X2,X3,
#                                    box=F,
#                                    pch=16,
#                                    color=alpha(X4,0.2)
#    ))
    
#    s3d.coords=s3d$xyz.convert(tapply(data$X1,data$X4,mean),
#                               tapply(data$X2,data$X4,mean),
#                               tapply(data$X3,data$X4,mean))
    
#    text(labels=unique(data$X4),s3d.coords,col='red')
    
#}

## find mutual information between IC coefficients

miDist=function(X,Y=NULL,bin=NULL){
    if(is.null(Y)){
        Y=X
    }
    n=dim(X)[1]
    p=dim(X)[2]
    print(n)
    
    if(is.null(bin)){
        bin=floor((log(p)+1+p^0.5)/2)
        print(paste('Bin size: ',as.character(bin),sep=''))
    }
    res=matrix(nrow=n,ncol=n)
    library(entropy)
    
    for(i in 1:n){
        for (j in 1:i){
            xy2d=discretize2d(X[i,],Y[j,],bin,bin)
            #mi=mi.empirical(xy2d)
            h1=entropy(rowSums(xy2d))
            h2=entropy(colSums(xy2d))
            h12=entropy(xy2d)
            mi=(h1+h2-h12)/max(h1,h2)
            res[i,j]=1-mi
            res[j,i]=1-mi
        }
    }
    
    return(res)
}


## the findCor function reports significant correlation between mixing scores and clinical data 
## if the clinical variable is binary, this function build an logistic regression model
## otherwise it takes the y as ordinal and build a linear regression model
## factor levels have to be in desired order before using this function

findCor=function(mix,clustering,clinical.data,p_value=T){
    
    
    clsName=names(table(clustering))
    res=matrix(nrow=dim(mix)[1],ncol=dim(clinical.data)[2])
    colnames(res)=colnames(clinical.data)
    
    for (i in 1:dim(mix)[1]){
        
        for (j in 1:dim(clinical.data)[2]){
            
            model.data=data.frame(y=clinical.data[,j],x=as.numeric(mix[i,]))
            #wprint(model.data)
            
            if (length(levels(clinical.data[,j]))==2){
                
                mod=glm(y~x,data=model.data,family='binomial')
                
            }
            else{
                model.data=data.frame(y=as.numeric(clinical.data[,j]),x=as.numeric(mix[i,]))
                mod=lm(y~x,data=model.data)
            }
            if(p_value){
                res[i,j]=coef(summary(mod))[2,4] ## output can be p-value or coeffcients
            }
            else{
                res[i,j]=coef(summary(mod))[2,1]
            }
        }
    }
    
    return(res)
}

## The fundCor2 function is similar to findCor, the only difference is it takes all y in the clinical data
## as ordinals, i.e., convert it into numerical

findCor2=function(mix,clustering,clinical.data,p_value=T){
    
    
    clsName=names(table(clustering))
    res=matrix(nrow=dim(mix)[1],ncol=dim(clinical.data)[2])
    colnames(res)=colnames(clinical.data)
    
    for (i in 1:dim(mix)[1]){
        
        for (j in 1:dim(clinical.data)[2]){
            
            #model.data=data.frame(y=clinical.data[,j],x=as.numeric(mix[i,]))
            #wprint(model.data)
            
            model.data=data.frame(y=as.numeric(clinical.data[,j]),x=as.numeric(mix[i,]))
            mod=lm(y~x,data=model.data)
            #print(colnames(clinical.data)[j])
            #print(summary(mod))
            if(p_value){
                res[i,j]=coef(summary(mod))[2,4] ## output can be p-value or coeffcients
            }
            else{
                res[i,j]=coef(summary(mod))[2,1]
            } 
        }
    }
    
    return(res)
}