core.sum <- function(data,intTr,prevalenceTr){
  d.bin <- data>=intTr
  prevalences <- rowSums(d.bin)
  nOTUs <- sum(prevalences>=prevalenceTr)# jos haluat titetää lajit, älä summaa!
  return(nOTUs)
}

core.which <- function(data,intTr,prevalenceTr){
  d.bin <- data>=intTr
  prevalences <- rowSums(d.bin)
  nOTUs <- as.numeric(prevalences>=prevalenceTr)
  return(nOTUs)
}

#creates a coreMatrix containing the info from the two vectors created
createCore=function(data,I.thr=1.8){
   ##Prevalence vector
   p.seq <- 1:ncol(data)
   #intensity vector
   #for TR set the min to static 1.8 - gives a better graph
   i.seq <- seq(I.thr, max(data), 0.08)
   #i.seq <- seq(min(data), max(data), 0.08)
   coreMat <- matrix(,nrow=length(i.seq), ncol= length(p.seq), dimnames=list(i.seq,p.seq))
   for(i in i.seq){for(p in p.seq){coreMat[as.character(i),as.character(p)] <- core.sum(data,i,p)}}
   return(coreMat)
}

#3D figure:
Core3D=function(coreMat,fname,writedir){
   ##Renaming the components for the figure
  coreSize <- coreMat
  MinimumPrevalence <- as.numeric(colnames(coreMat))
  MinimumLogIntensity <- as.numeric(rownames(coreMat))

  ##creating the figure
  #persp(MinimumLogIntensity, MinimumPrevalence, coreSize, theta=60, phi=5)
 
  #persp(MinimumLogIntensity, MinimumPrevalence, coreSize, theta=20, phi=5, main="Abundance core - Healthy", col="light blue")
  #persp(MinimumLogIntensity, MinimumPrevalence, coreSize, theta=60, phi=-5, main="Abundance core - Healthy", col="light blue")
  #persp(MinimumLogIntensity, MinimumPrevalence, coreSize, theta=40, phi=5, main="Abundance core - Healthy", col="light blue", axes=T, ticktype="detailed", nticks=9)
   pdf(paste(writedir,"Core3D_",fname,".pdf",sep=""))
   persp(MinimumLogIntensity, MinimumPrevalence, coreSize, theta=60, phi=5, main=paste("Abundance core:",fname), col="light blue", axes=T, ticktype="detailed", nticks=9, shade = 0.58,cex.axis=0.5)
   dev.off()
}

#2D figure
Core2D=function(coreMat,fname,writedir,colnum=NULL){
  i.seq=as.numeric(rownames(coreMat))
  if (is.null(colnum))
    colnum=ncol(coreMat)
  pdf(paste(writedir,"Core2D_",fname,".pdf",sep=""))
  plot(i.seq, coreMat[,colnum], main=paste("Common core: ",fname," (",colnum," subjects)",sep=""),  xlab="log intensity", ylab="number of Species", type="o", pch=16, cex.axis=1.5, cex.main=1.5, cex=1.5)
  dev.off()
}

bootstrap.microbes=function(D,Nsample,minPrev,Nboot=10000,I.thr=1.8,ncore=22){
   require(multicore)
   boot=replicate(Nboot,sample(ncol(D),Nsample,replace=T),simplify=F)

   # below:choose intensity such that there is at least one bacteria fulfilling prevalence criterion
   boot.which=mclapply(boot,function(x){ 
      Prev=round(runif(1,minPrev,length(x)));
      Imax=max(apply(D[,x],1,function(xx) quantile(xx,probs=(1-Prev/length(x)))))
      Insty=runif(1,I.thr,Imax)
      core.which(D[,x],Insty,Prev)
   },mc.cores=ncore)
   boot.prob=rowSums(as.data.frame(boot.which))/Nboot
   return(data.frame(Microbe=rownames(D),Frequency=boot.prob))
}

bootstrap.microbecount=function(D,Nsample,minprev=1,Nboot=10000,I.thr=1.8,ncore=22){
   require(multicore)
   boot=replicate(Nboot,sample(ncol(D),Nsample,replace=T),simplify=F)
   # below:choose intensity such that there is at least one bacteria fulfilling prevalence criterion
   if (Nsample>1)
     boot.which=mclapply(boot,function(x){ 
        Imax=max(apply(D[,x],1,min))
        Insty=runif(1,I.thr,Imax)
        sum(rowSums(D[,x]>Insty) >= minprev)
     },mc.cores=ncore)
   else{
     boot.which=lapply(boot,function(x){ 
        Imax=max(D[,x])
        Insty=runif(1,I.thr,Imax)
        return(sum(D[,x]>=Insty))
     })
   }
   boot.prob=as.matrix(as.data.frame(boot.which,check.names=F))
   t1=quantile(boot.prob,probs=c(0.05,0.5,0.95))
   #t1=vector("numeric",3)
   #t2=sd(as.numeric(boot.prob))
   t1[2]=mean(boot.prob)
   #t1[1]=t1[2]-t2
   #t1[3]=t1[2]+t2
   #t1[1]=t1[2]-1.96*t2
   #t1[3]=t1[2]+1.96*t2
   return(t1)
}


plot_cumulative <- function(d.sub, writedir, fname, i.set = NULL, type = "cumulative", ylim = NULL, phylogeny = NULL){

#d.sub <- healthy.microbes; fname = "healthy"; type = "cumulative"; ylim = c(-1,1)

   #PH.i <- read.delim("./blanket analysis/phylogenyinfo_008.tab")
   
   PH.i=unique(phylogeny[,1:3])
   PH.i=PH.i[order(PH.i[,3]),]
   d.sub$Microbe=PH.i[,1]
   d.sub=d.sub[order(d.sub[,2],decreasing=T),]
   pdf(paste(writedir,"cumsum_",fname,".pdf",sep=""))
   if (is.null(i.set))
      i.set=1:length(levels(d.sub[,1]))
   colmap=colorRampPalette(c("Red", "Green","Blue"),space="rgb")(length(levels(d.sub[,1])))
   cnt=1;
   i.accept=vector("logical",length(levels(d.sub[,1])))
   for (i in i.set){
     l.res=as.numeric(d.sub[,1]==levels(d.sub[,1])[i])
     if (type=="cumulative"){
#        l.res=l.res/sum(l.res)
#        l.res[l.res==0]=-1/sum(l.res==0)
        out=cumsum(l.res)
        #null=replicate(10000,cumsum(sample(l.res,length(l.res))))
        #null.med=apply(null,1,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
     }
     if (type=="gsea"){
#        l.res=l.res/sum(l.res)
#        l.res[l.res==0]=-1/sum(l.res==0)
        l.res[l.res==0]=-1
        out=l.res
        for (j in 2:length(l.res))
           out[j]=max(l.res[j]+out[j-1],0)
     }
     t1=seq(max(d.sub[,2]),min(d.sub[,2]),-0.01)
     out=vector("numeric",length(t1))
     null.cum=matrix(NA,length(t1),3)
     for (j in 1:length(t1)){
       out[j]=sum(l.res*(d.sub[,2]>t1[j]))
       null.cum[j,]=quantile(replicate(1000,sum(sample(l.res,length(l.res))*(d.sub[,2]>=t1[j]))),probs=c(0.025,0.5,0.975))
     }
     yplot=(out-null.cum[,2])/max(abs(out-null.cum[,2]))
     if (sum(out<null.cum[,1])>0 | sum(out>null.cum[,3])>0){
      if (cnt==1){
        if (is.null(ylim))
           plot(t1,yplot,type="l",main=paste(type,"prevalence of L1 taxa (",fname,")"),xlim=c(max(t1),min(t1)),col=colmap[i],ylab="proportion of total",xlab="Frequency")
        else
           plot(t1,yplot,type="l",main=paste(type,"prevalence of L1 taxa (",fname,")"),ylim=ylim,xlim=c(max(t1),min(t1)),col=colmap[i],ylab="proportion of total",xlab="Frequency")
      }else
         lines(t1,yplot,col=colmap[i])
      cnt=cnt+1;
      i.accept[i]=T
    }
   }
   legend(max(t1),1,levels(d.sub[,1])[which(i.accept==T)],fill=colmap[which(i.accept==T)],cex=0.5)
   dev.off()
}