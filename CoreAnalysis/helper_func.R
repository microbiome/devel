doPCA_RDA=function(y,d.i,dtype,prname,study,g1=NULL,wLines=F,wLegend=T,wText=NULL,maxRDA=20,xlim=NULL,ylim=NULL,nperm=5000){
  rdatestiHL=NULL
  colnames(d.i)=c("sample.ID","Group.info")
  ########### PCA ############
  mm=prcomp(t(y),retx=T)
  #screeplot
  pdf(paste(study,"_",prname,"_scree.pdf",sep=""))
  plot(mm,main=paste(study,prname,": Eigenvalues"))
  dev.off()

  pdf(paste(study,prname,dtype,"PCA_12.pdf",sep="_"))
  if (is.logical(wLegend) & (wLegend)){
    if (min(mm$x[,1])<0)
      xm=min(mm$x[,1])*1.7
    else
      xm=min(mm$x[,1])*0.7
  }else
    xm=min(mm$x[,1])

  plot(mm$x[,1],mm$x[,2],pch=20,xlab="PC 1",ylab="PC 2",main=paste(study,prname,": PCA"),xlim=c(xm,max(mm$x[,1])))
  RGBcol<-colorRampPalette(c("Red", "Green","Blue"),space="Lab")(length(levels(d.i$Group.info)))
  if (length(levels(d.i$Group.info))>1){
   for (i in 1:length(levels(d.i$Group.info))){
     t1=which(d.i$Group.info==levels(d.i$Group.info)[i])
     points(mm$x[t1,1],mm$x[t1,2],pch=20,col=RGBcol[i])
     if (wLines)
       lines(mm$x[t1,1],mm$x[t1,2],pch=20,col=RGBcol[i])
   }
   if (is.logical(wLegend)){
     if (wLegend)
       legend(x=xm,y=max(mm$x[,2]),pch=20,legend=levels(d.i$Group.info),col=RGBcol,ncol=2)
   }else{
      legend(x=xm,y=max(mm$x[,2]),pch=20,legend=wLegend,col=RGBcol,ncol=2)
   }
   if (!is.null(wText))
     text(mm$x[,1],mm$x[,2], labels=as.character(d.i[,wText]), cex=0.8)
  }
  dev.off()
  ########### RDA ############
  if (!is.null(g1)){
    require(vegan)
    for (i in 1:length(g1)){
      stmp=d.i$Group.info %in% g1[[i]]
      rdatesti <- rda(t(y) ~ Group.info,data=d.i,subset=stmp)
      ## piirtää kaikki lajit kuvaan; melko messy
      #plot(rdatesti, choices=c(1,2), display=c("sp","bp"), type="text", scaling=3,xlim=c(-1.5,2), ylim=c(-2,1.2)) 
      # arrows?
      ## piirtää maxRDA eniten kontribuoivaa lajia 
      g <- d.i$Group.info[stmp]
      if (!is.null(wText))
        l=as.character(d.i[,wText])
      else
        l <- as.character(d.i$Group.info[stmp])

      temp <- scores(rdatesti)
      loadings <- rowSums(temp$species^2)
      rdatestiHL <- rda(t(y[names(sort(loadings, decreasing=T))[1:maxRDA],stmp]) ~ g)
      print(permutest(rdatestiHL,permutations=nperm))
      # x11(width=10, height=10)
      #dev.new()
      pdf(paste(study,"_",prname,"_RDA_",dtype,"_",paste(g1[[i]],collapse="_"),".pdf",sep=""))
      plot(rdatestiHL, choices=c(1,2), type="text", scaling=3, display=c("sp","bp"),main=paste(study,prname,"RDA:",paste(g1[[i]],collapse=" vs ")),xlim=xlim,ylim=ylim)
      text(rdatestiHL, choices=c(1,2), scaling=3, labels=as.character(l), cex=0.8)
     ##yllä labels=ryhmäkoodi; alla labeleina näyteIDt 
     # text(rdatestiHL, choices=c(1,2), scaling=3, cex=0.7)
      dev.off()
    }
  }# RDA
  return(list(d.i=d.i,y=y,rda=rdatestiHL))
}

lambda.cv=function(x,y,l2,s,K=10){
  r1=list()
  for (i in 1:length(l2)){
     r1[[i]]=cv.enet(x, y,K=10,lambda=l2[i], s=s, mode="fraction",trace = F, plot.it = F, se = F,normalize=F)
  }
  minl2=which(unlist(lapply(r1,function(x) min(x$cv)))==min(unlist(lapply(r1,function(x) min(x$cv)))))[1]
  mins=which(r1[[minl2]]$cv==min(r1[[minl2]]$cv))
  return(data.frame(l2=l2[minl2],s=s[mins]))
}

thr.cv=function(xx,yy,ss,lambda,s,K=10){
  cvs=cv.folds(length(yy),folds=K)
  residmat=matrix(0,K,length(ss))
  for (i in seq(K)) {
     omit <- cvs[[i]]
     fit <- enet(xx[-omit, ], yy[-omit], lambda = lambda, normalize=F)
     fit <- predict(fit, xx[omit, , drop = FALSE], s = s, mode = "fraction")$fit
     residmat[i,] <- sapply(ss,function(x) sum(abs((fit>x)-y[omit]))/length(omit))
  }
  t1=colMeans(residmat)
  CC=ss[which(t1==min(t1))[1]]
  return(data.frame(c=CC,cverr=min(t1)))
}

do.enet.class=function(x,y,l2,s,SplitRatio=1/2,K=10,nmod=10){
  require(elasticnet)
  # cross validation: 
  require(caTools)
  test=which(sample.split(y, SplitRatio = 1/2))
  train=setdiff(1:length(y),test)
  tres=list()
  minC=1
  for (j in 1:nmod){
    rr=lambda.cv(x[train,],y[train],l2,s)
    while (rr$s==0)
      rr=lambda.cv(x[train,],y[train],l2,s)
    CC=thr.cv(x[train,],y[train],s,rr$l2,rr$s,K)

    trained_enet <- enet(x[train,],y[train],lambda=rr$l2,normalize=F)
    pred.train=predict.enet(trained_enet, x[train,], s=rr$s,type="fit",mode="fraction")
    coef=predict.enet(trained_enet, s=rr$s,type="coef",mode="fraction")
    # test data
    pred.test=predict.enet(trained_enet, x[test,], s=rr$s,type="fit",mode="fraction")
    fulld.enet <- enet(x,y,lambda=rr$l2,normalize=F)
    coef.full=predict.enet(fulld.enet, s=rr$s,type="coef",mode="fraction")
    pred.full=predict.enet(fulld.enet, x, s=rr$s,type="fit",mode="fraction")

    coef=coef$coefficients[which(coef$coefficients!=0)] 
    trainC=sum(abs((pred.train$fit>CC$c)-y[train]))/length(train)
    testC=sum(abs((pred.test$fit>CC$c)-y[test]))/length(test)
    dumb=(sum(y[train]==max(y[train]))/length(train))
    trainC.full=sum(abs((pred.full$fit>CC$c)-y))/length(y)
    coef.full=coef.full$coefficients[which(coef.full$coefficients!=0)]
    if (testC<minC){
      minC=testC
      tres=list(coef=coef,trainC=trainC,testC=testC,dumb=dumb,coef.full=coef.full,trainC.full=trainC.full,CC=CC,rr=rr,enet=trained_enet)
    }
  }
  return(tres)
}

do.enet.regress=function(x,y,l2,s,SplitRatio=1/2,K=10,nmod=10){
  require(elasticnet)
  # cross validation: 
  require(caTools)
  test=which(sample.split(y, SplitRatio = 1/2))
  train=setdiff(1:length(y),test)
  tres=list()
  minC=1
  for (j in 1:nmod){
    rr=lambda.cv(x[train,],y[train],l2,s)
    while (rr$s==0)
      rr=lambda.cv(x[train,],y[train],l2,s)
    CC=thr.cv(x[train,],y[train],s,rr$l2,rr$s,K)

    trained_enet <- enet(x[train,],y[train],lambda=rr$l2,normalize=F)
    pred.train=predict.enet(trained_enet, x[train,], s=rr$s,type="fit",mode="fraction")
    coef=predict.enet(trained_enet, s=rr$s,type="coef",mode="fraction")
    # test data
    pred.test=predict.enet(trained_enet, x[test,], s=rr$s,type="fit",mode="fraction")
    fulld.enet <- enet(x,y,lambda=rr$l2,normalize=F)
    coef.full=predict.enet(fulld.enet, s=rr$s,type="coef",mode="fraction")
    pred.full=predict.enet(fulld.enet, x, s=rr$s,type="fit",mode="fraction")

    coef=coef$coefficients[which(coef$coefficients!=0)] 
    trainC=sum(abs((pred.train$fit>CC$c)-y[train]))/length(train)
    testC=sum(abs((pred.test$fit>CC$c)-y[test]))/length(test)
    dumb=(sum(y[train]==max(y[train]))/length(train))
    trainC.full=sum(abs((pred.full$fit>CC$c)-y))/length(y)
    coef.full=coef.full$coefficients[which(coef.full$coefficients!=0)]
    if (testC<minC){
      minC=testC
      tres=list(coef=coef,trainC=trainC,testC=testC,dumb=dumb,coef.full=coef.full,trainC.full=trainC.full,CC=CC,rr=rr,enet=trained_enet)
    }
  }
  return(tres)
}

SLDA.CV=function(x,y,l2=c(0.05,0.1),SplitRatio=1/3,K=4){
   source("SLDAPathway.R")
   t1=sample.split(yy,SplitRatio=1/3)
   xxtr=as.matrix(x[t1==F,])
   yytr=y[t1==F]
   xxts=as.matrix(x[t1==T,])
   yyts=y[t1==T]
   n = length(yytr)
   n1 = sum(yytr)
   n2 = n-n1 
   f1 = split(sample(seq(n1)), rep(1:K, length = n1))
   f2 = split(sample(seq(n2)), rep(1:K, length = n2))
   minpred=0
   bestl=1
   for (j in 1:length(l2)){
     bloatpred=0
     for (i in seq(K)) {
       omit = c(f1[[i]], f2[[i]])
       ll=SLDA.PathwayTest(xxtr[-omit,],yytr[-omit],lambda2=l2[j])
       Xn=t(as.vector(ll$loadings))%*% t(xxtr[omit,])
       bloatpred=bloatpred+abs(colttests(t(Xn), as.factor(yytr[omit]),tstatOnly=T)$statistic)
     }
     if (bloatpred>minpred){
        bestl=j
        minpred=bloatpred
     }
   }
   ll=SLDA.PathwayTest(xxtr,yytr,lambda2=l2[bestl])
   Xn=t(as.vector(ll$loadings))%*% t(xxts)
   bloatpred=colttests(t(Xn), as.factor(yyts))
   list(ll=ll,pred=bloatpred,bestl=l2[bestl])
}

learn.glmnet=function(xx,yy, SplitRatio = 3/4){
  require(caTools) 
  require(glmnet)
  itrain=sample.split(yy, SplitRatio = SplitRatio, group = NULL )
  itest=(itrain==F)
  fitobj=list()
  lp=list()
  for (i in 1:10){
    fitobj[[i]]=cv.glmnet(t(xx[,itrain]),yy[itrain],family="binomial")
    lpr=predict(fitobj[[i]],newx=t(xx)[which(itest==T),])
    lp[[i]]=sum(1-abs((as.numeric(lpr>0))-yy[itest]))/sum(itest)
  }
  lp=unlist(lp)
  bi=which(lp==max(lp))
  minc=nrow(xx)
  mini=0
  for (i in bi){
    lr=coef(fitobj[[i]])
    fitcoef=lr[which(lr!=0),]
    fitcoef=fitcoef[2:length(fitcoef)]
    if (length(fitcoef)<minc){
      minc=length(fitcoef)
      mincoef=fitcoef
      mini=i
    }
  }
  return(list(fitobj=fitobj[[mini]],coef=mincoef,train=itest,test=itrain,pred=lp[[mini]]))
}

