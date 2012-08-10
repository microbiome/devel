# calculate core with given tresholds

# data = your datamatrix, intTr = intensity treshold
# Take Log for the CoV data!!

library(microbiome)
source("CoreFunctions.R")
source("helper_func.R")

# --------------------------------------

# I/O
writedir = "./results/"
load("Atlas.RData")
spe <- atlas[["species"]]$rpa

phylogeny <- read.delim("phylogenyinfo_008.tab")

nb <- 1e2 # was: 1e4 / LL
ncore <- 2

# --------------------------------------

# Define subsets of the data 
data <- list()

samples.healthy <- rownames(atlas.metadata[which(atlas.metadata$health.status == "healthy"),])
data$Healthy <- spe[, samples.healthy]

samples.compromised <- rownames(atlas.metadata[which(atlas.metadata$health.status == "compromised"),])
data$Compromised <- spe[, samples.compromised]

# --------------------------------------

# Estimate min threshold
low.thresh <- estimate.min.threshold(spe)

# Compute cores
Core.list <- lapply(data, createCore, I.thr = low.thresh)

# 3D core plots
mapply(Core3D, Core.list, names(Core.list), writedir)

# 2D core plots
mapply(Core2D, Core.list, names(Core.list), writedir)

# ----------------------------------------------------------

# bootstrap sample of Nhealthy, minimum prevalence 2
dat <- data$Healthy
healthy.microbes <- bootstrap.microbes(dat, Nsample = ncol(dat), minPrev = 2, I.thr = low.thresh, Nboot = nb, ncore = ncore)

# write results, 
write.table(healthy.microbes, file = paste(writedir, "SpeciesCoreFrequencies.txt", sep = ""), sep = "\t", row.names = F, quote = F)

# make a figure
plot_cumulative(healthy.microbes, writedir, "healthy", i.set, type = "cumulative", ylim = c(-1,1), phylogeny = phylogeny)

stop("HERE")

# ----------------------------------------------------------

# rarefaction curve - healthy group only
dat <- data$Healthy
count.microbes <- list()
for (i in 1:ncol(dat)) {
  count.microbes[[i]] <- bootstrap.microbecount(dat, minprev = 1, Nsample = i, I.thr = low.thresh, Nboot = nb)
}

# count.microbes <- lapply(ncol(dat), function (i) {bootstrap.microbecount(dat, minprev = 1, Nsample = i, I.thr = low.thresh, Nboot = nb)})

# -----------------------------------------------------------

pdf(paste(writedir,"rarefaction_curve_healthy.pdf",sep=""))
plot(1:ncol(dat), as.numeric(as.data.frame(count.microbes)[2,]), type = "l", ylim = c(30, 830), xlab = "Sample size", ylab = "# species detected")
lines(as.numeric(as.data.frame(count.microbes)[1,]), lty = "dashed")
lines(as.numeric(as.data.frame(count.microbes)[3,]), lty = "dashed")
dev.off()

# rarefaction curve - size of core
dat <- data$Healthy
core.microbes <- list()
for (i in 1:ncol(dat))
  core.microbes[[i]] <- bootstrap.microbecount(dat, Nsample = i, minprev = i, I.thr = low.thresh, Nboot = nb)

pdf(paste(writedir,"rarefaction_core_healthy.pdf",sep=""))
plot(1:ncol(dat), as.numeric(as.data.frame(core.microbes)[2,]), type = "l", ylim = c(0, 400), xlab = "Sample size", ylab = "# species present in all")
lines(as.numeric(as.data.frame(core.microbes)[1,]), lty = "dashed")
lines(as.numeric(as.data.frame(core.microbes)[3,]), lty = "dashed")
dev.off()

# ---------------------------------------------------------

# Comparisons between groups

pdf(paste(writedir,"Core_comparisons_UC_Healthy.pdf",sep=""))
plot(UC.microbes[,2], healthy.subset.microbes[,2], xlab = "UC core frequency", ylab = "Healthy core frequency", main = "Healthy vs UC Core Frequencies", pch = 20)
lines(c(0,1),c(0,1), col = "Red")
dev.off()

sub.diff <- healthy.subset.microbes
sub.diff[,2] <- healthy.subset.microbes[,2]-UC.microbes[,2]
plot_cumulative(sub.diff,writedir, "Healthy-UC", i.set, type="cumulative",ylim=c(-1,1))
plot_cumulative(sub.diff,writedir, "Healthy-UC all", type="cumulative",ylim=c(-1,1))

pdf(paste(writedir,"rarefaction_curve_healthy_UC.pdf",sep=""))
plot(1:ncol(data$UC),as.numeric(as.data.frame(count.microbes)[2,1:ncol(data$UC)]),type="l",ylim=c(0,640),xlab="Sample size",ylab="# species detected")
lines(as.numeric(as.data.frame(count.microbes)[1,1:ncol(data$UC)]),lty="dashed")
lines(as.numeric(as.data.frame(count.microbes)[3,1:ncol(data$UC)]),lty="dashed")
lines(1:ncol(data$UC),as.data.frame(UC.count.microbes)[2,],col="red")
dev.off()

pdf(paste(writedir,"rarefaction_core_healthy_UC.pdf",sep=""))
plot(1:ncol(data$UC),as.numeric(as.data.frame(core.microbes)[2,1:ncol(data$UC)]),type="l",ylim=c(0,400),xlab="Sample size",ylab="# species present in all")
lines(as.numeric(as.data.frame(core.microbes)[1,1:ncol(data$UC)]),lty="dashed")
lines(as.numeric(as.data.frame(core.microbes)[3,1:ncol(data$UC)]),lty="dashed")
lines(1:ncol(data$UC),as.data.frame(UC.core.microbes)[2,],col="red")
dev.off()

#UC vs. healthy - sample same number of patients from healthy
Nboot <- nb
Nsample <- ncol(data$UC)
boot <- replicate(Nboot,sample(ncol(data$Healthy),Nsample,replace=T),simplify=F)
boot.size <- mclapply(boot,function(x){ 
    Prev <- length(x);
    createCore(data$Healthy[,x], I.thr = low.thresh)[,as.character(Prev)]
}, mc.cores = ncore)

# mix of UC and Healthy
Nsample=ncol(data$UC)
boot=replicate(Nboot,sample(c(sample(grep("UC",colnames(data$All)),Nsample,replace=T),sample(grep("UC",colnames(data$All),invert=T),Nsample,replace=T)),Nsample,replace=T),simplify=F)
boot.uc.he.size=mclapply(boot,function(x){ 
    Prev=length(x); #Prev=round(runif(1,minPrev,length(x)));
    createCore(data$All[,x],I.thr=low.thresh)[,as.character(Prev)]
}, mc.cores = ncore)

#find out maximum vector length
ind.seq.2=names(boot.uc.he.size[[which(unlist(lapply(boot.uc.he.size,length))==max(unlist(lapply(boot.uc.he.size,length))))[1]]])
ind.seq=names(boot.size[[which(unlist(lapply(boot.size,length))==max(unlist(lapply(boot.size,length))))[1]]])

# confidence intervals
boot.size.ci=sapply(ind.seq,function(x){ q=vector("numeric",length(boot.size)); 
    for (i in 1:length(boot.size)) 
       q[i]=boot.size[[i]][x]
       #t1=quantile(q,probs=c(0.05,0.5,0.95),na.rm=T)
       t1=vector("numeric",3)
       t2=sd(as.numeric(q))
       t1[2]=mean(q)
       t1[1]=t1[2]-1.96*t2
       t1[3]=t1[2]+1.96*t2
     return(t1)}) 

boot.uc.he.size.ci=sapply(ind.seq.2,function(x){ q=vector("numeric",length(boot.uc.he.size)); 
    for (i in 1:length(boot.uc.he.size)) 
       q[i]=boot.uc.he.size[[i]][x]
       #t1=quantile(q,probs=c(0.05,0.5,0.95),na.rm=T)
       t1=vector("numeric",3)
       t2=sd(as.numeric(q))
       t1[2]=mean(q)
       t1[1]=t1[2]-1.96*t2
       t1[3]=t1[2]+1.96*t2
     return(t1)}) 

pdf(paste(writedir,"2d_core_confint.pdf",sep=""))
plot(as.numeric(colnames(boot.size.ci)),boot.size.ci[2,],type="l",xlab="log intensity",ylab="number of species",main="Common core: Healthy")
lines(as.numeric(colnames(boot.size.ci)),boot.size.ci[1,],type="l",lty=2)
lines(as.numeric(colnames(boot.size.ci)),boot.size.ci[3,],type="l",lty=2)
dev.off()

pdf(paste(writedir,"2d_core_confint_UC.pdf",sep=""))
plot(as.numeric(colnames(boot.size.ci)),boot.size.ci[2,],type="l",xlab="log intensity",ylab="number of species",main="Common core: Healthy vs. UC (red)")
lines(as.numeric(colnames(boot.size.ci)), boot.size.ci[1,], lty = 2)
lines(as.numeric(colnames(boot.size.ci)), boot.size.ci[3,], lty = 2)
lines(as.numeric(rownames(Core.list$UC)), Core.list$UC[,12], col = "red")
dev.off()

pdf(paste(writedir,"2d_core_confint_healthyUC.pdf",sep=""))
d1=1:which(boot.size.ci[3,]==0)[1]
plot(as.numeric(colnames(boot.size.ci))[d1],boot.size.ci[2,d1],type="l",xlab="log intensity",ylab="number of species",main="Common core: Healthy vs. UC")
lines(as.numeric(colnames(boot.size.ci))[d1],boot.size.ci[1,d1],lty=2)
lines(as.numeric(colnames(boot.size.ci))[d1],boot.size.ci[3,d1],lty=2)
lines(as.numeric(rownames(Core.list$UC))[d1],Core.list$UC[d1,12],col="red")
lines(as.numeric(colnames(boot.uc.he.size.ci))[d1],boot.uc.he.size.ci[2,d1],col="blue")
legend(3.2,400,c("Core from Healthy","Core from Healthy+UC","Core from UC"), fill=c("black","blue","red")) 
dev.off()

#### VENN
require(gplots)
source("../Rscripts/phylogeneticEnrichments_prof008.R")
PH.i <- unique(phylogeny[,1:3])

UC95=round(1.96*(UC.core.microbes[[12]][3]-UC.core.microbes[[12]][2])+UC.core.microbes[[12]][2])
HE95=round(1.96*(core.microbes[[12]][3]-core.microbes[[12]][2])+core.microbes[[12]][2])
coresets=list()
coresets[["UC"]]=as.character(UC.microbes$Microbe[order(UC.microbes$Frequency,decreasing=T)[1:UC95]])
coresets[["Healthy"]]=as.character(healthy.microbes$Microbe[order(healthy.microbes$Frequency,decreasing=T)[1:HE95]])
pdf(paste(writedir,"Venn_Healthy_vs_UC95pct.pdf",sep=""))
venn(coresets)
dev.off()
write.table(file=paste(writedir,"UC_core95.txt",sep=""),coresets[["UC"]],quote=F,row.names=F,col.names=F)
write.table(file=paste(writedir,"Healthy_core95.txt",sep=""),coresets[["Healthy"]],quote=F,row.names=F,col.names=F)

shared=coresets[["UC"]][coresets[["UC"]] %in% coresets[["Healthy"]]]
UCown=coresets[["UC"]][!(coresets[["UC"]] %in% coresets[["Healthy"]])]
HEown=coresets[["Healthy"]][!(coresets[["Healthy"]] %in% coresets[["UC"]])]
xx.UC=phylogeneticEnrichments(UCown,PH.i,origlevel="species",maplevel="level.2")[[1]]
xx.HE=phylogeneticEnrichments(HEown,PH.i,origlevel="species",maplevel="level.2")[[1]]
xx.SH=phylogeneticEnrichments(shared,PH.i,origlevel="species",maplevel="level.2")[[1]]
write.table(file=paste(writedir,"UC95_only_l2.txt",sep=""),xx.UC,quote=F)
write.table(file=paste(writedir,"HE95_only_l2.txt",sep=""),xx.HE,quote=F)
write.table(file=paste(writedir,"UC95_HE95_l2.txt",sep=""),xx.SH,quote=F)

# with expected core size
UC50 <- round(UC.core.microbes[[12]][2])
HE50 <- round(core.microbes[[12]][2])
coresets <- list()
coresets[["UC"]] <- as.character(UC.microbes$Microbe[order(UC.microbes$Frequency,decreasing=T)[1:UC50]])
coresets[["Healthy"]] <- as.character(healthy.microbes$Microbe[order(healthy.microbes$Frequency,decreasing=T)[1:HE50]])
pdf(paste(writedir,"Venn_Healthy_vs_UC50pct.pdf",sep=""))
venn(coresets)
dev.off()

write.table(file=paste(writedir,"UC_core50.txt",sep=""),coresets[["UC"]],quote=F,row.names=F,col.names=F)
write.table(file=paste(writedir,"Healthy_core50.txt",sep=""),coresets[["Healthy"]],quote=F,row.names=F,col.names=F)

shared <- coresets[["UC"]][coresets[["UC"]] %in% coresets[["Healthy"]]]
UCown <- coresets[["UC"]][!(coresets[["UC"]] %in% coresets[["Healthy"]])]
HEown <- coresets[["Healthy"]][!(coresets[["Healthy"]] %in% coresets[["UC"]])]
xx.UC <- phylogeneticEnrichments(UCown,PH.i, origlevel="species", maplevel = "level.2")[[1]]
xx.HE <- phylogeneticEnrichments(HEown,PH.i, origlevel="species", maplevel = "level.2")[[1]]
xx.SH <- phylogeneticEnrichments(shared,PH.i,origlevel="species", maplevel = "level.2")[[1]]
write.table(file = paste(writedir,"UC_only_l2.txt",sep=""), xx.UC, quote = F)
write.table(file = paste(writedir,"HE_only_l2.txt",sep=""), xx.HE, quote = F)
write.table(file = paste(writedir,"UC_HE_l2.txt",sep=""),   xx.SH, quote = F)







