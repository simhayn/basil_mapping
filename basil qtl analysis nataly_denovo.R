#*****************************
#Analyze QTL data of basil
#(BA17 Bc5F3 44-4-3 x Deep Purple)F2 population 
#SNP's from Lara
#GBS sequenced by Genewiz
#Nataly Simhayev
#*****************************
####       R=0       S=1      for FOB BDM 
####       R=1       S=9      for Cold

#### Setting up ####
#clean up R workspace
rm(list = ls())
#setwd("~/mydata/denovo basil")
#load("~/mydata/denovo basil/.RData")
setwd("C:/Users/Public/R nataly/denovo")
load("C:/Users/Public/R nataly/denovo/denovo.Rdata")


#load packages
library(qtl)
library(ASMap)
library(dplyr)

#import .csv file named "basil.csv", for bc error type genotypes=c("A","H","B")
basil<-read.cross("csv",".","basil.csv",map.function="kosambi",
                  genotypes = c("AA","AB","BB"))
phenames(basil)
cbind(1:14,phenames(basil))
#normal phenos
n<-11

#jitter the markers
basil<-jittermap(basil)

thcol<-c('blue','green','red')
##colorlist<-brewer.pal(8,"Set1")

#### Phenotype ####
#phenotypic distribution
parents<-read.csv("parents' pheno.csv")
row.names(parents)<-parents[,1];parents<-parents[,-1]
As<-rbind(matrix(0,12,3),matrix(5,12,3))
Fs<-rbind(matrix(0,12,3),matrix(5,12,3))
parents<-cbind.data.frame(parents[,-4],AG=As[,1],AYS=As[,2],AYL=As[,3],
                          FOB1=Fs[,1],FOB2=Fs[,2],FOB3=Fs[,3],Cold=parents[,4])
parents<-cbind.data.frame(parents[,1:9],F_AUDPC=rbind(matrix(0,12),matrix(c(90,142),12)),parents[,10])

#histograms or barplots
{tiff("phenoDistrib%d.tiff",700,650)
par(mfrow=c(2,3),cex=1.3)
#par(mfrow=c(2,3),cex=0.75)
l<-1
for(i in c(1:10,12,15,13,11)){
if(i<12){plotPheno(basil,i,ylab="Frequency",xlab=if(i<7){"Purple intensity"}else{"Disease intensity"})}
if(i==15){par(mfrow=c(2,3),cex=1.3)}
if(i>11&i!=15){
BP=plotPheno(basil,i,ylab="Frequency",xlab="Disease resistance",names.arg=c("R","S"))
text(BP, table(basil$pheno[, i]), labels=table(basil$pheno[, i]), pos=1)
points(0.7,1.6,pch=25,cex=1.3,bg="palegreen1")
points(1.9,1.6,pch=25,cex=1.3,bg="violet")}
if(i<10){points(if(i<4){0.1}else{0.7},1.6,pch=25,cex=1.3,bg="palegreen1")}
if(i<4){points((mean(parents[13:24,i])-0.4),1.6,pch=25,cex=1.3,bg="violet")}
if(i==10){points((mean(parents[1:12,i])-0.1),1.6,pch=25,cex=1.3,bg="palegreen1")
points(mean(parents[13:24,i]),1.6,pch=25,cex=1.3,bg="violet")}
if(i==11){points(mean(parents[1:12,i],na.rm=T),1.6,pch=25,cex=1.3,bg="palegreen1")
points(mean(parents[13:24,i],na.rm=T),1.6,pch=25,cex=1.3,bg="violet")}
if(i%in%c(5,6)){points(4.3,1.6,pch=25,cex=1.3,bg="violet")}
if(i%in%c(4,9)){points(6.7,1.6,pch=25,cex=1.3,bg="violet")}
if(i%in%c(7,8)){points(5.5,1.6,pch=25,cex=1.3,bg="violet")}
if(l>6)l<-1
if(i%in%c(1:10,12)){if(i==7){l<-1} 
  mtext(LETTERS[l],adj=-0.1,cex=1.4)}
l<-l+1
}
dev.off()
}



#### Correlation ####
library(GGally)
#fusarium
mat<-basil$pheno[c(7:10,12)]
names(mat)<-phenames(basil)[c(7:10,12)]
tiff("correlation fus.tiff",res=100)
ggcorr(mat, method=c("pairwise.complete.obs", "pearson"),label=T)
dev.off()
#anthocyanin
mat<-basil$pheno[1:6]
names(mat)<-phenames(basil)[1:6]
tiff("correlation ant.tiff",res=100)
ggcorr(mat, method=c("pairwise.complete.obs", "pearson"),label=T)
dev.off()
#all pheno but cold corr
mat<-basil$pheno[c(1:10,12:13)]
names(mat)<-phenames(basil)[c(1:10,12:13)]
tiff("correlation allpheno.tiff",res=100)
ggcorr(mat, method=c("pairwise.complete.obs", "pearson"),label=T)
dev.off()

####Scanone####
#calculate probabilities (necessary for scanone)
#Haley-knott regression
basil<-calc.genoprob(basil,2,map.function="kosambi")
#normal
#genome scan for Single-QTL
##scan1<-scanone(basil,pheno.col=c(1:n),method="hk")
load(file="sc1.hk.Rdata")
#permutation test
load(file="sc1perm.hk.Rdata")

####Not run!
operm <- vector("list", 100)
for(i in 1:100){
  operm[[i]]<-scanone(basil,pheno.col=c(1:n),method="hk",
                      n.perm=9,n.cluster=2)
}
scan1perm<-do.call("rbind", operm)
save(scan1perm,file="sc1perm.hk.Rdata")
####

#setting a threshold
(thresh1.hk<-summary(scan1perm,alpha=c(0.63,0.1,0.05)))

#EM algorithm maximum likelihood
#binary
basil<-calc.genoprob(basil,3,map.function="kosambi")
#genome scan for Single-QTL
##scan1.bin<-scanone(basil,pheno.col=c(12:13),method="em",model="binary")
load(file="sc1.bin.em3.Rdata")
#permutation test
load(file="sc1perm.bin.em3.Rdata")

####Not run!
scan1perm.bin1<-scanone(basil,pheno.col=c(12:13),method="em",
                        model="binary",n.perm=1000,n.cluster=2)
operm <- vector("list", 100)
for(i in 1:100){
  operm[[i]]<-scanone(basil,pheno.col=c(12:13),method="em",
            model="binary",n.perm=10,n.cluster=2)
  }
scan1perm.bin1<-do.call("rbind", operm)
save(scan1perm.bin,file="sc1perm.bin.em3.Rdata")
####

#setting a threshold
(thresh2.em<-summary(scan1perm.bin,alpha=c(0.63,0.1,0.05)))

####addqtl####
#normal
(qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",
                 alpha=0.63,ci.function="bayesint"))
basil<-calc.genoprob(basil,2,map.function="kosambi")
out.aq<-list();rqtl<-list()
for(i in 1:n){
  if(length(qtlist[[i]][,1])>0){
    p<-phenames(basil)[i]
    qtlobj<-makeqtl(basil,qtlist[[i]][,1],qtlist[[i]][,2],what="prob")
    rqtl[[p]]<-refineqtl(basil,p,qtlobj,method="hk")
    out.aq[[p]]<-addqtl(basil,qtlist[[i]][,1],p,rqtl[[p]],
                        method="hk",verbose=T)
  }
}
#binary
(qtlist.bin<-summary(scan1.bin,perms=scan1perm.bin,
                     format="tabByCol",alpha=0.63,ci.function="bayesint"))
basil<-calc.genoprob(basil,3,map.function="kosambi")
out.aq.bin<-list();rqtl.bin<-list()
for(i in 1:2){
  if(length(qtlist.bin[[i]][,1])>0){
    p<-phenames(basil)[i+n]
    qtlobj.bin<-makeqtl(basil,qtlist.bin[[i]][,1],qtlist.bin[[i]][,2],
                        what="prob")
    rqtl.bin[[p]]<-refineqtl(basil,p,qtlobj.bin,method="hk",model="binary")
    out.aq.bin[[p]]<-addqtl(basil,qtlist.bin[[i]][,1],p,rqtl.bin[[p]],
                            maxit=1e+9,tol=0.05,method="hk",model="binary",verbose=T)
  }
}

####plot addqtl####
plotAddqtl<-function(x,bin=F,list=qtlist,aq=out.aq,thresh=thresh1.hk,
                     mfrow=c(2,2),second=F,LETTERs=T,l=1){
  if(!second)par(mfrow=mfrow)
  par(cex.lab=1.5,cex.axis=1.3,cex.main=1.7,cex.sub=1.3)
  for(i in x){
    if(length(list[[i]][,1])>0){
      p<-phenames(basil)[i];if(bin)p<-phenames(basil)[i+n]
      plot(aq[[p]],alternate.chrid=nrow(list[[p]])>3,ylab="LOD")
      abline(h=thresh[,i],lty='dotted',lwd=2,col=thcol)
      for(j in 1:3){
        if(thresh[j,i]/par('usr')[4]<1){
          mtext(rownames(thresh)[j],cex=0.9,font=2,
                adj=thresh[j,i]/(par('usr')[4]-0.1),
                side=4,col=thcol[j])}
      }
      title(p)
      if(nrow(list[[p]])==1)title(sub=list[[p]][,1])
      if (LETTERs) mtext(LETTERS[l],adj=0,cex=1.3);l<-l+1;if(l>6)l<-1
    }
  }
}

tiff("addqtl.ant.tiff")
plotAddqtl(1:6,mfrow=c(2,3))
dev.off()
tiff("addqtl.fus.tiff")
plotAddqtl(7:10)
plotAddqtl(1,T,qtlist.bin,out.aq.bin,thresh2.em,second=T,l=4)
dev.off()
tiff("addqtl.FOB3_bin.tiff",850,550)
plotAddqtl(1,T,qtlist.bin,out.aq.bin,thresh2.em,l=6)
dev.off()
tiff("addqtl.BDM.tiff")
plotAddqtl(2,T,qtlist.bin,out.aq.bin,thresh2.em,LETTERs=F)
dev.off()

####plotScanone####
#plot qtl peaks of scanone
plot.sc<-function (x,bin=F,sc=scan1,thresh=thresh1.hk,LETTERs=T,l=1,
                   mf=c(2,2),second=F,first=F){
  if (!second) par(mfrow=mf)
  par(cex.lab=1.5,cex.main=1.5)
  for (i in x) {
    p<-phenames(basil)[i];if(bin)p<-phenames(basil)[i+n]
    plot(sc,lodcolumn=i,main=p,ylab="LOD",bandcol="gray80",
         ylim=c(0,max(sc[,2+i])+0.5),alternate.chrid = T)
    abline(h=thresh[,i],lty='dotted',lwd=2,col=thcol)
    for(j in 1:3){
      if(thresh[j,i]/par('usr')[4]<1){
        mtext(rownames(thresh)[j],side=4,font=2,
              adj=thresh[j,i]/(par('usr')[4]-0.2),col=thcol[j])
      }
    }
    if (LETTERs) mtext(LETTERS[l],adj=0,cex=1.3);l<-l+1;if(l>6)l<-1
  }
  if (first==F) par(mfrow=c(1,1))
}
tiff("plot.sc1.hk%d.tiff",1100,700,res=96)
plot.sc(1:6)
plot.sc(7:10,first=T)
plot.sc(1,T,sc=scan1.bin,thresh=thresh2.em,second=T,l=5)
plot.sc(2,T,sc=scan1.bin,thresh=thresh2.em,LETTERs=F,first=T)
plot.sc(11,LETTERs=F,second=T)
dev.off()

####QTLs list and their LGs####
#normal
{intpPhen<-vector();effpPhen<-vector()
for (i in 1:n) {
  (sc1<-summary(scan1,perms=scan1perm,
                    alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
  p<-phenames(basil)[i]
  cat(nrow(sc1),"QTLs in",p,"\n")
  if(nrow(sc1)>1) intpPhen<-c(intpPhen,p)
  if(nrow(sc1)==1) effpPhen<-c(effpPhen,p)}
#binary
intpPhen.bin<-vector();effpPhen.bin<-vector()
for (i in 1:2) {
  (sc1<-summary(scan1.bin,perms=scan1perm.bin,
                    alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
  p<-phenames(basil)[i+n]
  cat(nrow(sc1),"QTLs in",p,"\n")
  if(nrow(sc1)>1) intpPhen.bin<-c(intpPhen.bin,p)
  if(nrow(sc1)==1) effpPhen.bin<-c(effpPhen.bin,p)}
}

#sig and almost-sig qtls sum table
#normal
{qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",
                alpha=0.95,ci.function="bayesint",pvalues=T)
basil<-calc.genoprob(basil,2,map.function="kosambi")
qtlist.aq<-list()
s.aq<-list()
for (i in 1:n){
  p<-phenames(basil)[i]
  if(!is.null(out.aq[[p]])){
    s<-summary(out.aq[[p]],format="tabByCol",perms=scan1perm[,p],
               alpha=0.95,ci.function="bayesint",pvalues=T)
    if(nrow(s[[1]])>0){
      qtlist.aq[p]<-s
      qtlist.aq[[p]]<-cbind.data.frame(Trait=p,qtlist.aq[[p]])
      s.aq[[p]]<-summary(out.aq[[p]],perms=scan1perm[,p],alpha=0.63)
      if(nrow(s.aq[[p]])>0){
        rqtl[[p]]<-addtoqtl(basil,rqtl[[p]],s.aq[[p]][,1],s.aq[[p]][,2])
      }
    }
  }
}
for (i in 1:length(qtlist)){
  if(colnames(qtlist[[i]])[1]!="Trait"){
      qtlist[[i]]<-cbind.data.frame(Trait=names(qtlist[i]),qtlist[[i]])
  }
}
qtldf<-do.call(rbind.data.frame,c(qtlist,make.row.names=F))

for (i in 1:n){
  if(names(qtlist[i])%in%names(qtlist.aq)){
    qtlist[[i]]<-rbind(qtlist[[i]],qtlist.aq[[phenames(basil)[i]]])
  }
}
qtldf.aq<-do.call(rbind.data.frame,c(qtlist,make.row.names=F))

#binary
qtlist<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",
                alpha=0.95,ci.function="bayesint",pvalues=T)
basil<-calc.genoprob(basil,3,map.function="kosambi")
for (i in 1:2){
  p<-phenames(basil)[i+n]
  if(!is.null(out.aq.bin[[p]])){
    s<-summary(out.aq.bin[[p]],perms=scan1perm.bin[,p],alpha=0.95,
               format="tabByCol",ci.function="bayesint",pvalues=T)
    if(nrow(s[[1]])>0){
      qtlist.aq[p]<-s
      qtlist.aq[[p]]<-cbind.data.frame(Trait=p,qtlist.aq[[p]])
      s.aq[[p]]<-summary(out.aq.bin[[p]],perms=scan1perm.bin[,p],alpha=0.63)
      if(nrow(s.aq[[p]])>0){
        rqtl.bin[[p]]<-addtoqtl(basil,rqtl.bin[[p]],s.aq[[p]][,1],s.aq[[p]][,2])
      }
    }
  }
}
for (i in 1:length(qtlist)){
  qtlist[[i]]<-cbind.data.frame(Trait=names(qtlist[i]),qtlist[[i]])
}
for (i in 1:length(qtlist)){
  qtldf<-rbind.data.frame(qtldf,qtlist[[i]],make.row.names=F)
}
for (i in 1:2){
  if(names(qtlist[i])%in%names(qtlist.aq)){
    qtlist[[i]]<-rbind(qtlist[[i]],qtlist.aq[[phenames(basil)[i+n]]])
  }
}
qtldf.aq.bin<-do.call(rbind.data.frame,c(qtlist,make.row.names=F))
qtldf.aq<-rbind(qtldf.aq,qtldf.aq.bin)

#final df of qtls at alp= 0.99 and sig *** levels
su<-1-summary(basil)$missing.phe
qtldf.aq<-qtldf.aq%>%
  mutate("Len of LG"=round(chrlen(basil)[chr],1),.after=chr)%>%
  mutate("Len of QTL"=round(ci.high-ci.low,1),.after="Len of LG")%>%
  mutate("Flanking markers"=
           paste0(chr,"_m",find.marker(basil,chr,ci.low),"-",
                  chr,"_m",find.marker(basil,chr,ci.high)))%>%
  mutate("Central marker"=paste0(chr,"_m",find.marker(basil,chr,pos)))%>%
  mutate("Pval"=paste0(pval,if_else(pval<0.63,"*",""), 
                       if_else(pval<0.1,"*",""), if_else(pval<0.05,"*","")))%>%
  select(!c(pval,ci.low,ci.high))%>%
  rename("QTL's LG"=chr)%>%
  mutate("No. Inds/% phenotyped"=paste0(nind(basil)*su[find.pheno(basil,Trait)],
                                        " ind / ",round(100*su[find.pheno(basil,Trait)],1),"%"),
         .after=Trait)%>%
  mutate("No. QTLs"="",.before="QTL's LG")%>%
  mutate("pos"=round(pos,1))%>%
  mutate("lod"=round(lod,1))
}
write.csv(qtldf.aq,file="qtldf.aq.csv")

####Effectplot####
#Simulate genotypes given observed marker data 
basil<-sim.geno(basil,n.draws=128,step=2,map.function="kosambi")
#get the  qtl maximum pick for a threshold in pheno i
{tiff(filename="effectPlot.tiff",750,700)
  par(mfrow=c(2,3),cex=1.3,cex.lab=1.2)
  #normal
  for (i in find.pheno(basil,effpPhen)) {
    (sc1<-summary(scan1,perms=scan1perm,
                  alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
    (mn<-find.marker(basil,chr=sc1[,1],pos=sc1[,2]))
    (pos.nam<-paste0(sc1[,1],"_m",mn))
    effectplot(basil,pheno.col=i,mname1=mn,xlab=pos.nam,
               main="",ylab=phenames(basil)[i])
  }
  #binary
  if(length(effpPhen.bin)>0){
    for (i in 1:2){
      (sc1<-summary(scan1.bin,perms=scan1perm.bin,
                    alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
      (mn<-find.marker(basil,chr=sc1[,1],pos=sc1[,2]))
      (pos.nam<-paste0(sc1[,1],"_m",mn))
      effectplot(basil,pheno.col=i+n,mname1=mn,xlab=pos.nam,
                 main="",ylab=phenames(basil)[i+n])
    }
  }
  dev.off()
}

####Interaction plot####
#merged interaction plot
mergedIntp<-function(x=1,bin=F,sc=scan1,perm=scan1perm,
                     rit.inx=NULL,mf=c(2,3),l=1,LETTERs=T,second=F,first=F){
  rit<-1;if (!second) par(mfrow=mf);par(cex=1.1,cex.lab=1.4)
  for (i in x) {
    (sc1<-summary(sc,perm=perm,alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
    (r<-nrow(sc1))
    sorted<-sc1%>%arrange(desc(across(3)))
    (mn<-find.marker(basil,chr=sorted[,1],pos=sorted[,2]))
    for(j in 1:(r-1)){
      for(k in (j+1):r){
        if (bin) i<-i+n
        effectplot(basil,pheno.col=i,mname1=mn[j],mname2=mn[k],main="",
                   ylab=paste("Ave. phenotype:",phenames(basil)[i]),
                   xlab=paste0(sorted[k,1],"_m",mn[k]),add.legend = F)
        if(rit%in%rit.inx) lpos<-"topright" else lpos<-"topleft"
        legend(lpos,c("AA","AB","BB"),lty=1,pch=1,col=c("black","red","blue"),
               bty="n",inset=c(0.01,0))
        a<-par("usr")
        x.leg <- a[1] * 0.05 + a[2] * 0.75
        y.leg <- a[4] - diff(a[3:4]) * 0.05
        if(rit%in%rit.inx) tpos<-2 else tpos<-NULL
        text(x.leg,y.leg,pos=tpos,paste0(sorted[j,1],"_m",mn[j]))
        if (LETTERs)mtext(LETTERS[l],adj=0,cex=1.3);l<-l+1;if(l>6)l<-1
        rit<-rit+1}
    }
  }
  if (!first) par(mfrow=c(1,1))
}
#normal
{tiff("interaction plot%d.tiff",750,700)
x<-(1:6)[1:6%in%find.pheno(basil,intpPhen)]
mergedIntp(x,rit.inx=c(3,6,8,12,15,21,22,26,35))
x<-(7:11)[7:11%in%find.pheno(basil,intpPhen)]
mergedIntp(x,first=T,rit.inx=1,mf=c(2,2))
#binary
x<-((12:13)[12:13%in%find.pheno(basil,intpPhen.bin)])-n
mergedIntp(x[1],T,scan1.bin,scan1perm.bin,second=T,l=2)
mergedIntp(x[2],T,scan1.bin,scan1perm.bin,LETTERs=F)
dev.off()}

####Scantwo####
#normal
#scan genome for Two-QTL model
basil<-clean(basil)
basil<-calc.genoprob(basil,2,map.function="kosambi")
##scan2<-scantwo(basil,pheno.col=c(1:n),method="hk")
load(file="sc2.hk.Rdata")
#permutation test
load(file="sc2perm.hk.Rdata")

####Not run!
scan2perm<-scantwo(basil,pheno.col=c(1:n),method="hk",n.perm=1000)
operm2 <- vector("list", 100)
for(i in 1:100 ){
  operm2[[i]]<-scantwo(basil,pheno.col=c(1:6),method="hk",
                       n.perm=9,n.cluster=2)}
res1<-do.call("rbind", operm2)
#res2
operm2<-vector("list", 100)
for(i in 1:100 ){
  operm2[[i]]<-scantwo(basil,pheno.col=c(7:8),method="hk",
                       n.perm=9,n.cluster=2)}
res2<-do.call("rbind", operm2)
#res3
operm2 <- vector("list", 100)
for(i in 1:100 ){
  operm2[[i]]<-scantwo(basil,pheno.col=c(9:10),method="hk",
                       n.perm=9,n.cluster=2)}
res3<-do.call("rbind", operm2)
#res4
operm2 <- vector("list", 100)
for(i in 1:100 ){
  operm2[[i]]<-scantwo(basil,pheno.col=11,method="hk",
                       n.perm=9,n.cluster=2)}
res4<-do.call("rbind", operm2)
#cbind the 11 phenos 
scan2perm<-res1
for(i in 1:6){
  scan2perm[[i]]<-cbind(res1[[i]],res2[[i]],res3[[i]],res4[[i]])
}
save(scan2perm,file="sc2perm.hk.Rdata")
####

#scan2 results for every phenotype as a df
c.thr1<-list()
for(i in 1:n){
  (thr1<-summary(scan2, perms=scan2perm, alpha=0.2, 
                 lodcolumn=i,pvalues=T))
  if(i==1){c.thr1[[phenames(basil)[i]]]<-thr1
  }else c.thr1[[phenames(basil)[i]]]<-thr1
}
for (i in 1:length(c.thr1)){
  if(nrow(c.thr1[[i]])>0&&colnames(c.thr1[[i]])[1]!="Trait"){
    c.thr1[[i]]<-cbind.data.frame(Trait=names(c.thr1[i]),c.thr1[[i]])
  }
}
thr1df<-do.call(rbind.data.frame,c(c.thr1,make.row.names=F))

#binary
#genome scan for two-QTL
basil<-calc.genoprob(basil,10,map.function="kosambi")
##scan2.bin<-scantwo(basil,pheno.col=c(12,13),method="em",
##                   model="binary",verbose=T)
load("sc2.bin.em10.Rdata")
#permutation test
load("sc2perm.bin.em10.Rdata")

####Not run!
scan2perm.bin<-scantwo(basil,pheno.col=c(12:13),method="em",
                       model="binary",n.perm=1000)
operm2 <- vector("list", 200)
for(i in 1:100 ){
  operm2[[i]]<-scantwo(basil,pheno.col=c(12:13),method="em",
                       model="binary",n.perm=5,n.cluster=2)
}
for(i in 101:200 ){
  operm2[[i]]<-scantwo(basil,pheno.col=c(12:13),method="em",
                       model="binary",n.perm=5,n.cluster=2)
}
scan2perm.bin<-do.call("rbind", operm2[1:200])
save(scan2perm.bin,file="sc2perm.bin.em10.Rdata")
####

#scan2 results for every phenotype as a df
c.thr2<-list()
for(i in 1:2){
  (thr2<-summary(scan2.bin, perms=scan2perm.bin, alpha=0.2, 
                 lodcolumn=i,pvalues=T))
  if(i==1){c.thr2[[phenames(basil)[i+n]]]<-thr2
  }else c.thr2[[phenames(basil)[i+n]]]<-thr2
}
for (i in 1:length(c.thr2)){
  if(nrow(c.thr2[[i]])>0&&colnames(c.thr2[[i]])[1]!="Trait"){
    c.thr2[[i]]<-cbind.data.frame(Trait=names(c.thr2[i]),c.thr2[[i]])
  }
}
thr2df<-do.call(rbind.data.frame,c(c.thr2,make.row.names=F))
thr<-rbind(thr1df,thr2df)
#sc2thr1+2.hk+em10
write.csv(thr,"qtlpairs.csv")

####addpair+PVE####
#normal
{
basil<-calc.genoprob(basil,2,map.function="kosambi")
(qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",alpha=0.63,
                 ci.function="bayesint"))
(sc2thr1<-summary(scan2perm,alpha=0.2))
#rearrange the threshold list
th<-vector('list',5)
for(j in 1:5){th[[j]]<-t(sc2thr1[[j]])}
m<-do.call('cbind',th)
dimnames(m)<-list(phenames(basil)[1:n],names(sc2thr1)[1:5])
out.ap<-list();qtlpairs<-list();s.fq<-list();out.fq<-list()
for(i in 1:n){
  p<-phenames(basil)[i]
  if(length(qtlist[[p]][,1])>0){
    out.ap[[p]]<-addpair(basil,qtlist[[p]][,1],p,rqtl[[p]],
                         method="hk",verbose=T)
    qtlpairs[[p]]<-summary(out.ap[[p]],thresholds=m[p,])
    s.fq[[p]]<-summary(out.fq[[p]]<- fitqtl(basil,p,rqtl[[p]],
                                            method="hk",get.ests=T))
  }
}
}
#binary
{
basil<-calc.genoprob(basil,10,map.function="kosambi")
(qtlist.bin<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",
                 alpha=0.63,ci.function="bayesint"))
(sc2thr2<-summary(scan2perm.bin,alpha=0.2))
#rearrange the threshold list
th<-vector('list',5)
for(j in 1:5){th[[j]]<-t(sc2thr2[[j]])}
m.bin<-do.call('cbind',th)
dimnames(m.bin)<-list(phenames(basil)[n+1:2],names(sc2thr2)[1:5])
out.ap.bin<-list();rqtl2.bin<-list()
for(i in 1:2){
  p<-phenames(basil)[i+n]
  if(length(qtlist.bin[[p]][,1])>0){
    q<-rbind(qtlist.bin[[p]][,-c(3,4)],s.aq[[p]])
    rqtl2.bin[[p]]<-refineqtl(basil,p,makeqtl(basil,q[,1],q[,2],
      what="prob"),maxit.fitqtl=1e+6,tol=0.05,method="hk",model="binary")
    out.ap.bin[[p]]<-addpair(basil,q[,1],p,rqtl2.bin[[p]],maxit=1e+6,
                             tol=0.2,method="hk",model="binary",verbose=T)
    qtlpairs[[p]]<-summary(out.ap.bin[[p]],thresholds=m.bin[p,])
    s.fq[[p]]<-summary(out.fq[[p]]<- fitqtl(basil,p,rqtl2.bin[[p]],
              maxit=1e+6,tol=0.01,method="hk",model="binary",get.ests=T))
  }
}
}

for (i in 1:length(qtlpairs)){
  if(nrow(qtlpairs[[i]])>0 && names(qtlpairs[[i]])[1]!="Trait"){
    qtlpairs[[i]]<-cbind.data.frame(Trait=names(qtlpairs[i]),qtlpairs[[i]])
  }
}
qtlpairsdf<-do.call(rbind.data.frame,c(qtlpairs,make.row.names=F))
#interacting QTL detected for AFL
qtlpairsdf<-c(qtlpairsdf,thr=m[qtlpairsdf[,1],])

#sc2thr1+2.hk+em10.ap.hk
write.csv(qtlpairsdf,"qtlpairs.ap.Alp0.2.csv")

#plot(scan2,contours=T,upper="int",lower="add",lodcolumn=i,
#     alternate.chrid=T)
#mtext("int",adj=-2)
#mtext("add",side=4,line=-4,las=1)

save.image("denovo.Rdata")