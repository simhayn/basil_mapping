#*****************************
#Genetic map construction of basil
#(BA17 Bc5F3 44-4-3 x Deep Purple)F2 population 
#Using r/qtl and r/asmap packages 
#SNP's from Lara
#GBS sequenced by Genewiz
#Nataly Simhayev
#*****************************
#       R=0       S=1      for FOB BDM 
#       R=1       S=9      for Cold

# Install and load qtl package
if (!requireNamespace("qtl", quietly = TRUE)) install.packages("qtl"); library(qtl)
# Install and load ASMap package
if (!requireNamespace("ASMap", quietly = TRUE)) install.packages("ASMap"); library(ASMap)


rm(list=ls())
#load("C:/Users/Public/R nataly/denovo/denovo.Rdata")

#import .csv file named "cknat.csv" #warning message: long chromosomes. It is ok cause all markers are on 1 chr 
cknat<-read.cross("csvr",".","bdenovo.csv",genotypes=c("a","h","b"),
                  map.function="kosambi")
#convert cross to bcsft type f2
cknat<-convert2bcsft(cknat,F.gen=2,estimate.map=F)
phenames(cknat)

#look at the pattern of missing data. Black pixels indicate missing genotypes
{tiff("plotMissing.tiff")
par(cex=1.3)
plotMissing(cknat)
dev.off()}
#omit the individuals with more than 672 missing genotype data, 85% of total markers(791)
{cknat<-quickEst(cknat,map.function="kosambi")
sg<-statGen(cknat,bychr=F,stat.type="miss",id='index')
cknat1<-subset(cknat,ind=sg$miss<(totmar(cknat))*0.85)
cat(which(!sg$miss<(totmar(cknat))*0.85) ,"ind omitted for miss > 85%")}

#plot the number of genotyped markers for each individual aside the number of genotyped individuals for each marker
{tiff("ntyped.tiff",900,400,res=100)
par(mfrow=c(1,2), las=1,cex=1.3)
plot(ntyped(cknat1), ylim=c(0,totmar(cknat1)+50),
     ylab="No. typed markers",main="Typed Markers by Individual")
mtext("A",adj=0)
plot(ntyped(cknat1, "mar"), ylim=c(0,(nind(cknat1)+15)),
     ylab="No. typed individuals",main="Typed Individuals by Marker")
mtext("B",adj=0)
dev.off()}

#plot the genotype frequencies by individual
{g <- pull.geno(cknat1)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))}
{tiff("genoByind.tiff",900,500,res=100)
par(mfrow=c(1,3), las=1,cex=1.4)
for(i in 1:3){
  plot(gfreq[i,], ylab="Genotype frequency", 
       main="",ylim=c(0,1))
  abline(h=mean(gfreq[i,]),lty=3,col="red",lwd=3)
  mtext(c("AA", "AB", "BB")[i],cex=1.4)
}
par(mfrow=c(1,1))
title("Genotypes' Frequency and Segregation Ratio",cex.main=1.7)
dev.off()}

#compare the genotypes for all pairs of individuals
cg<-comparegeno(cknat1)
{tiff("matchingGeno.tiff",700,500)
par(cex=1.3)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101),
     xlab="No. matching genotypes", main="Matching Pairs of Individuals")
rug(cg[lower.tri(cg)])
#mark the outlier with a red arrow 
x<-max(cg[lower.tri(cg)])
arrow.plot(x,50,0,-1,true.angle = T,arrow.ex=30, length=.1,
           col='red', lwd=2)
text(x,70,paste0(round(100*x,1),"%"))
dev.off()}

#omit the ind with higher than 90% identical markers 
wh<-which(cg>0.9,arr=T)
cknat2<-subset(cknat1,ind=-wh[,2])
cat(nind(cknat1)-nind(cknat2),"ind omitted for",paste0(round(100*x,1),
    "%"),"identical geno\n",cknat1$pheno$index[wh[,2]],"\n")

#pull out markers from cross
totmar(cknat2)
cknat3<-pullCross(cknat2,type="missing",
                  pars=list(miss.thresh=0.1))
totmar(cknat3)
cat(totmar(cknat2)-totmar(cknat3),"mar pulled for missing")

cknat4<-pullCross(cknat3,type="seg.distortion",
                  pars=list(seg.thresh=0.001))
totmar(cknat4)
cat(totmar(cknat3)-totmar(cknat4),"mar pulled for seg. distortion")

cknat5<-pullCross(cknat4,type="co.located")
totmar(cknat5)
cat(totmar(cknat4)-totmar(cknat5),"mar pulled for co. located")

#plot p. value options to determine threshold for marker clustering
nind(cknat5)
tiff("pvalue.tiff",690,520,res=96)
pValue(dist=seq(25,40,by=5),pop.size=110:190)
dev.off()

#form linkage groups with LOD=7
cknat5<-mstmap(cknat5,bychr=F,p.value=1e-7,id='index')
plotMap(cknat5,alternate.chrid = T)

#profile ind genotype statistics
pg<-profileGen(cknat5,bychr=F,stat.type=c("xo","dxo","miss"),
               id="index",xo.lambda=20,layout=c(1,3),lty=2,cex=0.7)
dev.off()
tiff("profileGen.tiff",700,620,res=100)
pg<-profileGen(cknat5,bychr=F,stat.type=c("xo","dxo","miss"),
               id="index",xo.lambda=median(pg$stat$xo),layout=c(1,3),
               lty=2,cex=0.8)
dev.off()

#omit missing/ double xo/ xo statistics outlier individuals
cknat6<-subsetCross(cknat5,ind=!pg$xo.lambda)
cat(nind(cknat5)-nind(cknat6),"ind omitted by profileGen")

#double check DXO
pg1<-profileGen(cknat6,bychr=F,stat.type=c("xo","dxo","miss"),
               id="index",xo.lambda=median(pg$stat$xo),layout=c(1,3),lty=2,cex=0.7)
pg1<-profileGen(cknat6,bychr=F,stat.type=c("xo","dxo","miss"),
                id="index",xo.lambda=median(pg1$stat$xo),layout=c(1,3),lty=2,cex=0.7)

#re-construct map
cknat7<-mstmap(cknat6,bychr=F,anchor=T,p.value=1e-7,id='index')

#push back markers to the map
totmar(cknat7)
cknat8<-pushCross(cknat7,type="co.located")
totmar(cknat8)
cat(totmar(cknat8)-totmar(cknat7),"mar pushed for co. located")

#re-construct final map by adding markers to existing LGs
cknat9<-mstmap(cknat8,anchor=T,p.value=2,id='index')

#drop LGs with less than 2 markers
mndrop<-markernames(cknat9,nmar(cknat9)<2)
basil<-drop.markers(cknat9,as.character(mndrop))
totmar(basil)
cat(totmar(cknat9)-totmar(basil),"mar omitted for LG < 2")

#rename chr by numerical order
{x<-1:nchr(basil)
for (i in x) {
  names(basil$geno)[i]<-paste0("LG",i)
}
tiff("plotMap.tiff",800,550)
par(cex=1.1,cex.lab=1.4,cex.main=2)
plot.map(basil,alternate.chrid=T)
dev.off()
}

#LGs summary table
LGtab<-summaryMap(basil)
write.csv(LGtab,"summaryMap.csv")

#heatmap of LOD and Rf
basil<-est.rf(basil)
{tiff("heatMap.tiff",1050,590,res=105)
heatMap(basil,lmax=50,main="")
mtext("Pairwise Recombination Fractions and LOD Scores",
      cex=1.8,line=2.5,adj=0.4,font=2)
dev.off()}

#plot Pairwise LOD vs. Rf
rf<-pull.rf(basil)
lod<-pull.rf(basil,what="lod")
tiff("LODvsRf.tiff",800,550)
par(cex=1.5)
plot(as.numeric(rf),as.numeric(lod),xlab="Recombination fraction",
     ylab="LOD score",main=paste("Pairwise LOD vs. Rf for",
                                 totmar(basil),"Markers"))
dev.off()

#save map as csv file 
write.cross(basil,"csv","basil")
save.image("denovo.Rdata")