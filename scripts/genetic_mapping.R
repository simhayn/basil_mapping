#*****************************
#Genetic map construction of basil
#(BA17 Bc5F3 44-4-3 x Deep Purple)F2 population 
#Using r/qtl and r/asmap packages 
#SNP's from Lara
#GBS sequenced by Genewiz
#Nataly Yakobov
#*****************************
#       R=0       S=1      for FOB BDM 
#       R=1       S=9      for Cold

#recommended to source, added breakpoints where the user is asked whether to run the current analysis or plot

#clear workspace
rm(list=ls())
#install and load packages
sapply(c("qtl", "ASMap"), function(pkg) if (!requireNamespace(pkg, quietly = T)) { install.packages(pkg); library(pkg, character.only = T) })

#import genotypic data #warning message: long chromosomes. It is ok cause all markers are on 1 chr 
data<-read.cross("csvr",".","bdenovo.csv",genotypes=c("a","h","b"),map.function="kosambi")

#convert cross to bcsft type f2
data<-convert2bcsft(data,F.gen=2,estimate.map=F)
cat("Phenos:",head(phenames(data),-1),"\n\n\n")
cat("\n\n\n---Data Pre-processing\n\n")
#look at the pattern of missing data. 
plotMissing(data)
cat("\nBlack pixels indicate missing genotypes\n\n\n")

#omit the individuals with more than 672 missing genotype data, 85% of total markers(791)
{
  cat("Briefly estimating map\n")
  data<-quickEst(data,map.function="kosambi")
  sg<-statGen(data,bychr=F,stat.type="miss",id='index')
  data1<-subset(data,ind=om<-sg$miss<(totmar(data))*0.85)
  cat("\n",if(F%in%om) which(!om)else "no","ind omitted for missing > 85% mar\n\n\n")
}

#plot the number of genotyped markers for each individual aside the number of genotyped individuals for each marker
{
  par(mfrow=c(1,2), las=1)
  plot(ntyped(data1), ylim=c(0,totmar(data1)+50),ylab="No. typed markers",main="Typed Markers by Individual")
  mtext("A",adj=0)
  plot(ntyped(data1, "mar"), ylim=c(0,(nind(data1)+15)),ylab="No. typed individuals",main="Typed Individuals by Marker")
  mtext("B",adj=0)
}

#plot the genotype frequencies by individual
{
  g <- pull.geno(data1)
  gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
  gfreq <- t(t(gfreq) / colSums(gfreq))
  par(mfrow=c(1,3), las=1)
  for(i in 1:3){
    plot(gfreq[i,], ylab="Genotype frequency",ylim=c(0,1))
    abline(h=mean(gfreq[i,]),lty=3,col="red",lwd=3)
    mtext(c("AA", "AB", "BB")[i])
  }
  par(mfrow=c(1,1));title("Genotypes' Frequency and Segregation Ratio",line = 2.5)
}

#compare the genotypes for all pairs of individuals
{
  cg<-comparegeno(data1);cgr<-cg[lower.tri((cg))]
  hist(cgr, breaks=seq(0, 1, len=101),xlab="No. matching genotypes", main="Matching Pairs of Individuals")
  rug(cgr)
  #mark the outlier with a red arrow 
  x<-max(cgr)
  arrow.plot(x,50,0,-1,true.angle = T,arrow.ex=30, length=.1,col='red', lwd=2)
  text(x,70,paste0(round(100*x,1),"%"),adj=c(0.5,0.2))
}

#omit the ind with more than 90% identical markers 
{
  wh<-which(cg>0.9,arr=T)
  data2<-subset(data1,ind=-wh[,2])
  cat("\n",nind(data1)-nind(data2),"ind omitted for",paste0(round(100*x,1),"%"),"identical geno\n",paste0("#", data1$pheno$index[wh[,2]]),"\n\n\n")
}

#pull out markers from cross
if(ask("Press <RETURN> to pull out markers from cross")==""){
totmar(data2)
data3<-pullCross(data2,type="missing",pars=list(miss.thresh=0.1))
totmar(data3)
cat(totmar(data2)-totmar(data3),"mar pulled for missing\n")

data4<-pullCross(data3,type="seg.distortion",pars=list(seg.thresh=0.001))
totmar(data4)
cat(totmar(data3)-totmar(data4),"mar pulled for seg. distortion\n")

data5<-pullCross(data4,type="co.located")
totmar(data5)
cat(totmar(data4)-totmar(data5),"mar pulled for co. located\n\n\n")
}
#plot p. value options to determine threshold for marker clustering
cat(nind(data5),"individuals are used to cluster markers\n\n\n")
{
  pValue(dist=seq(25,40,by=5),pop.size=110:190)
  cat("\nI chose pValue= 1e-7 (on y axis), meaning->> split linkage groups with more than 30cM gap between markers, according to 150 ind population (on x axis)\n\n\n")
}
#LOD(Logarithm of the Odds)= statistical measure of the likelihood that two loci (positions on a chromosome) are linked and therefore inherited together, rather than assorting independently
#form linkage groups with LOD=7
{
  data5<-mstmap(data5,bychr=F,p.value=1e-7,id='index')
  plotMap(data5,alternate.chrid = T)
  
  #profile ind genotype statistics
  pg<-profileGen(data5,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=20,layout=c(1,3),lty=2,cex=0.7)
  pg<-profileGen(data5,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg$stat$xo),layout=c(1,3),lty=2,cex=0.8)
  
  #omit missing/ double crossover(dxo)/ xo statistics outlier individuals. 
  data6<-subsetCross(data5,ind=!pg$xo.lambda)
  cat("\n\n",nind(data5)-nind(data6),"ind omitted by profileGen\n\n\n")
  
  #double check dxo. An unusually high rate of double crossovers might indicate genotyping errors
  pg1<-profileGen(data6,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg$stat$xo),layout=c(1,3),lty=2,cex=0.7)
  pg1<-profileGen(data6,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg1$stat$xo),layout=c(1,3),lty=2,cex=0.7)
  
  #re-construct map.The genotyping errors can distort the distances between markers, the order of the inputted markers is respected though
  data7<-mstmap(data6,bychr=F,anchor=T,p.value=1e-7,id='index')
  cat("\n\n")
}
  #push back markers to the map
  if(ask("Press <RETURN> to push back markers")==""){
  totmar(data7)
  data8<-pushCross(data7,type="co.located")
  totmar(data8)
  cat("\n",totmar(data8)-totmar(data7),"mar pushed for co. located\n\n\n")
  
  #re-construct final map by adding markers to existing LGs(linkage groups)
  data9<-mstmap(data8,anchor=T,p.value=2,id='index')
  cat("\n\n")
  #drop LGs with less than 2 markers
  mndrop<-markernames(data9,nmar(data9)<2)
  data10<-drop.markers(data9,as.character(mndrop))
  totmar(data10)
  cat(totmar(data9)-totmar(data10),"mar omitted for LG < 2 mar\n\n\n")
}
#rename chr by numerical order
{
  x<-1:nchr(data10)
  for (i in x) {
    names(data10$geno)[i]<-paste0("LG",i)
  }
  tiff("plotMap.tiff",800,550)
  par(cex=1.1,cex.lab=1.4,cex.main=2)
  plot.map(data10,alternate.chrid=T)
  dev.off()
}

#LGs summary table
LGtab<-summaryMap(data10)
if(ask("Press <RETURN> to save map summary to csv file")=="")
  write.csv(LGtab,"summaryMap.csv")

#heatmap of LOD and Rf
cat("\nestimating recombination fraction\n\n")
data10<-est.rf(data10)
if(ask("Press <RETURN> to plot rf&LOD heatmap to tiff file, it takes a long time")==""){
  tiff("heatMap.tiff",1050,590,res=105)
  heatMap(data10,lmax=50,main="")
  mtext("Pairwise Recombination Fractions and LOD Scores",cex=1.8,line=2.5,adj=0.4,font=2)
  dev.off()
}

#plot Pairwise LOD vs. Rf
rf<-pull.rf(data10);lod<-pull.rf(data10,what="lod")
if(ask("Press <RETURN> to plot LOD over rf to tiff file, it takes a long time")==""){
  tiff("LODvsRf.tiff",800,550);par(cex=1.5)
  plot(as.numeric(rf),as.numeric(lod),xlab="Recombination fraction",ylab="LOD score",main=paste("Pairwise LOD vs. Rf for",totmar(data10),"Markers"))
  dev.off()
}
#save map as csv file 
cat("\nsaving cross 'data10' to csv file\n")
write.cross(data10,"csv")
