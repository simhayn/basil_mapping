---
title: "Genetic Mapping and QTL Analysis"
output: github_document
---

<img src="https://octodex.github.com/images/labtocat.png" width="200" height="200"/>

Load packages

```{r}
library(qtl)
library(qtlcharts)
library(ASMap)
library(dplyr)
library(corrplot)
library(LinkageMapView)
```

Import genotypic data. for now, all markers are on 1 chr

```{r}
data<-read.cross("csvr",".","raw_data.csv",genotypes=c("a","h","b"),map.function="kosambi")
```

Convert cross to bcsft type f2 to use ASMap package

```{r}
data<-convert2bcsft(data,F.gen=2,estimate.map=F)
cat("Phenos:",(head(phenames(data),-1)))
n<-11 #normal phenotypes count
```

# Data Pre-processing

Look at the pattern of missing data. Black pixels indicate missing genotypes.

```{r}
plotMissing(data)
```

Omit the individuals with more than 85% of total markers -\> less than 672 markers. Briefly estimate map first.

```{r}
data<-quickEst(data,map.function="kosambi")
sg<-statGen(data,bychr=F,stat.type="miss",id='index')
data1<-subset(data,ind=om<-sg$miss<(totmar(data))*0.85)
cat(if(F%in%om) which(!om)else "no","ind omitted for missing > 85% mar")
```

Plot the number of genotyped markers per individual as well as the number of genotyped individuals per marker

```{r}
par(mfrow=c(1,2), las=1)
plot(ntyped(data1), ylim=c(0,totmar(data1)+50),ylab="No. typed markers",main="Markers by Individual")
mtext("A",adj=0)
plot(ntyped(data1, "mar"), ylim=c(0,(nind(data1)+15)),ylab="No. typed individuals",main="Individuals by Marker")
mtext("B",adj=0)
```

Plot the genotype frequencies per individual

```{r}
g <- pull.geno(data1)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3){
  plot(gfreq[i,], ylab="Genotype frequency",ylim=c(0,1))
  abline(h=mean(gfreq[i,]),lty=3,col="red",lwd=3)
  mtext(c("AA", "AB", "BB")[i])}
par(mfrow=c(1,1));title("Genotypes' Frequency and Segregation Ratio",line = 2.5)
```

Compare the genotypes for all pairs of individuals

```{r}
cg<-comparegeno(data1);cgr<-cg[lower.tri((cg))]
hist(cgr, breaks=seq(0, 1, len=101),xlab="No. matching genotypes", main="Matching Pairs of Individuals")
rug(cgr)
#mark the outlier with a red arrow 
x<-max(cgr)
arrow.plot(x,50,0,-1,true.angle = T,arrow.ex=30, length=.1,col='red', lwd=2)
text(x,70,paste0(round(100*x,1),"%"),adj=c(0.5,0.2))
```

Omit individuals with more than 90% identical markers

```{r}
wh<-which(cg>0.9,arr=T)
data2<-subset(data1,ind=-wh[,2])
cat(nind(data1)-nind(data2),"ind omitted for",paste0(round(100*x,1),"%"),"identical geno\n",paste0("#", data1$pheno$index[wh[,2]]))
```

Pull out markers from cross temporarily

```{r}
cat(totmar(data2),'total mar\n')
data3<-pullCross(data2,type="missing",pars=list(miss.thresh=0.1))
cat(totmar(data3),'total mar\n')
cat(totmar(data2)-totmar(data3),"mar pulled for missing\n")

data4<-pullCross(data3,type="seg.distortion",pars=list(seg.thresh=0.001))
cat(totmar(data4),'total mar\n')
cat(totmar(data3)-totmar(data4),"mar pulled for seg. distortion\n")

data5<-pullCross(data4,type="co.located")
cat(totmar(data5),'total mar\n')
cat(totmar(data4)-totmar(data5),"mar pulled for co. located")
```

Plot optional p. values to determine distance threshold for marker clustering. I chose pValue= 1e-7 (on y axis), meaning-\>\> split linkage groups with more than 30cM gap between markers, according to 150 ind population (on x axis)

```{r}
cat(nind(data5),"individuals are used to cluster markers")
pValue(dist=seq(25,40,by=5),pop.size=110:190)
```

LOD(Logarithm of the Odds)= statistical measure of the likelihood that two loci (positions on a chromosome) are linked and therefore inherited together, rather than assorting independently.

# Map Construction

Form linkage groups with LOD=7, raw map

```{r}
data5<-mstmap(data5,bychr=F,p.value=1e-7,id='index')
plotMap(data5,alternate.chrid = T)
```

Profile individuals' genotype statistics

```{r include=FALSE}
pg<-profileGen(data5,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=20,layout=c(1,3),lty=2,cex=0.7)
pg<-profileGen(data5,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg$stat$xo),layout=c(1,3),lty=2,cex=0.8)
```

Omit missing/ double crossover (dxo)/ xo statistics outlier individuals

```{r}
data6<-subsetCross(data5,ind=!pg$xo.lambda)
cat(nind(data5)-nind(data6),"ind omitted by profileGen")
```

Double check dxo. An unusually high rate of double crossovers might indicate genotyping errors

```{r}
pg1<-profileGen(data6,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg$stat$xo),layout=c(1,3),lty=2,cex=0.7)
pg1<-profileGen(data6,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg1$stat$xo),layout=c(1,3),lty=2,cex=0.7)
```

Re-construct map. The genotyping errors can distort the distances between markers, the order of the inputted markers is respected

```{r}
data7<-mstmap(data6,bychr=F,anchor=T,p.value=1e-7,id='index')
```

Push back markers to the map

```{r}
cat(totmar(data7),'total mar\n ')
data8<-pushCross(data7,type="co.located")
cat(totmar(data8),'total mar\n ')
cat(totmar(data8)-totmar(data7),"mar pushed for co. located")
```

Re-construct final map by adding markers to existing LGs (linkage groups)

```{r include=FALSE}
data9<-mstmap(data8,anchor=T,p.value=2,id='index')
```

Drop LGs with less than 2 markers

```{r}
mndrop<-markernames(data9,nmar(data9)<2)
data10<-drop.markers(data9,as.character(mndrop))
cat(totmar(data10),'total mar\n')
cat(totmar(data9)-totmar(data10),"mar omitted for LG < 2 mar")
```

Rename chr by numerical order

```{r}
x<-1:nchr(data10)
for (i in x) {names(data10$geno)[i]<-paste0("LG",i)}
```

Plot genetic map illustration, final map

```{r}
plot.map(data10,alternate.chrid=T)
```

LGs summary table

```{r}
summaryMap(data10)
```

# Map Evaluation

Estimate recombination fraction

```{r}
data10<-est.rf(data10)
```

Heatmap of LOD and Rf

```{r}
heatMap(data10,lmax = 50,main='')
mtext("Pairwise Recombination Fractions and LOD Scores",cex=1.1,line=3.2,adj=0.4,font=2)
```

Plot Pairwise LOD vs. Rf

```{r}
rf<-pull.rf(data10);lod<-pull.rf(data10,what="lod")
plot(as.numeric(rf),as.numeric(lod),xlab="Recombination fraction",ylab="LOD score",main=paste("Pairwise LOD vs. Rf for",totmar(data10),"Markers"))
```

The evaluation looks good:

The heatmap is continuous with low gradient along the chr meaning that the markers within a chromosome are gradually distant from each other.

On the LOD/rf scatter plot there is a trend line with no outliers.

# Exploratory Data Analysis

Jitter map- to avoid marker overlaping by slightly adding gaps between them

```{r}
data<-jittermap(data10)
```

Create a df of parents phenotype

```{r}
cold_values <- c(5, 5, 4, 6, 7, 9, 5, 6, 6, 7, NA, NA, 7, 8, 7, 8, 8, 8, 9, 9, 9, 8, NA, NA)
faudpc_values <- c(rep(0, 12), rep(c(90, 142), each = 6))
parents <- data.frame(rbind(matrix(0,12,9), matrix(5,12,9)),faudpc_values, Cold = cold_values,row.names = c(paste0("P1_", 1:12), paste0("DP_", 1:12)))
ant_values <-c(5, 4.5, 4.5, 4.5, 4.5, 3.5, 4, 4.5, 4, 4.5, 4, 4)
parents[13:24,1:3] <-matrix(ant_values,12,3)
colnames(parents) <-phenames(data)[1:11]
head(parents)
```

Check the phenotypic distribution

Plot histograms or barplots

```{r}
cbind(phenames(data)) 
```

```{r}
par(mfrow = c(2, 3))
l <- 1
for (i in c(1:10, 12, 15, 13, 11)) {
  if (i < 12) {
    plotPheno(data, i, ylab = "Frequency", xlab = if (i < 7) "Purple intensity" else "Disease intensity")}
  if (i == 15) {
    par(mfrow = c(2, 3))}
  if (i > 11 & i != 15) {
    BP <- plotPheno(data, i, ylab = "Frequency",xlab = "Disease resistance",names.arg = c("R", "S"))
    text(BP, table(data$pheno[, i]),labels = table(data$pheno[, i]),pos = 1)
    points(c(0.7, 1.9), rep(1.6, 2), pch = 25,cex = 1.3,bg = c("palegreen1", "violet"))}
   if (i < 10) {
     points(if (i < 4) 0.1 else 0.7, 1.6, pch = 25, cex = 1.3, bg = "palegreen1")}
   if (i < 4) {
     points(mean(parents[13:24, i]) - 0.4, 1.6, pch = 25, cex = 1.3, bg = "violet")}
   if (i == 10) {
    points(c(mean(parents[1:12, i]) - 0.1, mean(parents[13:24, i])), rep(1.6, 2), pch = 25, cex = 1.3, bg = c("palegreen1", "violet"))}
  if (i == 11) {
    points(c(mean(parents[1:12, i], na.rm = TRUE), mean(parents[13:24, i], na.rm = TRUE)), rep(1.6, 2), pch = 25, cex = 1.3, bg = c("palegreen1", "violet"))}
  if (i %in% c(4, 5, 6, 9)) {
    points(if (i == 4) 6.7 else if (i == 9) 6.7 else if (i == 5 || i == 6) 4.3, 1.6, pch = 25, cex = 1.3, bg = "violet")}
   if (i %in% c(7, 8)) {
     points(5.5, 1.6, pch = 25, cex = 1.3, bg = "violet")}
   if (l > 6) l <- 1
   if (i %in% c(1:10, 12)) {
     if (i == 7) l <- 1
     mtext(LETTERS[l], adj = -0.1, cex = 1.1)}
   l <- l + 1
}
```

Plot the correlation matrices using corrplot

```{r}
mat<-cor(pull.pheno(data,1:13),use = "complete.obs")
corrplot(mat, method = "color", addCoef.col = "orange", tl.col = "black", tl.srt = 35,tl.cex=0.7,number.cex=0.5) 
```

Now, for each group of phenotypes

Fusarium

```{r}
mat<-cor(pull.pheno(data,c(7:10,12)),use = "complete.obs")
corrplot(mat, method = "color", addCoef.col = "white", tl.col = "black", tl.srt = 35)
```

Anthocyanin

```{r}
mat<-cor(pull.pheno(data,1:6),use = "complete.obs")
corrplot(mat, method = "color", addCoef.col = "white", tl.col = "black", tl.srt = 35)
```

There is a strong correlation between the phenotypes in the same group. the different conditions (in anthocyanin) or repetitions (in fusarium) had a limited effect on the resulting phenotype.

# QTL Analysis

Perform genome scans to identify QTL

```{r}
if (all(file.exists("genome_scans.Rdata","scan2_part1.Rdata","scan2_lod_part2.Rdata"))) {
  message("Loading precomputed genome scans from genome_scans.Rdata to save time.")
  load("genome_scans.Rdata")
  # Load the parts (scan2 is a large file so i split it to 2 parts)
  load("scan2_part1.Rdata")
  load("scan2_lod_part2.Rdata")
  lod_part1<-scan2.partial$lod
  # Combine the 'lod' matrix parts back into a single matrix
  combined_lod <- array(NA, dim = c(dim(lod_part1)[1] + dim(lod_part2)[1], dim(lod_part1)[2], dim(lod_part1)[3]))
  # Copy the data from the parts into the combined array
  combined_lod[1:dim(lod_part1)[1], , ] <- lod_part1
  combined_lod[(dim(lod_part1)[1] + 1):dim(combined_lod)[1], , ] <- lod_part2
  # Set the dimnames attribute correctly
  dimnames(combined_lod) <- dimnames(scan2.partial$lod)
  # Reconstruct the scantwo object
  scan2 <- scan2.partial  # Copy the original object structure
  scan2$lod <- combined_lod  # Replace the 'lod' part with the combined matrix
  # cleanup
  rm(list=c("scan2.partial","combined_lod","lod_part1","lod_part2"))
} else {
  cat("The precomputed genome scans were not found\n")
  message("Performing genome scans... This may take a long time.")
  
  ###ScanOne###
  #calculate probabilities (necessary for scanone)
  #normal #Haley-knott regression 
  data<-calc.genoprob(data,2,map.function="kosambi")
  #genome scan for Single-QTL
  scan1<-scanone(data,pheno.col=c(1:n),method="hk")
  #permutation test
  operm <- vector("list", 100)
  for(i in 1:100){operm[[i]]<-scanone(data,pheno.col=c(1:n),method="hk",n.perm=9,n.cluster=2)}
  scan1perm<-do.call("rbind", operm)
  
  #binary #EM algorithm maximum likelihood
  data<-calc.genoprob(data,3,map.function="kosambi")
  #genome scan for Single-QTL
  scan1.bin<-scanone(data,pheno.col=c(12:13),method="em",model="binary")
  #permutation test
  scan1perm.bin1<-scanone(data,pheno.col=c(12:13),method="em",model="binary",n.perm=1000,n.cluster=2)
  operm <- vector("list", 100)
  for(i in 1:100){operm[[i]]<-scanone(data,pheno.col=c(12:13),method="em",model="binary",n.perm=10,n.cluster=2)}
   scan1perm.bin1<-do.call("rbind", operm)
   
   save(scan1,scan1.bin,scan1perm,scan1perm.bin,scan2.bin,scan2perm,scan2perm.bin,file="genome_scans.Rdata")
   
   ###ScanTwo###
   #normal #scan genome for Two-QTL model
   data<-clean(data)
   data<-calc.genoprob(data,2,map.function="kosambi")
   scan2<-scantwo(data,pheno.col=c(1:n),method="hk")
   #permutation test
   ##this proccess is heavy so I've split it into 4 parts and then merged the results. with better CPU performance U can run the following line instead:
   ##scan2perm.bin<-scantwo(data,pheno.col=c(1:n),method="hk",n.perm=1000)
   #res1
   operm2 <- vector("list", 100)
   for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=c(1:6),method="hk",n.perm=9,n.cluster=2)}
   res1<-do.call("rbind", operm2)
   #res2
   operm2 <-vector("list", 100)
   for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=c(7:8),method="hk",n.perm=9,n.cluster=2)}
   res2<-do.call("rbind", operm2)
   #res3
   operm2 <- vector("list", 100)
   for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=c(9:10),method="hk",n.perm=9,n.cluster=2)}
   res3<-do.call("rbind", operm2)
   #res4
   operm2 <- vector("list", 100)
   for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=11,method="hk",n.perm=9,n.cluster=2)}
   res4<-do.call("rbind", operm2)
   #cbind the 11 phenos 
   scan2perm<-res1
   for(i in 1:6){scan2perm[[i]]<-cbind(res1[[i]],res2[[i]],res3[[i]],res4[[i]])}
   
   #binary #genome scan for two-QTL
   data<-calc.genoprob(data,10,map.function="kosambi")
   scan2.bin<-scantwo(data,pheno.col=c(12,13),method="em",model="binary",verbose=T)
   #permutation test
   ##this proccess is heavy so I've split it into 4 parts and then merged the results. with better CPU performance U can run the following line instead:
  ##scan2perm.bin<-scantwo(data,pheno.col=c(12:13),method="em",model="binary",n.perm=1000)
  operm2 <- vector("list", 200)
  for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=c(12:13),method="em",model="binary",n.perm=5,n.cluster=2)}
   for(i in 101:200 ){operm2[[i]]<-scantwo(data,pheno.col=c(12:13),method="em", model="binary",n.perm=5,n.cluster=2)}
   scan2perm.bin<-do.call("rbind", operm2[1:200])
}
```

Setting a QTL detection threshold according to permutation tests. The lower the precentage (5%), the better significance of the QTL.

```{r}
(thresh1.hk<-summary(scan1perm,alpha=c(0.63,0.1,0.05))) 
(thresh2.em<-summary(scan1perm.bin,alpha=c(0.63,0.1,0.05)))
```

Scan for additional QTL after reducing the masking effect of QTL with major peaks.

```{r}
#normal #add additional QTL
qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",alpha=0.63,ci.function="bayesint")
message("Calculating genoprob step= 2 cM")
data<-calc.genoprob(data,2,map.function="kosambi")
out.aq<-list();rqtl<-list()
for(i in 1:n){
  if(length(qtlist[[i]][,1])>0){
    p<-phenames(data)[i]
    qtlobj<-makeqtl(data,qtlist[[i]][,1],qtlist[[i]][,2],what="prob")
    rqtl[[p]]<-refineqtl(data,p,qtlobj,method="hk")
    out.aq[[p]]<-addqtl(data,qtlist[[i]][,1],p,rqtl[[p]],method="hk",verbose=T)
  }
}
#binary #add QTL
qtlist.bin<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=0.63,ci.function="bayesint")
data<-calc.genoprob(data,3,map.function="kosambi")
out.aq.bin<-list();rqtl.bin<-list()
for(i in 1:2){
  if(length(qtlist.bin[[i]][,1])>0){
    p<-phenames(data)[i+n]
    qtlobj.bin<-makeqtl(data,qtlist.bin[[i]][,1],qtlist.bin[[i]][,2],what="prob")
    rqtl.bin[[p]]<-refineqtl(data,p,qtlobj.bin,method="hk",model="binary")
    out.aq.bin[[p]]<-addqtl(data,qtlist.bin[[i]][,1],p,rqtl.bin[[p]],maxit=1e+9,tol=0.05,method="hk",model="binary",verbose=T)
    #maxit controls the trade-off between computational resources and the precision of the optimization algorithm. Increasing maxit may improve accuracy but also increases computation time.
  }
}
```

Threshold colors

```{r}
thcol<-c('blue','green','red')
```

Plot the LOD peaks, from the original scan and from the additional scan

```{r}
#plot qtl peaks of scanone
plot.sc<-function (x,bin=F,sc=scan1,thresh=thresh1.hk,LETTERs=T,l=1,mf=c(2,2),second=F,first=F){
  if (!second) par(mfrow=mf)
  for (i in x) {
    p<-phenames(data)[i];if(bin)p<-phenames(data)[i+n]
    plot(sc,lodcolumn=i,main=p,ylab="LOD",bandcol="gray80",ylim=c(0,max(sc[,2+i])+0.5),alternate.chrid = T)
    abline(h=thresh[,i],lty='dotted',lwd=2,col=thcol)
    for(j in 1:3){
      if(thresh[j,i]/par('usr')[4]<1){
        mtext(rownames(thresh)[j],side=4,font=2,adj=thresh[j,i]/(par('usr')[4]-0.2),col=thcol[j])
      }
    }
    if (LETTERs) mtext(LETTERS[l],adj=0,cex=1.2);l<-l+1;if(l>6)l<-1
  }
  if (first==F) par(mfrow=c(1,1))
}
```

Plot LOD score by phenotype

```{r}
plot.sc(1:6)
plot.sc(7:10,first=T)
plot.sc(1,T,sc=scan1.bin,thresh=thresh2.em,second=T,l=5)
plot.sc(2,T,sc=scan1.bin,thresh=thresh2.em,LETTERs=F,first=T)
plot.sc(11,LETTERs=F,second=T)
```

Customized QTL plotting function, combining original scan and addqtl scan

```{r}
plotAddqtl<-function(x,bin=F,list=qtlist,aq=out.aq,thresh=thresh1.hk,mfrow=c(2,2),second=F,LETTERs=T,l=1){
  if(!second)par(mfrow=mfrow)
  par(cex.lab=1.5,cex.axis=1.3,cex.main=1.7,cex.sub=1.3)
  for(i in x){
    if(length(list[[i]][,1])>0){
      p<-phenames(data)[i];if(bin)p<-phenames(data)[i+n]
      plot(aq[[p]],alternate.chrid=nrow(list[[p]])>3,ylab="LOD")
      abline(h=thresh[,i],lty='dotted',lwd=2,col=thcol)
      for(j in 1:3){
        if(thresh[j,i]/par('usr')[4]<1){
          mtext(rownames(thresh)[j],cex=0.9,font=2,adj=thresh[j,i]/(par('usr')[4]-0.1),side=4,col=thcol[j])}
      }
      title(p)
      if(nrow(list[[p]])==1)title(sub=list[[p]][,1])
      if (LETTERs) mtext(LETTERS[l],adj=0,cex=1.3);l<-l+1;if(l>6)l<-1
    }
  }
}
```

Plot the QTL, LOD graphs

```{r}
plotAddqtl(1:6,mfrow=c(2,3))
plotAddqtl(7:10)
plotAddqtl(1,T,qtlist.bin,out.aq.bin,thresh2.em,second=T,l=4)
plotAddqtl(1,T,qtlist.bin,out.aq.bin,thresh2.em,l=6)
plotAddqtl(2,T,qtlist.bin,out.aq.bin,thresh2.em,LETTERs=F)
```

Classify the phenotypes for the presence of joint interaction of markers.

```{r}
intpPhen<-vector();effpPhen<-vector()
for (i in 1:n) {
  (sc1<-summary(scan1,perms=scan1perm,alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
  p<-phenames(data)[i]
  cat(nrow(sc1),"QTL in",p,"\n")
  if(nrow(sc1)>1) intpPhen<-c(intpPhen,p)
  if(nrow(sc1)==1) effpPhen<-c(effpPhen,p)
  }
  #binary
  intpPhen.bin<-vector();effpPhen.bin<-vector()
  for (i in 1:2) {
    (sc1<-summary(scan1.bin,perms=scan1perm.bin,alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
    p<-phenames(data)[i+n]
    cat(nrow(sc1),"QTL in",p,"\n")
    if(nrow(sc1)>1) intpPhen.bin<-c(intpPhen.bin,p)
    if(nrow(sc1)==1) effpPhen.bin<-c(effpPhen.bin,p)
  }
```

```{r}
options(warn=0)
#normal
qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",alpha=0.95,ci.function="bayesint",pvalues=T)
  data<-calc.genoprob(data,2,map.function="kosambi")
  qtlist.aq<-list()
  s.aq<-list()
  for (i in 1:n){
    p<-phenames(data)[i]
    if(!is.null(out.aq[[p]])){
      s<-summary(out.aq[[p]],format="tabByCol",perms=scan1perm[,p],alpha=0.95,ci.function="bayesint",pvalues=T)
      if(nrow(s[[1]])>0){
        qtlist.aq[p]<-s
        qtlist.aq[[p]]<-cbind.data.frame(Trait=p,qtlist.aq[[p]])
        s.aq[[p]]<-summary(out.aq[[p]],perms=scan1perm[,p],alpha=0.63)
        if(nrow(s.aq[[p]])>0){
          rqtl[[p]]<-addtoqtl(data,rqtl[[p]],s.aq[[p]][,1],s.aq[[p]][,2])
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
      qtlist[[i]]<-rbind(qtlist[[i]],qtlist.aq[[phenames(data)[i]]])
    }
  }
  qtldf.aq<-do.call(rbind.data.frame,c(qtlist,make.row.names=F))
  
  #binary
  qtlist<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=0.95,ci.function="bayesint",pvalues=T)
  data<-calc.genoprob(data,3,map.function="kosambi")
  for (i in 1:2){
    p<-phenames(data)[i+n]
    if(!is.null(out.aq.bin[[p]])){
      s<-summary(out.aq.bin[[p]],perms=scan1perm.bin[,p],alpha=0.95,format="tabByCol",ci.function="bayesint",pvalues=T)
      if(nrow(s[[1]])>0){
        qtlist.aq[p]<-s
        qtlist.aq[[p]]<-cbind.data.frame(Trait=p,qtlist.aq[[p]])
        s.aq[[p]]<-summary(out.aq.bin[[p]],perms=scan1perm.bin[,p],alpha=0.63)
        if(nrow(s.aq[[p]])>0){
          rqtl.bin[[p]]<-addtoqtl(data,rqtl.bin[[p]],s.aq[[p]][,1],s.aq[[p]][,2])
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
      qtlist[[i]]<-rbind(qtlist[[i]],qtlist.aq[[phenames(data)[i+n]]])
    }
  }
  qtldf.aq.bin<-do.call(rbind.data.frame,c(qtlist,make.row.names=F))
  qtldf.aq<-rbind(qtldf.aq,qtldf.aq.bin)

```

QTL summary as dataframe for final report. QTL at alp= 0.99 and sig \*\*\* levels

```{r}
su<-1-summary(data)$missing.phe
qtldf.aq<-qtldf.aq%>%
  mutate("Len of LG"=round(chrlen(data)[chr],1),.after=chr)%>%
  mutate("Len of QTL"=round(ci.high-ci.low,1),.after="Len of LG")%>%
  mutate("Flanking markers"=paste0(chr,"_m",find.marker(data,chr,ci.low),"-",chr,"_m",find.marker(data,chr,ci.high)))%>%
  mutate("Central marker"=paste0(chr,"_m",find.marker(data,chr,pos)))%>%
  mutate("Pval"=paste0(pval,if_else(pval<0.63,"*",""),if_else(pval<0.1,"*",""), if_else(pval<0.05,"*","")))%>%
  select(!c(pval,ci.low,ci.high))%>%
  rename("QTL's LG"=chr)%>%
  mutate("No. Inds/% phenotyped"=paste0(nind(data)*su[find.pheno(data,Trait)]," ind / ",round(100*su[find.pheno(data,Trait)],1),"%"),.after=Trait)%>%
  mutate("pos"=round(pos,1))%>%
  mutate("lod"=round(lod,1))
qtldf.aq
```

Customized merged interaction plots

```{r}
mergedIntp<-function(x=1,bin=F,sc=scan1,perm=scan1perm,rit.inx=NULL,mf=c(2,3),l=1,LETTERs=T,second=F,first=F){
  rit<-1 #legend on the right
  if (!second) par(mfrow=mf) #if plot is second then letter is continuous
  for (i in x) {
    (sc1<-summary(sc,perm=perm,alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
    (r<-nrow(sc1))
    sorted<-sc1%>%arrange(desc(across(3)))
    (mn<-find.marker(data,chr=sorted[,1],pos=sorted[,2]))
    for(j in 1:(r-1)){
      for(k in (j+1):r){
        if (bin) i<-i+n
        effectplot(data,pheno.col=i,mname1=mn[j],mname2=mn[k],main="",ylab=paste("Ave. phenotype:",phenames(data)[i]),xlab=paste0(sorted[k,1],"_m",mn[k]),add.legend = F)
        if(rit%in%rit.inx) lpos<-"topright" else lpos<-"topleft"
        legend(lpos,c("AA","AB","BB"),lty=1,pch=1,col=c("black","red","blue"),bty="n",inset=c(0.01,0))
        a<-par("usr")
        x.leg <- a[1] * 0.05 + a[2] * 0.75
        y.leg <- a[4] - diff(a[3:4]) * 0.05
        if(rit%in%rit.inx) tpos<-2 else tpos<-NULL
        text(x.leg,y.leg,pos=tpos,paste0(sorted[j,1],"_m",mn[j]))
        if (LETTERs)mtext(LETTERS[l],adj=0,cex=1.1);l<-l+1;if(l>6)l<-1
        rit<-rit+1}
    }
  }
  if (!first) par(mfrow=c(1,1))
}
```

Plot marker interactions

```{r}
data<-sim.geno(data,step=3)
x<-(1:6)[1:6%in%find.pheno(data,intpPhen)]
mergedIntp(x,rit.inx=c(3,6,8,12,15,21,22,26,35))
```

```{r eval=FALSE, include=FALSE}
x<-(7:11)[7:11%in%find.pheno(data,intpPhen)]
mergedIntp(x,first=T,rit.inx=1,mf=c(2,2))
#binary
x<-((12:13)[12:13%in%find.pheno(data,intpPhen.bin)])-n
mergedIntp(x[1],T,scan1.bin,scan1perm.bin,second=T,l=2)
mergedIntp(x[2],T,scan1.bin,scan1perm.bin,LETTERs=F)
```

QTL pairs summary

```{r}
#normal
c.thr1<-list()
for(i in 1:n){
  (thr1<-summary(scan2, perms=scan2perm, alpha=0.2,lodcolumn=i,pvalues=T))
  if(i==1){c.thr1[[phenames(data)[i]]]<-thr1
  }else c.thr1[[phenames(data)[i]]]<-thr1
}
for (i in 1:length(c.thr1)){
  if(nrow(c.thr1[[i]])>0&&colnames(c.thr1[[i]])[1]!="Trait"){
    c.thr1[[i]]<-cbind.data.frame(Trait=names(c.thr1[i]),c.thr1[[i]])
  }
}
thr1df<-do.call(rbind.data.frame,c(c.thr1,make.row.names=F))
#binary
c.thr2<-list()
for(i in 1:2){
  (thr2<-summary(scan2.bin, perms=scan2perm.bin, alpha=0.2,lodcolumn=i,pvalues=T))
  if(i==1){c.thr2[[phenames(data)[i+n]]]<-thr2
  }else c.thr2[[phenames(data)[i+n]]]<-thr2
}
for (i in 1:length(c.thr2)){
  if(nrow(c.thr2[[i]])>0&&colnames(c.thr2[[i]])[1]!="Trait"){
    c.thr2[[i]]<-cbind.data.frame(Trait=names(c.thr2[i]),c.thr2[[i]])
  }
}
thr2df<-do.call(rbind.data.frame,c(c.thr2,make.row.names=F))
thr<-rbind(thr1df,thr2df)
```

QTL pairs

```{r}
#normal
data<-calc.genoprob(data,2,map.function="kosambi")
qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",alpha=0.63,ci.function="bayesint")
sc2thr1<-summary(scan2perm,alpha=0.2)
#rearrange the threshold list
th<-vector('list',5)
for(j in 1:5){th[[j]]<-t(sc2thr1[[j]])}
m<-do.call('cbind',th)
dimnames(m)<-list(phenames(data)[1:n],names(sc2thr1)[1:5])
out.ap<-list();qtlpairs<-list();s.fq<-list();out.fq<-list()
for(i in 1:n){
  p<-phenames(data)[i]
  if(length(qtlist[[p]][,1])>0){
    out.ap[[p]]<-addpair(data,qtlist[[p]][,1],p,rqtl[[p]],method="hk",verbose=T)
    qtlpairs[[p]]<-summary(out.ap[[p]],thresholds=m[p,])
    s.fq[[p]]<-summary(out.fq[[p]]<- fitqtl(data,p,rqtl[[p]],method="hk",get.ests=T))
  }
}
#binary
data<-calc.genoprob(data,10,map.function="kosambi")
qtlist.bin<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=0.63,ci.function="bayesint")
sc2thr2<-summary(scan2perm.bin,alpha=0.2)
#rearrange the threshold list
th<-vector('list',5)
for(j in 1:5){th[[j]]<-t(sc2thr2[[j]])}
m.bin<-do.call('cbind',th)
dimnames(m.bin)<-list(phenames(data)[n+1:2],names(sc2thr2)[1:5])
out.ap.bin<-list();rqtl2.bin<-list()
for(i in 1:2){
  p<-phenames(data)[i+n]
  if(length(qtlist.bin[[p]][,1])>0){
    q<-rbind(qtlist.bin[[p]][,-c(3,4)],s.aq[[p]])
    rqtl2.bin[[p]]<-refineqtl(data,p,makeqtl(data,q[,1],q[,2],what="prob"),maxit.fitqtl=1e+6,tol=0.05,method="hk",model="binary")
    out.ap.bin[[p]]<-addpair(data,q[,1],p,rqtl2.bin[[p]],maxit=1e+6,tol=0.2,method="hk",model="binary",verbose=T)
    qtlpairs[[p]]<-summary(out.ap.bin[[p]],thresholds=m.bin[p,])
    s.fq[[p]]<-summary(out.fq[[p]]<- fitqtl(data,p,rqtl2.bin[[p]],maxit=1e+6,tol=0.01,method="hk",model="binary",get.ests=T))
  }
}
for (i in 1:length(qtlpairs)){
if(nrow(qtlpairs[[i]])>0 && names(qtlpairs[[i]])[1]!="Trait"){
  qtlpairs[[i]]<-cbind.data.frame(Trait=names(qtlpairs[i]),qtlpairs[[i]])
  }
}
qtlpairsdf<-do.call(rbind.data.frame,c(qtlpairs,make.row.names=F))
#interacting QTL
qtlpairsdf<-c(qtlpairsdf,thr=m[qtlpairsdf[,1],])
```

No interactive QTL pair was found in the add pair scan.

# Linkage Map view

```{r}
alp<-0.63
colorlist<-RColorBrewer::brewer.pal(8,"Set1")
```

Genetic map to pdf

```{r}
qtldf_initial<-\(){ 
  # make a df to pass qtl info
  qtldf <- data.frame(
    chr = character(),
    qtl = character(),
    so = numeric(),
    si = numeric(),
    ei = numeric(),
    eo = numeric(),
    col = character(),
    stringsAsFactors = F
  )
  return(qtldf)
}
outfile<-file.path("results/basil_linkage_map.pdf")
main<-"Basil Genetic Map"
qtldf<-qtldf_initial()

setting<-list(mapthis=data,outfile=outfile,main=main,ruler=T,maxnbrcolsfordups=2,dupnbr=T,lg.col='lightblue1',lgw=0.15,labdist=0.15,lgperrow=3)
do.call(lmv.linkage.plot,setting)
```

Anthocyanin QTL map

```{r}
qtldf<-qtldf_initial()
for (i in 1:6) {
  (qtls<-summary(scan1,perms=scan1perm,alpha=alp,lodcolumn=i)[,c(1:2,2+i)])
  if(nrow(qtls)>0){
    for (j in 1:nrow(qtls)) {
      (bay <-bayesint(scan1[,c(1:2,2+i)],chr=qtls$chr[j]))
      qtldf <- rbind(qtldf,
                     data.frame(
                       chr = qtls$chr[j],
                       qtl = colnames(bay)[3],
                       so = bay$pos[1],
                       si = bay$pos[2],
                       ei = bay$pos[2],
                       eo = bay$pos[3],
                       col=colorlist[(i+1)]))
    }
  }
}
outfile<-file.path("results/basil_QTLs.anthocyanin.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Anthocyanin (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main,qtldf=qtldf))
do.call(lmv.linkage.plot,setting)
```

Fusarium QTL map

```{r}
qtldf<-qtldf_initial()
for (i in 7:10) {
  qtls<-summary(scan1,perms=scan1perm,alpha=alp,lodcolumn=i)[,c(1:2,2+i)]
  if(nrow(qtls)>0){
    for (j in 1:nrow(qtls)) {
      (bay <-bayesint(scan1[,c(1:2,2+i)],chr=qtls$chr[j]))
      qtldf <- rbind(qtldf,
                     data.frame(
                       chr = qtls$chr[j],
                       qtl = colnames(bay)[3],
                       so = bay$pos[1],
                       si = bay$pos[2],
                       ei = bay$pos[2],
                       eo = bay$pos[3],
                       col=colorlist[(i-5)]))
    }
  }
}
i<-1
(qtls<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=alp,ci.function="bayesint"))
(p<-phenames(data)[i+n])
for (j in 1:(nrow(qtls[[p]])+1)) {
  if(j==3){
    qtls[[p]]<-rbind(qtls[[p]],cbind(s.aq[[p]][,-3],`ci.low`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[1,2],`ci.high`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[3,2],`lod`=s.aq[[p]][,3]))
  }
  qtldf<-rbind(qtldf,
               data.frame(
                 chr = qtls[[p]]$chr[j],
                 qtl = p,
                 so = qtls[[p]]$`ci.low`[j],
                 si = qtls[[p]]$pos[j],
                 ei = qtls[[p]]$pos[j],
                 eo = qtls[[p]]$`ci.high`[j],
                 col=colorlist[7]))
}
outfile<-file.path("results/basil_QTLs.fusarium.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Fusarium (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main))
setting$qtldf<-qtldf
do.call(lmv.linkage.plot,setting)
```

Downy Mildew QTL map

```{r}
qtldf<-qtldf_initial()
i<-2
(qtls<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=alp,ci.function="bayesint"))
p<-phenames(data)[i+n]
for (j in 1:(nrow(qtls[[p]])+1)) {
  if(j==3)qtls[[p]]<-rbind(qtls[[p]],cbind(s.aq[[p]][,-3],`ci.low`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[1,2],`ci.high`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[3,2],`lod`=s.aq[[p]][,3]))
  qtldf<-rbind(qtldf,
               data.frame(
                 chr = qtls[[p]]$chr[j],
                 qtl = p,
                 so = qtls[[p]]$`ci.low`[j],
                 si = qtls[[p]]$pos[j],
                 ei = qtls[[p]]$pos[j],
                 eo = qtls[[p]]$`ci.high`[j],
                 col=colorlist[8]))
}
outfile<-file.path("results/basil_QTLs.BDM.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Downy Mildew (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main))
setting$qtldf<-qtldf
do.call(lmv.linkage.plot,setting)
```

The pdf files are in the results folder.

There is more for that project that I have investigated, such as heatmap of scantwo output, multiple qtl model analysis, larger genome-based dataset (5k markers) with the same workflow but takes more computational resources, alignment of the genetic and genomic (large) map side by side anchored by shared markers, comparing the qtl from both maps. I didn't add it to this notebook because u can see it's quiet much already.

---
title: "Genetic Mapping and QTL Analysis"
output: github_document
---

<img src="https://octodex.github.com/images/labtocat.png" width="200" height="200"/>

Load packages

```{r}
library(qtl)
library(qtlcharts)
library(ASMap)
library(dplyr)
library(corrplot)
library(LinkageMapView)
```

Import genotypic data. for now, all markers are on 1 chr

```{r}
data<-read.cross("csvr",".","raw_data.csv",genotypes=c("a","h","b"),map.function="kosambi")
```

Convert cross to bcsft type f2 to use ASMap package

```{r}
data<-convert2bcsft(data,F.gen=2,estimate.map=F)
cat("Phenos:",(head(phenames(data),-1)))
n<-11 #normal phenotypes count
```

# Data Pre-processing

Look at the pattern of missing data. Black pixels indicate missing genotypes.

```{r}
plotMissing(data)
```

Omit the individuals with more than 85% of total markers -\> less than 672 markers. Briefly estimate map first.

```{r}
data<-quickEst(data,map.function="kosambi")
sg<-statGen(data,bychr=F,stat.type="miss",id='index')
data1<-subset(data,ind=om<-sg$miss<(totmar(data))*0.85)
cat(if(F%in%om) which(!om)else "no","ind omitted for missing > 85% mar")
```

Plot the number of genotyped markers per individual as well as the number of genotyped individuals per marker

```{r}
par(mfrow=c(1,2), las=1)
plot(ntyped(data1), ylim=c(0,totmar(data1)+50),ylab="No. typed markers",main="Markers by Individual")
mtext("A",adj=0)
plot(ntyped(data1, "mar"), ylim=c(0,(nind(data1)+15)),ylab="No. typed individuals",main="Individuals by Marker")
mtext("B",adj=0)
```

Plot the genotype frequencies per individual

```{r}
g <- pull.geno(data1)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3){
  plot(gfreq[i,], ylab="Genotype frequency",ylim=c(0,1))
  abline(h=mean(gfreq[i,]),lty=3,col="red",lwd=3)
  mtext(c("AA", "AB", "BB")[i])}
par(mfrow=c(1,1));title("Genotypes' Frequency and Segregation Ratio",line = 2.5)
```

Compare the genotypes for all pairs of individuals

```{r}
cg<-comparegeno(data1);cgr<-cg[lower.tri((cg))]
hist(cgr, breaks=seq(0, 1, len=101),xlab="No. matching genotypes", main="Matching Pairs of Individuals")
rug(cgr)
#mark the outlier with a red arrow 
x<-max(cgr)
arrow.plot(x,50,0,-1,true.angle = T,arrow.ex=30, length=.1,col='red', lwd=2)
text(x,70,paste0(round(100*x,1),"%"),adj=c(0.5,0.2))
```

Omit individuals with more than 90% identical markers

```{r}
wh<-which(cg>0.9,arr=T)
data2<-subset(data1,ind=-wh[,2])
cat(nind(data1)-nind(data2),"ind omitted for",paste0(round(100*x,1),"%"),"identical geno\n",paste0("#", data1$pheno$index[wh[,2]]))
```

Pull out markers from cross temporarily

```{r}
cat(totmar(data2),'total mar\n')
data3<-pullCross(data2,type="missing",pars=list(miss.thresh=0.1))
cat(totmar(data3),'total mar\n')
cat(totmar(data2)-totmar(data3),"mar pulled for missing\n")

data4<-pullCross(data3,type="seg.distortion",pars=list(seg.thresh=0.001))
cat(totmar(data4),'total mar\n')
cat(totmar(data3)-totmar(data4),"mar pulled for seg. distortion\n")

data5<-pullCross(data4,type="co.located")
cat(totmar(data5),'total mar\n')
cat(totmar(data4)-totmar(data5),"mar pulled for co. located")
```

Plot optional p. values to determine distance threshold for marker clustering. I chose pValue= 1e-7 (on y axis), meaning-\>\> split linkage groups with more than 30cM gap between markers, according to 150 ind population (on x axis)

```{r}
cat(nind(data5),"individuals are used to cluster markers")
pValue(dist=seq(25,40,by=5),pop.size=110:190)
```

LOD(Logarithm of the Odds)= statistical measure of the likelihood that two loci (positions on a chromosome) are linked and therefore inherited together, rather than assorting independently.

# Map Construction

Form linkage groups with LOD=7, raw map

```{r}
data5<-mstmap(data5,bychr=F,p.value=1e-7,id='index')
plotMap(data5,alternate.chrid = T)
```

Profile individuals' genotype statistics

```{r include=FALSE}
pg<-profileGen(data5,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=20,layout=c(1,3),lty=2,cex=0.7)
pg<-profileGen(data5,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg$stat$xo),layout=c(1,3),lty=2,cex=0.8)
```

Omit missing/ double crossover (dxo)/ xo statistics outlier individuals

```{r}
data6<-subsetCross(data5,ind=!pg$xo.lambda)
cat(nind(data5)-nind(data6),"ind omitted by profileGen")
```

Double check dxo. An unusually high rate of double crossovers might indicate genotyping errors

```{r}
pg1<-profileGen(data6,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg$stat$xo),layout=c(1,3),lty=2,cex=0.7)
pg1<-profileGen(data6,bychr=F,stat.type=c("xo","dxo","miss"),id="index",xo.lambda=median(pg1$stat$xo),layout=c(1,3),lty=2,cex=0.7)
```

Re-construct map. The genotyping errors can distort the distances between markers, the order of the inputted markers is respected

```{r}
data7<-mstmap(data6,bychr=F,anchor=T,p.value=1e-7,id='index')
```

Push back markers to the map

```{r}
cat(totmar(data7),'total mar\n ')
data8<-pushCross(data7,type="co.located")
cat(totmar(data8),'total mar\n ')
cat(totmar(data8)-totmar(data7),"mar pushed for co. located")
```

Re-construct final map by adding markers to existing LGs (linkage groups)

```{r include=FALSE}
data9<-mstmap(data8,anchor=T,p.value=2,id='index')
```

Drop LGs with less than 2 markers

```{r}
mndrop<-markernames(data9,nmar(data9)<2)
data10<-drop.markers(data9,as.character(mndrop))
cat(totmar(data10),'total mar\n')
cat(totmar(data9)-totmar(data10),"mar omitted for LG < 2 mar")
```

Rename chr by numerical order

```{r}
x<-1:nchr(data10)
for (i in x) {names(data10$geno)[i]<-paste0("LG",i)}
```

Plot genetic map illustration, final map

```{r}
plot.map(data10,alternate.chrid=T)
```

LGs summary table

```{r}
summaryMap(data10)
```

# Map Evaluation

Estimate recombination fraction

```{r}
data10<-est.rf(data10)
```

Heatmap of LOD and Rf

```{r}
heatMap(data10,lmax = 50,main='')
mtext("Pairwise Recombination Fractions and LOD Scores",cex=1.1,line=3.2,adj=0.4,font=2)
```

Plot Pairwise LOD vs. Rf

```{r}
rf<-pull.rf(data10);lod<-pull.rf(data10,what="lod")
plot(as.numeric(rf),as.numeric(lod),xlab="Recombination fraction",ylab="LOD score",main=paste("Pairwise LOD vs. Rf for",totmar(data10),"Markers"))
```

The evaluation looks good:

The heatmap is continuous with low gradient along the chr meaning that the markers within a chromosome are gradually distant from each other.

On the LOD/rf scatter plot there is a trend line with no outliers.

# Exploratory Data Analysis

Jitter map- to avoid marker overlaping by slightly adding gaps between them

```{r}
data<-jittermap(data10)
```

Create a df of parents phenotype

```{r}
cold_values <- c(5, 5, 4, 6, 7, 9, 5, 6, 6, 7, NA, NA, 7, 8, 7, 8, 8, 8, 9, 9, 9, 8, NA, NA)
faudpc_values <- c(rep(0, 12), rep(c(90, 142), each = 6))
parents <- data.frame(rbind(matrix(0,12,9), matrix(5,12,9)),faudpc_values, Cold = cold_values,row.names = c(paste0("P1_", 1:12), paste0("DP_", 1:12)))
ant_values <-c(5, 4.5, 4.5, 4.5, 4.5, 3.5, 4, 4.5, 4, 4.5, 4, 4)
parents[13:24,1:3] <-matrix(ant_values,12,3)
colnames(parents) <-phenames(data)[1:11]
head(parents)
```

Check the phenotypic distribution

Plot histograms or barplots

```{r}
cbind(phenames(data)) 
```

```{r}
par(mfrow = c(2, 3))
l <- 1
for (i in c(1:10, 12, 15, 13, 11)) {
  if (i < 12) {
    plotPheno(data, i, ylab = "Frequency", xlab = if (i < 7) "Purple intensity" else "Disease intensity")}
  if (i == 15) {
    par(mfrow = c(2, 3))}
  if (i > 11 & i != 15) {
    BP <- plotPheno(data, i, ylab = "Frequency",xlab = "Disease resistance",names.arg = c("R", "S"))
    text(BP, table(data$pheno[, i]),labels = table(data$pheno[, i]),pos = 1)
    points(c(0.7, 1.9), rep(1.6, 2), pch = 25,cex = 1.3,bg = c("palegreen1", "violet"))}
   if (i < 10) {
     points(if (i < 4) 0.1 else 0.7, 1.6, pch = 25, cex = 1.3, bg = "palegreen1")}
   if (i < 4) {
     points(mean(parents[13:24, i]) - 0.4, 1.6, pch = 25, cex = 1.3, bg = "violet")}
   if (i == 10) {
    points(c(mean(parents[1:12, i]) - 0.1, mean(parents[13:24, i])), rep(1.6, 2), pch = 25, cex = 1.3, bg = c("palegreen1", "violet"))}
  if (i == 11) {
    points(c(mean(parents[1:12, i], na.rm = TRUE), mean(parents[13:24, i], na.rm = TRUE)), rep(1.6, 2), pch = 25, cex = 1.3, bg = c("palegreen1", "violet"))}
  if (i %in% c(4, 5, 6, 9)) {
    points(if (i == 4) 6.7 else if (i == 9) 6.7 else if (i == 5 || i == 6) 4.3, 1.6, pch = 25, cex = 1.3, bg = "violet")}
   if (i %in% c(7, 8)) {
     points(5.5, 1.6, pch = 25, cex = 1.3, bg = "violet")}
   if (l > 6) l <- 1
   if (i %in% c(1:10, 12)) {
     if (i == 7) l <- 1
     mtext(LETTERS[l], adj = -0.1, cex = 1.1)}
   l <- l + 1
}
```

Plot the correlation matrices using corrplot

```{r}
mat<-cor(pull.pheno(data,1:13),use = "complete.obs")
corrplot(mat, method = "color", addCoef.col = "orange", tl.col = "black", tl.srt = 35,tl.cex=0.7,number.cex=0.5) 
```

Now, for each group of phenotypes

Fusarium

```{r}
mat<-cor(pull.pheno(data,c(7:10,12)),use = "complete.obs")
corrplot(mat, method = "color", addCoef.col = "white", tl.col = "black", tl.srt = 35)
```

Anthocyanin

```{r}
mat<-cor(pull.pheno(data,1:6),use = "complete.obs")
corrplot(mat, method = "color", addCoef.col = "white", tl.col = "black", tl.srt = 35)
```

There is a strong correlation between the phenotypes in the same group. the different conditions (in anthocyanin) or repetitions (in fusarium) had a limited effect on the resulting phenotype.

# QTL Analysis

Perform genome scans to identify QTL

```{r}
if (all(file.exists("genome_scans.Rdata","scan2_part1.Rdata","scan2_lod_part2.Rdata"))) {
  message("Loading precomputed genome scans from genome_scans.Rdata to save time.")
  load("genome_scans.Rdata")
  # Load the parts (scan2 is a large file so i split it to 2 parts)
  load("scan2_part1.Rdata")
  load("scan2_lod_part2.Rdata")
  lod_part1<-scan2.partial$lod
  # Combine the 'lod' matrix parts back into a single matrix
  combined_lod <- array(NA, dim = c(dim(lod_part1)[1] + dim(lod_part2)[1], dim(lod_part1)[2], dim(lod_part1)[3]))
  # Copy the data from the parts into the combined array
  combined_lod[1:dim(lod_part1)[1], , ] <- lod_part1
  combined_lod[(dim(lod_part1)[1] + 1):dim(combined_lod)[1], , ] <- lod_part2
  # Set the dimnames attribute correctly
  dimnames(combined_lod) <- dimnames(scan2.partial$lod)
  # Reconstruct the scantwo object
  scan2 <- scan2.partial  # Copy the original object structure
  scan2$lod <- combined_lod  # Replace the 'lod' part with the combined matrix
  # cleanup
  rm(list=c("scan2.partial","combined_lod","lod_part1","lod_part2"))
} else {
  cat("The precomputed genome scans were not found\n")
  message("Performing genome scans... This may take a long time.")
  
  ###ScanOne###
  #calculate probabilities (necessary for scanone)
  #normal #Haley-knott regression 
  data<-calc.genoprob(data,2,map.function="kosambi")
  #genome scan for Single-QTL
  scan1<-scanone(data,pheno.col=c(1:n),method="hk")
  #permutation test
  operm <- vector("list", 100)
  for(i in 1:100){operm[[i]]<-scanone(data,pheno.col=c(1:n),method="hk",n.perm=9,n.cluster=2)}
  scan1perm<-do.call("rbind", operm)
  
  #binary #EM algorithm maximum likelihood
  data<-calc.genoprob(data,3,map.function="kosambi")
  #genome scan for Single-QTL
  scan1.bin<-scanone(data,pheno.col=c(12:13),method="em",model="binary")
  #permutation test
  scan1perm.bin1<-scanone(data,pheno.col=c(12:13),method="em",model="binary",n.perm=1000,n.cluster=2)
  operm <- vector("list", 100)
  for(i in 1:100){operm[[i]]<-scanone(data,pheno.col=c(12:13),method="em",model="binary",n.perm=10,n.cluster=2)}
   scan1perm.bin1<-do.call("rbind", operm)
   
   save(scan1,scan1.bin,scan1perm,scan1perm.bin,scan2.bin,scan2perm,scan2perm.bin,file="genome_scans.Rdata")
   
   ###ScanTwo###
   #normal #scan genome for Two-QTL model
   data<-clean(data)
   data<-calc.genoprob(data,2,map.function="kosambi")
   scan2<-scantwo(data,pheno.col=c(1:n),method="hk")
   #permutation test
   ##this proccess is heavy so I've split it into 4 parts and then merged the results. with better CPU performance U can run the following line instead:
   ##scan2perm.bin<-scantwo(data,pheno.col=c(1:n),method="hk",n.perm=1000)
   #res1
   operm2 <- vector("list", 100)
   for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=c(1:6),method="hk",n.perm=9,n.cluster=2)}
   res1<-do.call("rbind", operm2)
   #res2
   operm2 <-vector("list", 100)
   for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=c(7:8),method="hk",n.perm=9,n.cluster=2)}
   res2<-do.call("rbind", operm2)
   #res3
   operm2 <- vector("list", 100)
   for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=c(9:10),method="hk",n.perm=9,n.cluster=2)}
   res3<-do.call("rbind", operm2)
   #res4
   operm2 <- vector("list", 100)
   for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=11,method="hk",n.perm=9,n.cluster=2)}
   res4<-do.call("rbind", operm2)
   #cbind the 11 phenos 
   scan2perm<-res1
   for(i in 1:6){scan2perm[[i]]<-cbind(res1[[i]],res2[[i]],res3[[i]],res4[[i]])}
   
   #binary #genome scan for two-QTL
   data<-calc.genoprob(data,10,map.function="kosambi")
   scan2.bin<-scantwo(data,pheno.col=c(12,13),method="em",model="binary",verbose=T)
   #permutation test
   ##this proccess is heavy so I've split it into 4 parts and then merged the results. with better CPU performance U can run the following line instead:
  ##scan2perm.bin<-scantwo(data,pheno.col=c(12:13),method="em",model="binary",n.perm=1000)
  operm2 <- vector("list", 200)
  for(i in 1:100 ){operm2[[i]]<-scantwo(data,pheno.col=c(12:13),method="em",model="binary",n.perm=5,n.cluster=2)}
   for(i in 101:200 ){operm2[[i]]<-scantwo(data,pheno.col=c(12:13),method="em", model="binary",n.perm=5,n.cluster=2)}
   scan2perm.bin<-do.call("rbind", operm2[1:200])
}
```

Setting a QTL detection threshold according to permutation tests. The lower the precentage (5%), the better significance of the QTL.

```{r}
(thresh1.hk<-summary(scan1perm,alpha=c(0.63,0.1,0.05))) 
(thresh2.em<-summary(scan1perm.bin,alpha=c(0.63,0.1,0.05)))
```

Scan for additional QTL after reducing the masking effect of QTL with major peaks.

```{r}
#normal #add additional QTL
qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",alpha=0.63,ci.function="bayesint")
message("Calculating genoprob step= 2 cM")
data<-calc.genoprob(data,2,map.function="kosambi")
out.aq<-list();rqtl<-list()
for(i in 1:n){
  if(length(qtlist[[i]][,1])>0){
    p<-phenames(data)[i]
    qtlobj<-makeqtl(data,qtlist[[i]][,1],qtlist[[i]][,2],what="prob")
    rqtl[[p]]<-refineqtl(data,p,qtlobj,method="hk")
    out.aq[[p]]<-addqtl(data,qtlist[[i]][,1],p,rqtl[[p]],method="hk",verbose=T)
  }
}
#binary #add QTL
qtlist.bin<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=0.63,ci.function="bayesint")
data<-calc.genoprob(data,3,map.function="kosambi")
out.aq.bin<-list();rqtl.bin<-list()
for(i in 1:2){
  if(length(qtlist.bin[[i]][,1])>0){
    p<-phenames(data)[i+n]
    qtlobj.bin<-makeqtl(data,qtlist.bin[[i]][,1],qtlist.bin[[i]][,2],what="prob")
    rqtl.bin[[p]]<-refineqtl(data,p,qtlobj.bin,method="hk",model="binary")
    out.aq.bin[[p]]<-addqtl(data,qtlist.bin[[i]][,1],p,rqtl.bin[[p]],maxit=1e+9,tol=0.05,method="hk",model="binary",verbose=T)
    #maxit controls the trade-off between computational resources and the precision of the optimization algorithm. Increasing maxit may improve accuracy but also increases computation time.
  }
}
```

Threshold colors

```{r}
thcol<-c('blue','green','red')
```

Plot the LOD peaks, from the original scan and from the additional scan

```{r}
#plot qtl peaks of scanone
plot.sc<-function (x,bin=F,sc=scan1,thresh=thresh1.hk,LETTERs=T,l=1,mf=c(2,2),second=F,first=F){
  if (!second) par(mfrow=mf)
  for (i in x) {
    p<-phenames(data)[i];if(bin)p<-phenames(data)[i+n]
    plot(sc,lodcolumn=i,main=p,ylab="LOD",bandcol="gray80",ylim=c(0,max(sc[,2+i])+0.5),alternate.chrid = T)
    abline(h=thresh[,i],lty='dotted',lwd=2,col=thcol)
    for(j in 1:3){
      if(thresh[j,i]/par('usr')[4]<1){
        mtext(rownames(thresh)[j],side=4,font=2,adj=thresh[j,i]/(par('usr')[4]-0.2),col=thcol[j])
      }
    }
    if (LETTERs) mtext(LETTERS[l],adj=0,cex=1.2);l<-l+1;if(l>6)l<-1
  }
  if (first==F) par(mfrow=c(1,1))
}
```

Plot LOD score by phenotype

```{r}
plot.sc(1:6)
plot.sc(7:10,first=T)
plot.sc(1,T,sc=scan1.bin,thresh=thresh2.em,second=T,l=5)
plot.sc(2,T,sc=scan1.bin,thresh=thresh2.em,LETTERs=F,first=T)
plot.sc(11,LETTERs=F,second=T)
```

Customized QTL plotting function, combining original scan and addqtl scan

```{r}
plotAddqtl<-function(x,bin=F,list=qtlist,aq=out.aq,thresh=thresh1.hk,mfrow=c(2,2),second=F,LETTERs=T,l=1){
  if(!second)par(mfrow=mfrow)
  par(cex.lab=1.5,cex.axis=1.3,cex.main=1.7,cex.sub=1.3)
  for(i in x){
    if(length(list[[i]][,1])>0){
      p<-phenames(data)[i];if(bin)p<-phenames(data)[i+n]
      plot(aq[[p]],alternate.chrid=nrow(list[[p]])>3,ylab="LOD")
      abline(h=thresh[,i],lty='dotted',lwd=2,col=thcol)
      for(j in 1:3){
        if(thresh[j,i]/par('usr')[4]<1){
          mtext(rownames(thresh)[j],cex=0.9,font=2,adj=thresh[j,i]/(par('usr')[4]-0.1),side=4,col=thcol[j])}
      }
      title(p)
      if(nrow(list[[p]])==1)title(sub=list[[p]][,1])
      if (LETTERs) mtext(LETTERS[l],adj=0,cex=1.3);l<-l+1;if(l>6)l<-1
    }
  }
}
```

Plot the QTL, LOD graphs

```{r}
plotAddqtl(1:6,mfrow=c(2,3))
plotAddqtl(7:10)
plotAddqtl(1,T,qtlist.bin,out.aq.bin,thresh2.em,second=T,l=4)
plotAddqtl(1,T,qtlist.bin,out.aq.bin,thresh2.em,l=6)
plotAddqtl(2,T,qtlist.bin,out.aq.bin,thresh2.em,LETTERs=F)
```

Classify the phenotypes for the presence of joint interaction of markers.

```{r}
intpPhen<-vector();effpPhen<-vector()
for (i in 1:n) {
  (sc1<-summary(scan1,perms=scan1perm,alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
  p<-phenames(data)[i]
  cat(nrow(sc1),"QTL in",p,"\n")
  if(nrow(sc1)>1) intpPhen<-c(intpPhen,p)
  if(nrow(sc1)==1) effpPhen<-c(effpPhen,p)
  }
  #binary
  intpPhen.bin<-vector();effpPhen.bin<-vector()
  for (i in 1:2) {
    (sc1<-summary(scan1.bin,perms=scan1perm.bin,alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
    p<-phenames(data)[i+n]
    cat(nrow(sc1),"QTL in",p,"\n")
    if(nrow(sc1)>1) intpPhen.bin<-c(intpPhen.bin,p)
    if(nrow(sc1)==1) effpPhen.bin<-c(effpPhen.bin,p)
  }
```

```{r}
options(warn=0)
#normal
qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",alpha=0.95,ci.function="bayesint",pvalues=T)
  data<-calc.genoprob(data,2,map.function="kosambi")
  qtlist.aq<-list()
  s.aq<-list()
  for (i in 1:n){
    p<-phenames(data)[i]
    if(!is.null(out.aq[[p]])){
      s<-summary(out.aq[[p]],format="tabByCol",perms=scan1perm[,p],alpha=0.95,ci.function="bayesint",pvalues=T)
      if(nrow(s[[1]])>0){
        qtlist.aq[p]<-s
        qtlist.aq[[p]]<-cbind.data.frame(Trait=p,qtlist.aq[[p]])
        s.aq[[p]]<-summary(out.aq[[p]],perms=scan1perm[,p],alpha=0.63)
        if(nrow(s.aq[[p]])>0){
          rqtl[[p]]<-addtoqtl(data,rqtl[[p]],s.aq[[p]][,1],s.aq[[p]][,2])
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
      qtlist[[i]]<-rbind(qtlist[[i]],qtlist.aq[[phenames(data)[i]]])
    }
  }
  qtldf.aq<-do.call(rbind.data.frame,c(qtlist,make.row.names=F))
  
  #binary
  qtlist<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=0.95,ci.function="bayesint",pvalues=T)
  data<-calc.genoprob(data,3,map.function="kosambi")
  for (i in 1:2){
    p<-phenames(data)[i+n]
    if(!is.null(out.aq.bin[[p]])){
      s<-summary(out.aq.bin[[p]],perms=scan1perm.bin[,p],alpha=0.95,format="tabByCol",ci.function="bayesint",pvalues=T)
      if(nrow(s[[1]])>0){
        qtlist.aq[p]<-s
        qtlist.aq[[p]]<-cbind.data.frame(Trait=p,qtlist.aq[[p]])
        s.aq[[p]]<-summary(out.aq.bin[[p]],perms=scan1perm.bin[,p],alpha=0.63)
        if(nrow(s.aq[[p]])>0){
          rqtl.bin[[p]]<-addtoqtl(data,rqtl.bin[[p]],s.aq[[p]][,1],s.aq[[p]][,2])
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
      qtlist[[i]]<-rbind(qtlist[[i]],qtlist.aq[[phenames(data)[i+n]]])
    }
  }
  qtldf.aq.bin<-do.call(rbind.data.frame,c(qtlist,make.row.names=F))
  qtldf.aq<-rbind(qtldf.aq,qtldf.aq.bin)

```

QTL summary as dataframe for final report. QTL at alp= 0.99 and sig \*\*\* levels

```{r}
su<-1-summary(data)$missing.phe
qtldf.aq<-qtldf.aq%>%
  mutate("Len of LG"=round(chrlen(data)[chr],1),.after=chr)%>%
  mutate("Len of QTL"=round(ci.high-ci.low,1),.after="Len of LG")%>%
  mutate("Flanking markers"=paste0(chr,"_m",find.marker(data,chr,ci.low),"-",chr,"_m",find.marker(data,chr,ci.high)))%>%
  mutate("Central marker"=paste0(chr,"_m",find.marker(data,chr,pos)))%>%
  mutate("Pval"=paste0(pval,if_else(pval<0.63,"*",""),if_else(pval<0.1,"*",""), if_else(pval<0.05,"*","")))%>%
  select(!c(pval,ci.low,ci.high))%>%
  rename("QTL's LG"=chr)%>%
  mutate("No. Inds/% phenotyped"=paste0(nind(data)*su[find.pheno(data,Trait)]," ind / ",round(100*su[find.pheno(data,Trait)],1),"%"),.after=Trait)%>%
  mutate("pos"=round(pos,1))%>%
  mutate("lod"=round(lod,1))
qtldf.aq
```

Customized merged interaction plots

```{r}
mergedIntp<-function(x=1,bin=F,sc=scan1,perm=scan1perm,rit.inx=NULL,mf=c(2,3),l=1,LETTERs=T,second=F,first=F){
  rit<-1 #legend on the right
  if (!second) par(mfrow=mf) #if plot is second then letter is continuous
  for (i in x) {
    (sc1<-summary(sc,perm=perm,alpha=0.63,lodcolumn=i)[,c(1:2,2+i)])
    (r<-nrow(sc1))
    sorted<-sc1%>%arrange(desc(across(3)))
    (mn<-find.marker(data,chr=sorted[,1],pos=sorted[,2]))
    for(j in 1:(r-1)){
      for(k in (j+1):r){
        if (bin) i<-i+n
        effectplot(data,pheno.col=i,mname1=mn[j],mname2=mn[k],main="",ylab=paste("Ave. phenotype:",phenames(data)[i]),xlab=paste0(sorted[k,1],"_m",mn[k]),add.legend = F)
        if(rit%in%rit.inx) lpos<-"topright" else lpos<-"topleft"
        legend(lpos,c("AA","AB","BB"),lty=1,pch=1,col=c("black","red","blue"),bty="n",inset=c(0.01,0))
        a<-par("usr")
        x.leg <- a[1] * 0.05 + a[2] * 0.75
        y.leg <- a[4] - diff(a[3:4]) * 0.05
        if(rit%in%rit.inx) tpos<-2 else tpos<-NULL
        text(x.leg,y.leg,pos=tpos,paste0(sorted[j,1],"_m",mn[j]))
        if (LETTERs)mtext(LETTERS[l],adj=0,cex=1.1);l<-l+1;if(l>6)l<-1
        rit<-rit+1}
    }
  }
  if (!first) par(mfrow=c(1,1))
}
```

Plot marker interactions

```{r}
data<-sim.geno(data,step=3)
x<-(1:6)[1:6%in%find.pheno(data,intpPhen)]
mergedIntp(x,rit.inx=c(3,6,8,12,15,21,22,26,35))
```

```{r eval=FALSE, include=FALSE}
x<-(7:11)[7:11%in%find.pheno(data,intpPhen)]
mergedIntp(x,first=T,rit.inx=1,mf=c(2,2))
#binary
x<-((12:13)[12:13%in%find.pheno(data,intpPhen.bin)])-n
mergedIntp(x[1],T,scan1.bin,scan1perm.bin,second=T,l=2)
mergedIntp(x[2],T,scan1.bin,scan1perm.bin,LETTERs=F)
```

QTL pairs summary

```{r}
#normal
c.thr1<-list()
for(i in 1:n){
  (thr1<-summary(scan2, perms=scan2perm, alpha=0.2,lodcolumn=i,pvalues=T))
  if(i==1){c.thr1[[phenames(data)[i]]]<-thr1
  }else c.thr1[[phenames(data)[i]]]<-thr1
}
for (i in 1:length(c.thr1)){
  if(nrow(c.thr1[[i]])>0&&colnames(c.thr1[[i]])[1]!="Trait"){
    c.thr1[[i]]<-cbind.data.frame(Trait=names(c.thr1[i]),c.thr1[[i]])
  }
}
thr1df<-do.call(rbind.data.frame,c(c.thr1,make.row.names=F))
#binary
c.thr2<-list()
for(i in 1:2){
  (thr2<-summary(scan2.bin, perms=scan2perm.bin, alpha=0.2,lodcolumn=i,pvalues=T))
  if(i==1){c.thr2[[phenames(data)[i+n]]]<-thr2
  }else c.thr2[[phenames(data)[i+n]]]<-thr2
}
for (i in 1:length(c.thr2)){
  if(nrow(c.thr2[[i]])>0&&colnames(c.thr2[[i]])[1]!="Trait"){
    c.thr2[[i]]<-cbind.data.frame(Trait=names(c.thr2[i]),c.thr2[[i]])
  }
}
thr2df<-do.call(rbind.data.frame,c(c.thr2,make.row.names=F))
thr<-rbind(thr1df,thr2df)
```

QTL pairs

```{r}
#normal
data<-calc.genoprob(data,2,map.function="kosambi")
qtlist<-summary(scan1,perms=scan1perm,format="tabByCol",alpha=0.63,ci.function="bayesint")
sc2thr1<-summary(scan2perm,alpha=0.2)
#rearrange the threshold list
th<-vector('list',5)
for(j in 1:5){th[[j]]<-t(sc2thr1[[j]])}
m<-do.call('cbind',th)
dimnames(m)<-list(phenames(data)[1:n],names(sc2thr1)[1:5])
out.ap<-list();qtlpairs<-list();s.fq<-list();out.fq<-list()
for(i in 1:n){
  p<-phenames(data)[i]
  if(length(qtlist[[p]][,1])>0){
    out.ap[[p]]<-addpair(data,qtlist[[p]][,1],p,rqtl[[p]],method="hk",verbose=T)
    qtlpairs[[p]]<-summary(out.ap[[p]],thresholds=m[p,])
    s.fq[[p]]<-summary(out.fq[[p]]<- fitqtl(data,p,rqtl[[p]],method="hk",get.ests=T))
  }
}
#binary
data<-calc.genoprob(data,10,map.function="kosambi")
qtlist.bin<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=0.63,ci.function="bayesint")
sc2thr2<-summary(scan2perm.bin,alpha=0.2)
#rearrange the threshold list
th<-vector('list',5)
for(j in 1:5){th[[j]]<-t(sc2thr2[[j]])}
m.bin<-do.call('cbind',th)
dimnames(m.bin)<-list(phenames(data)[n+1:2],names(sc2thr2)[1:5])
out.ap.bin<-list();rqtl2.bin<-list()
for(i in 1:2){
  p<-phenames(data)[i+n]
  if(length(qtlist.bin[[p]][,1])>0){
    q<-rbind(qtlist.bin[[p]][,-c(3,4)],s.aq[[p]])
    rqtl2.bin[[p]]<-refineqtl(data,p,makeqtl(data,q[,1],q[,2],what="prob"),maxit.fitqtl=1e+6,tol=0.05,method="hk",model="binary")
    out.ap.bin[[p]]<-addpair(data,q[,1],p,rqtl2.bin[[p]],maxit=1e+6,tol=0.2,method="hk",model="binary",verbose=T)
    qtlpairs[[p]]<-summary(out.ap.bin[[p]],thresholds=m.bin[p,])
    s.fq[[p]]<-summary(out.fq[[p]]<- fitqtl(data,p,rqtl2.bin[[p]],maxit=1e+6,tol=0.01,method="hk",model="binary",get.ests=T))
  }
}
for (i in 1:length(qtlpairs)){
if(nrow(qtlpairs[[i]])>0 && names(qtlpairs[[i]])[1]!="Trait"){
  qtlpairs[[i]]<-cbind.data.frame(Trait=names(qtlpairs[i]),qtlpairs[[i]])
  }
}
qtlpairsdf<-do.call(rbind.data.frame,c(qtlpairs,make.row.names=F))
#interacting QTL
qtlpairsdf<-c(qtlpairsdf,thr=m[qtlpairsdf[,1],])
```

No interactive QTL pair was found in the add pair scan.

# Linkage Map view

```{r}
alp<-0.63
colorlist<-RColorBrewer::brewer.pal(8,"Set1")
```

Genetic map to pdf

```{r}
qtldf_initial<-\(){ 
  # make a df to pass qtl info
  qtldf <- data.frame(
    chr = character(),
    qtl = character(),
    so = numeric(),
    si = numeric(),
    ei = numeric(),
    eo = numeric(),
    col = character(),
    stringsAsFactors = F
  )
  return(qtldf)
}
outfile<-file.path("results/basil_linkage_map.pdf")
main<-"Basil Genetic Map"
qtldf<-qtldf_initial()

setting<-list(mapthis=data,outfile=outfile,main=main,ruler=T,maxnbrcolsfordups=2,dupnbr=T,lg.col='lightblue1',lgw=0.15,labdist=0.15,lgperrow=3)
do.call(lmv.linkage.plot,setting)
```

Anthocyanin QTL map

```{r}
qtldf<-qtldf_initial()
for (i in 1:6) {
  (qtls<-summary(scan1,perms=scan1perm,alpha=alp,lodcolumn=i)[,c(1:2,2+i)])
  if(nrow(qtls)>0){
    for (j in 1:nrow(qtls)) {
      (bay <-bayesint(scan1[,c(1:2,2+i)],chr=qtls$chr[j]))
      qtldf <- rbind(qtldf,
                     data.frame(
                       chr = qtls$chr[j],
                       qtl = colnames(bay)[3],
                       so = bay$pos[1],
                       si = bay$pos[2],
                       ei = bay$pos[2],
                       eo = bay$pos[3],
                       col=colorlist[(i+1)]))
    }
  }
}
outfile<-file.path("results/basil_QTLs.anthocyanin.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Anthocyanin (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main,qtldf=qtldf))
do.call(lmv.linkage.plot,setting)
```

Fusarium QTL map

```{r}
qtldf<-qtldf_initial()
for (i in 7:10) {
  qtls<-summary(scan1,perms=scan1perm,alpha=alp,lodcolumn=i)[,c(1:2,2+i)]
  if(nrow(qtls)>0){
    for (j in 1:nrow(qtls)) {
      (bay <-bayesint(scan1[,c(1:2,2+i)],chr=qtls$chr[j]))
      qtldf <- rbind(qtldf,
                     data.frame(
                       chr = qtls$chr[j],
                       qtl = colnames(bay)[3],
                       so = bay$pos[1],
                       si = bay$pos[2],
                       ei = bay$pos[2],
                       eo = bay$pos[3],
                       col=colorlist[(i-5)]))
    }
  }
}
i<-1
(qtls<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=alp,ci.function="bayesint"))
(p<-phenames(data)[i+n])
for (j in 1:(nrow(qtls[[p]])+1)) {
  if(j==3){
    qtls[[p]]<-rbind(qtls[[p]],cbind(s.aq[[p]][,-3],`ci.low`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[1,2],`ci.high`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[3,2],`lod`=s.aq[[p]][,3]))
  }
  qtldf<-rbind(qtldf,
               data.frame(
                 chr = qtls[[p]]$chr[j],
                 qtl = p,
                 so = qtls[[p]]$`ci.low`[j],
                 si = qtls[[p]]$pos[j],
                 ei = qtls[[p]]$pos[j],
                 eo = qtls[[p]]$`ci.high`[j],
                 col=colorlist[7]))
}
outfile<-file.path("results/basil_QTLs.fusarium.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Fusarium (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main))
setting$qtldf<-qtldf
do.call(lmv.linkage.plot,setting)
```

Downy Mildew QTL map

```{r}
qtldf<-qtldf_initial()
i<-2
(qtls<-summary(scan1.bin,perms=scan1perm.bin,format="tabByCol",alpha=alp,ci.function="bayesint"))
p<-phenames(data)[i+n]
for (j in 1:(nrow(qtls[[p]])+1)) {
  if(j==3)qtls[[p]]<-rbind(qtls[[p]],cbind(s.aq[[p]][,-3],`ci.low`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[1,2],`ci.high`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[3,2],`lod`=s.aq[[p]][,3]))
  qtldf<-rbind(qtldf,
               data.frame(
                 chr = qtls[[p]]$chr[j],
                 qtl = p,
                 so = qtls[[p]]$`ci.low`[j],
                 si = qtls[[p]]$pos[j],
                 ei = qtls[[p]]$pos[j],
                 eo = qtls[[p]]$`ci.high`[j],
                 col=colorlist[8]))
}
outfile<-file.path("results/basil_QTLs.BDM.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Downy Mildew (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main))
setting$qtldf<-qtldf
do.call(lmv.linkage.plot,setting)
```

The pdf files are in the results folder.

There is more for that project that I have investigated, such as heatmap of scantwo output, multiple qtl model analysis, larger genome-based dataset (5k markers) with the same workflow but takes more computational resources, alignment of the genetic and genomic (large) map side by side anchored by shared markers, comparing the qtl from both maps. I didn't add it to this notebook because u can see it's quiet much already.