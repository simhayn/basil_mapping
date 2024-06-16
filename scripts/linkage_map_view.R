####linkage map view####
#clean workspace
rm(list=ls())
# Install and load the LinkageMapView package if not already installed
if (!requireNamespace("LinkageMapView", quietly = T)) {
  install.packages("LinkageMapView")
}
library(LinkageMapView)
library(qtl)

alp<-0.63
n<-11
colorlist<-RColorBrewer::brewer.pal(8,"Set1")

data<-read.cross("csv",".","data.csv")
load("./genome_scans.Rdata")

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

#plot map without qtl to pdf file
outfile<-file.path("results/basil_linkage_map.pdf")
main<-"Basil Genetic Map"
qtldf<-qtldf_initial()

setting<-list(mapthis=data,outfile=outfile,main=main,ruler=T,maxnbrcolsfordups=2,dupnbr=T,lg.col='lightblue1',lgw=0.15,labdist=0.15,lgperrow=3)
do.call(lmv.linkage.plot,setting)

#normal
####anthocyanin####
{
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
}

outfile<-file.path("results/basil_QTLs.anthocyanin.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Anthocyanin (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main,qtldf=qtldf))
do.call(lmv.linkage.plot,setting)


####fusarium####
{
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
                         col=colorlist[(i-5)]
                       ))
      }
    }
  }
}

#binary
#FOB3_bin
{
  # Load the add qtl output  
  load("aq.Rdata")
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
}
outfile<-file.path("results/basil_QTLs.fusarium.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Fusarium (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main))
setting$qtldf<-qtldf
do.call(lmv.linkage.plot,setting)


####BDM####
{
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
                   col=colorlist[8])
    )
  }
}

outfile<-file.path("results/basil_QTLs.BDM.pdf")
(mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
(main<-paste0("Basil Genetic Map + QTLs for Downy Mildew (",paste0(mapthese,collapse = ","),")")) 
setting<-modifyList(setting,list(outfile=outfile,mapthese=mapthese,main=main))
setting$qtldf<-qtldf
do.call(lmv.linkage.plot,setting)
