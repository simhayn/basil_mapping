####linkage map view####
library(LinkageMapView)
alp<-0.63
colorlist<-brewer.pal(8,"Set1")

outfile = file.path("basil_denovo.pdf")
lmv.linkage.plot(
  basil,
  outfile,
  main="Basil Denovo Genetic Map",
  ruler = T,
  maxnbrcolsfordups = 2,
  dupnbr = T,
  lg.col = "lightblue1",
  lgw=0.15,
  labdist=0.15,
  lgperrow = 3
  )

#normal
#antho
{
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
  for (i in 1:6) {
      (qtls<-summary(scan1,perms=scan1perm,alpha=alp,lodcolumn=i)[,c(1:2,2+i)])
      if(nrow(qtls)>0){
        for (j in 1:nrow(qtls)) {
          (bay <-bayesint(scan1[,c(1:2,2+i)],chr=qtls$chr[j]))
          qtldf <- rbind(qtldf,
                         data.frame(
                           chr = qtls$chr[j],
                           qtl = phenames(basil)[i],
                           so = bay$pos[1],
                           si = bay$pos[2],
                           ei = bay$pos[2],
                           eo = bay$pos[3],
                           col=colorlist[(i+1)]
                         ))
        }}}
  }

{  
  outfile = file.path("basil_denovo_QTLs.anthocyanin1.pdf")
  (mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
  (mapthese1<-mapthese[5:6])
  (tit<-paste0("Basil Denovo QTLs Anthocyanin (",paste0(mapthese,collapse = ","),")")) 
  lmv.linkage.plot(
    basil,
    outfile,
    mapthese1,
    ruler = T,
    maxnbrcolsfordups = 2,
    dupnbr = T,
    lg.col = "lightblue1",
    main =  tit,
    qtldf = qtldf,
    lgw=0.15,
    labdist=0.15,
    lgperrow = 2
  )}

#####

#fusarium
{
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
  for (i in 7:10) {
    qtls<-summary(scan1,perms=scan1perm,alpha=alp,lodcolumn=i)[,c(1:2,2+i)]
    if(nrow(qtls)>0){
      for (j in 1:nrow(qtls)) {
        (bay <-bayesint(scan1[,c(1:2,2+i)],chr=qtls$chr[j]))
        qtldf <- rbind(qtldf,
                       data.frame(
                         chr = qtls$chr[j],
                         qtl = phenames(basil)[i],
                         so = bay$pos[1],
                         si = bay$pos[2],
                         ei = bay$pos[2],
                         eo = bay$pos[3],
                         col=colorlist[(i-5)]
                       ))
      }}}
}
  
  #binary
  #FOB3_bin
{
  i<-1
  (qtls<-summary(scan1.bin,perms=scan1perm.bin,
                 format="tabByCol",alpha=alp,ci.function="bayesint"))
  p<-phenames(basil)[i+n]
  for (j in 1:(nrow(qtls[[p]])+1)) {
    if(j==3)qtls[[p]]<-rbind(qtls[[p]],
                             cbind(s.aq[[p]][,-3],
                                   `ci.low`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[1,2],
                                   `ci.high`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[3,2],
                                   `lod`=s.aq[[p]][,3])
    )
    qtldf<-rbind(qtldf,
                 data.frame(
                   chr = qtls[[p]]$chr[j],
                   qtl = p,
                   so = qtls[[p]]$`ci.low`[j],
                   si = qtls[[p]]$pos[j],
                   ei = qtls[[p]]$pos[j],
                   eo = qtls[[p]]$`ci.high`[j],
                   col=colorlist[7]
                 ))
  }
}



{  
  outfile = file.path("basil_denovo_QTLs.fusarium.pdf")
  (mapthese<-paste0("LG",sort(unique(as.numeric(qtldf$chr)))))
  (tit<-paste0("Basil Denovo QTLs Fusarium (",paste0(mapthese,collapse = ","),")"))
  lmv.linkage.plot(
    basil,
    outfile,
    mapthese,
    ruler = T,
    maxnbrcolsfordups = 2,
    dupnbr = T,
    lg.col = "lightblue1",
    main =  tit,
    qtldf = qtldf,
    lgw=0.15,
    labdist=0.15,
    lgperrow = 2
  )}

#BDM
{
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
  i<-3
  (qtls<-summary(scan1.bin,perms=scan1perm.bin,
                      format="tabByCol",alpha=alp,ci.function="bayesint"))
  p<-phenames(basil)[i+n]
  for (j in 1:(nrow(qtls[[p]])+1)) {
    if(j==3)qtls[[p]]<-rbind(qtls[[p]],
                             cbind(s.aq[[p]][,-3],
                                   `ci.low`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[1,2],
                                   `ci.high`=bayesint(out.aq.bin[[p]],s.aq[[p]][,1])[3,2],
                                   `lod`=s.aq[[p]][,3])
                             )
    qtldf<-rbind(qtldf,
                 data.frame(
                   chr = qtls[[p]]$chr[j],
                   qtl = p,
                   so = qtls[[p]]$`ci.low`[j],
                   si = qtls[[p]]$pos[j],
                   ei = qtls[[p]]$pos[j],
                   eo = qtls[[p]]$`ci.high`[j],
                   col=colorlist[8]
                 ))
  }
}

{
  outfile = file.path("basil_denovo_QTLs.BDM1.pdf")
  (mapthese<-unique(qtldf$chr)[1])
  (tit<-paste0("Basil Denovo QTLs BDM (",paste0(mapthese,collapse = ", "),")"))
  lmv.linkage.plot(
    basil,
    outfile,
    mapthese,
    ruler = T,
    maxnbrcolsfordups = 2,
    dupnbr = T,
    lg.col = "lightblue1",
    main =  tit,
    qtldf = qtldf,
    lgw=0.15,
    labdist=0.15,
    lgperrow = 2
  )}

