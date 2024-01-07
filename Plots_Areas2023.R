##############################
#### Additional plots ######

#working directory
dire <- "~/GitHub/ENSO_analysis"
setwd(dire)

y.folder <- "C_Isopara95-19New"

##New segment areas ##

#reading data
dat_mat <- read.csv(paste0(dire,"/",y.folder,"/Spcomp_ALL9520.csv"))

#loading function to segment areas
source(paste0(dire,"/toAreaIsoparaSpComp.R")) #for isoparalitoras areas

dat_segm <- toAreaIsoparaSpComp(data=dat_mat[,-(1:6)],lon=dat_mat$Lon_1,
                                           lat=dat_mat$Lat_1,dire=dire)
dat_segm$coord[,2:3] %>% plot
dat_segm$coord$IsoCode[1]

segm_code.area <- dat_segm$coord
segm_code.area$SegmCode <- dat_segm$geog$SiteLD

#isoparalitoras areas 
isopara <- read.csv(paste0(dire,"/areas-isoparalitorales.csv"))
#information of isoparalitoral area
load(paste0(dire,"/infoAreasIsop.RData"))

in_code <- dat_segm$coord$IsoCode
isopara$segm <- NA
for (i in seq_along(in_code)) {
  a1 <- strsplit(dat_segm$coord$IsoCode[i],"-") %>% unlist
  out.iso <- which(isopara$area %in% a1)
  isopara$segm[out.iso] <- rep(dat_segm$geog$SiteLD[i],length(out.iso)) 
}

#eliminate NA values - only segments in study area
segm_area <- isopara[-which(is.na(isopara$segm)),]

seg1 <- unique(segm_area$segm)

lat.5 <- segm_area %>% filter(lat %in% -5)
lat.6 <- segm_area %>% filter(lat %in% -6)
lat.7 <- segm_area %>% filter(lat %in% -7)
lat.8 <- segm_area %>% filter(lat %in% -8)
lat.9 <- segm_area %>% filter(lat %in% -9)
lat.10 <- segm_area %>% filter(lat %in% -10)
lat.11 <- segm_area %>% filter(lat %in% -11)
lat.12 <- segm_area %>% filter(lat %in% -12)
lat.13 <- segm_area %>% filter(lat %in% -13)
lat.14 <- segm_area %>% filter(lat %in% -14)
lat.15 <- segm_area %>% filter(lat %in% -15)
lat.16 <- segm_area %>% filter(lat %in% -16)

require(mapdata)
require(ggspatial)
require(sf)
#peruMap = map_data('worldHires', region = 'Peru')
peruMap <- st_as_sf(maps::map("world",plot=FALSE,fill=TRUE,region="Peru"))

#puertos
ports <- read.csv(paste0(dire,"/puertos1.csv"))
#loading coast line
coast_line <- read.table(file=paste0(dire,"/linea_costa.txt"))
coast_line1 <- coast_line %>% filter(V2 >= -16 & V2 <= -5)
#loading shelfbreak line data
shelf_line <- read.csv("ShelfBreakPosition.csv")

segm_plot <- ggplot() +
  geom_point(data=segm_area ,aes(x=lon,y=lat),
             size=0.025,shape=16,colour="black") +
  geom_point(data=shelf_line,aes(x=lon,y=lat),size=0.025,colour="blue",shape=16)+
  geom_line(data=lat.5,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.6,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.7,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.8,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.9,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.10,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.11,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.12,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.13,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.14,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.15,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_line(data=lat.16,aes(x=lon,y=lat),linewidth=.35,colour="black") +
  geom_sf(data=peruMap) +
  coord_sf(ylim=c(-16,-4),xlim=c(-83,-74)) +
  #blank() +
 # north(peruMap) +
  annotation_scale() +
  # geom_polygon(data = peruMap, aes(long, lat, group = group),
  #              fill = 8, color="black") +
  # coord_cartesian(ylim=c(-16,-3),xlim=c(-83,-74)) +
  scale_y_continuous(breaks=seq(-16,-4,4),
                     labels = c("16°S","12°S","8°S","4°S")) +
  scale_x_continuous(breaks=seq(-83,-74,3),
                     labels = c("83°W","80°W","77°W","74°W")) +
  geom_text(data=ports[c(7,9,12,22),],aes(x=lon+0.7,y=lat,label=Puertos),
            size=2) +
  geom_text(data=ports[26,],aes(x=lon+0.4,y=lat,label=Puertos),
            size=2) +
  geom_text(data=ports[3,],aes(x=lon+0.3,y=lat,label=Puertos),
            size=2) +
  theme(panel.background = element_rect(fill = NA,
                                        color="black"),
        panel.border = element_rect(fill = NA,
                                    color="black"),
        axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.title.y = element_blank())

#save(segm_plot,file=paste0(dire,"/",y.folder,"/SegmPlotVal.RData"))


jpeg(filename = paste0(dire,"/",y.folder,"/Plots/","Segments_areas.jpeg"),
     width = 85, height = 100, units = "mm", res = 500)
segm_plot 
dev.off()

#####
### Plotting sites by ENSO phases ####
##adding LCBD values####

require(mapdata)

peruMap = map_data('worldHires', region = 'Peru')
# phasesDat <- list.files(path=paste0(dire,"/New_Data/ENSO_phases/",y.folder,
#                                     "/Dat_surveys"))
#ordering phases
phasesDat <- c("summ_surv","spr_surv","en_9798_surv","en_1516_surv",
               "wk.en_surv","st.ln_surv")
#anch - no anch
anch.out <- c("Anch","NoAnch","paAnch")


sites_ENSOphase <- list()
sites_ENSOphase.1 <- list()
BD.total <- list()
BD.total1 <- list()
AddLon_plot <- c(0.2,3.7,7.2,10.7,14.2,17.7)

for (a in seq_along(anch.out)) {
  out <- 1
  
  #name folder
  name.beta <- paste0(dire,"/",y.folder,"/BetaDiv/",
                             anch.out[a])
  name.cclust <- paste0(dire,"/",y.folder,"/SpatialClust/",
                        anch.out[a]) 
  for (b in seq_along(phasesDat)) {
    #loading data
    betadivDat <- load(paste0(name.beta,"_",phasesDat[b],".RData")) %>% get()
    cclustDat <- load(paste0(name.cclust,"_",phasesDat[b],".RData")) %>% get()
    
    #segment areas data
    segmDat <- load(paste0(dire,"/",y.folder,"/Segm_",
                           phasesDat[b],".RData")) %>% get()
    
    lonlat_sites <- cclustDat$clust$coords %>% as.data.frame()
    lonlat_sites$lonNot <- lonlat_sites$lon - AddLon_plot[b] 
    lonlat_sites$LCBD <- betadivDat$LCBD
    lonlat_sites$pLCBD1 <- betadivDat$pLCBD1
    lonlat_sites$Rich <- betadivDat$Rich
    lonlat_sites$Phase <- rep(phasesDat[b],dim(lonlat_sites)[1])
    lonlat_sites$Site <- segmDat$geog$SiteLD
    lonlat_sites$latG <- segmDat$geog$latG
    lonlat_sites$dc <- segmDat$geog$dc
    lonlat_sites$n_trawl <- segmDat$geog$n_trawl
    lonlat_sites$shelf <- segmDat$geog$shelf
    
    #intervals for Lat and DC
    latGint <- seq(5,15,3)
    dcint <- seq(10,80,20)
    latGIntv <- matrix(NA,nrow=nrow(lonlat_sites),ncol=1)
    dcIntv <- matrix(NA,nrow=nrow(lonlat_sites),ncol=1)
    for (c in 1:4) {
      lat.out <- which(lonlat_sites$latG >= latGint[c] & lonlat_sites$latG < (latGint[c]+3))
      dc.out <- which(lonlat_sites$dc >= dcint[c] & lonlat_sites$dc < (dcint[c]+15))
      latGIntv[lat.out] <- paste0(latGint[c],"-",(latGint[c]+2))
      dcIntv[dc.out] <- paste0(dcint[c],"-",(dcint[c]+10))
      }
    lonlat_sites$latGInt <- latGIntv
    lonlat_sites$dcInt <- dcIntv
    
    sites_ENSOphase[[out]] <- lonlat_sites
    BD.total[[out]] <- betadivDat$beta %>% as.matrix()
    out <- out + 1
  }
  sites_ENSOphase.1[[a]] <- bind_rows(sites_ENSOphase)
  rownames(sites_ENSOphase.1[[a]]) <- NULL
  
  BD.total1[[a]] <- bind_cols(BD.total)
  colnames(BD.total1[[a]]) <- phasesDat
  row.names(BD.total1[[a]]) <- c("SStotal","BDtotal")
}


# sites_ENSOphase.1$Phase <- factor(sites_ENSOphase.1$Phase,
#                                   levels=c("summ_surv","spr_surv","en_9798_surv",
#                                           "en_1516_surv","wk.en_surv",
#                                           "st.ln_surv"))

save(sites_ENSOphase.1,file=paste0(dire,"/",
                     y.folder,"/sites_ENSOphaseAll.RData"))

#Plotting sites with LCBD values #####
#will be for each Anchovy scenario

require(ggpubr)

sites_ENSOphase.1 <- load(paste0(dire,"/",
            y.folder,"/sites_ENSOphaseAll.RData")) %>% get()

#anch - no anch
anch.out <- c("Anch","NoAnch","paAnch")

col1 <- c("#336633","#CC9966","#CC0000","#FF3333","#FF6666","#0033FF")

for (a in seq_along(anch.out)) {
  
  dat <- sites_ENSOphase.1[[a]]
  
  dat$Phase <- factor(dat$Phase,levels=c("summ_surv","spr_surv","en_9798_surv",
                        "en_1516_surv","wk.en_surv","st.ln_surv"),
                      labels = c("N_Summer","N_Spring","EN_9798","EN_1516",
                                 "EN_weak","LN_strong"))
  dat$latGInt <- factor(dat$latGInt,levels=c("5-7","8-10","11-13","14-16"))
  dat$dcInt <- factor(dat$dcInt,levels=c("10-20","30-40","50-60","70-80"))
  
  sitesLCBDpval <-  dat %>% filter(pLCBD1 <= 0.05)
  #plot
  lcbd.plot <-
    ggplot(data = dat,aes(x=shelf,y=LCBD)) +
    facet_grid(cols=vars(Phase),vars(latGInt)) +
    geom_point() +
    geom_smooth() +
    scale_color_manual(values=col1) +
    geom_point(data=sitesLCBDpval,aes(x=dc,y=LCBD),
               color="red") +
    labs(x="Distance to shelf break (NM)",y="LCBD value",
         shape="Latitude (°S)",color="ENSO phase") +
    #guides(color="none") +
    theme(panel.background = element_rect(fill = NA,
                                          color="black"),
          panel.grid.minor = element_line(colour = "gray83", linewidth= 0.1),
          panel.grid.major = element_line(colour = "gray83", linewidth= 0.1),
          axis.text.x = element_text(angle=90,size=5),
          axis.text.y = element_text(size=5),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))
  
  
  jpeg(filename = paste0(dire,"/",y.folder,"/Plots/",
                         anch.out[a],"_LCBDval1.jpeg"),
       width = 170, height = 100, units = "mm", res = 500)
  annotate_figure(lcbd.plot,right="Latitude (°S)",top="ENSO phase")
  dev.off()  
  
}

#saving significant LCBD values - without anchovy
dat <- sites_ENSOphase.1[[1]]
sitesLCBDpval <-  dat %>% filter(pLCBD1 <= 0.05)

write.csv(sitesLCBDpval,file=paste0(dire,"/New_Data/ENSO_phases/",y.folder,
                                    "/sitesLCBDpvalAnch.csv"),
          row.names = FALSE)



##### Plotting spatial clusters ####
#considering anchovy condition

#puertos
ports <- read.csv(paste0(dire,"/puertos1.csv"))
#color spatial clusters
colors=rainbow(5)

for (i in seq_along(anch.out)) {
  name.beta <- paste0(dire,"/",y.folder,"/BetaDiv/",
                      anch.out[i])
  name.cclust <- paste0(dire,"/",y.folder,"/SpatialClust/",
                        anch.out[i]) 
  for (j in seq_along(phasesDat)) {
    #loading data
    betadivDat <- load(paste0(name.beta,"_",phasesDat[j],".RData")) %>% get()
    cclustDat <- load(paste0(name.cclust,"_",phasesDat[j],".RData")) %>% get()
    
    #segment areas data
    segmDat <- load(paste0(dire,"/",y.folder,"/Segm_",
                           phasesDat[j],".RData")) %>% get()
    
    jpeg(filename = paste0(dire,"/",y.folder,"/Plots/",
                           anch.out[i],"_",phasesDat[j],
                           "_spat_cluster.jpeg"),
         width = 170, height = 235, units = "mm", res = 500)
    par(mar=c(1,3,0.5,0.5),mfrow=c(2,2),font=2)
    for (l in 2:5) {
    constr.hclust:::plot.constr.hclust(cclustDat$clust, k=l, links=T, xlab="", 
                                       ylab="",cex=0.1,lwd=2,pch=0,
                                       col=colors,ylim=c(-15.5,-4.5),
                                       xlim=c(-83,-74), axes=FALSE)
    axis(2,las=1,seq(-16,-6,2),labels = c("16°S","14°S","12°S","10°S","8°S","6°S"))
    map("worldHires", add=T, fill=T, col=8)
    legend("bottomleft", lty=1L, lwd=3, col=colors[1:l],
           legend=sprintf("Group %d",1:l), cex=0.9, bty="n")
    #legend(-79.25,-4.5,legend=phase[l],bty="n",cex=1.2)
    text(ports[c(7,9,12,22),2]+0.85,ports[c(7,9,12,22),3],labels=ports[c(7,9,12,22),1],cex=0.75)
    text(ports[26,2]+0.5,ports[26,3],labels=ports[26,1],cex=0.75)
    text(ports[3,2]+0.5,ports[3,3],labels=ports[3,1],cex=0.75)
    box()
    }
    dev.off()
    
    
  }
}

###IndVal values for spatial clusters ####
##adding anchovy condition

require(indicspecies)

#ordering phases
phasesDat <- c("summ_surv","spr_surv","en_9798_surv","en_1516_surv",
               "wk.en_surv","st.ln_surv")
#anch - no anch
anch.out <- c("Anch","NoAnch","paAnch")

#number clusters
kclust <- list(c(3,2,3,3,3,2),c(rep(3,6)),c(2,3,2,3,3,3))
#selected species
sp_sel <- read.csv(paste0("./table_species1.csv"))

for (i in seq_along(anch.out)) {
  #cluster data
  name.cclust <- paste0(dire,"/",y.folder,"/SpatialClust/",anch.out[i])
  #number of clusters
  k <- kclust[[i]]
  
  for (j in seq_along(phasesDat)) {
    #loading cluster data
    cclustDat <- load(paste0(name.cclust,"_",phasesDat[j],".RData")) %>% get()
    #segment areas data
    MatDat <- load(paste0(dire,"/",y.folder,"/",phasesDat[j],
                          "_MatDat.RData")) %>% get()
    
    #indval  
    indval_cluster <- multipatt(MatDat[[i]], cutree(cclustDat$clust,k=k[j]),
                                control = how(nperm = 99999),max.order = k[j])
    ##saving results of significant species
    #significant sp
    (pval.adj2 <- p.adjust(indval_cluster$sign$p.value,"holm"))
    #sp_sign <- indval_cluster$sign[which((indval_cluster$sign$p.value) <= 0.05),]
    sp_sign <- indval_cluster$sign[which(pval.adj2 <= 0.05),]
    sele_A <- indval_cluster$A[which(pval.adj2 <= 0.05),] %>% matrix(nrow=nrow(sp_sign)) %>% 
                          as.data.frame()
    colnames(sele_A) <- paste0("A.",colnames(indval_cluster$A))
    sele_B <- indval_cluster$B[which(pval.adj2 <= 0.05),] %>% matrix(nrow =nrow(sp_sign)) %>% 
                        as.data.frame()
    colnames(sele_B) <- paste0("B.",colnames(indval_cluster$B))
    sel_sp_indval <- mutate(sp_sign,A=sele_A,B=sele_B)
    
    write.csv(sel_sp_indval,file=paste0(dire,"/",y.folder,
                                        "/IndVal/",anch.out[i],"/IndVal_",
                                        phasesDat[j],".csv"),row.names = T)
      }
  }





###Plotting FINAL SPATIAL PLOT FOR ALL PHASES #####
##adding anchovy condition

name.phase <- c("N_Summer","N_Spring","EN_9798","EN_1516",
                "EN_weak","LN_strong")

for (i in seq_along(anch.out)) {
  name.beta <- paste0(dire,"/New_Data/ENSO_phases/",y.folder,"/BetaDiv/",
                      anch.out[i])
  name.cclust <- paste0(dire,"/New_Data/ENSO_phases/",y.folder,"/SpatialClust/",
                        anch.out[i]) 
 
  jpeg(filename = paste0(dire,"/New_Data/ENSO_phases/",y.folder,"/Plots/",
                         anch.out[i],"_spat_cluster-SillVal.jpeg"),
       width = 270, height = 235, units = "mm", res = 500)
  par(mar=c(1,3,0.5,0.5),mfrow=c(2,3),font=2)
  for (j in seq_along(phasesDat)) {
    #loading data
    cclustDat <- load(paste0(name.cclust,"_",phasesDat[j],".RData")) %>% get()
    
    #plot
    constr.hclust:::plot.constr.hclust(cclustDat$clust, k=cclustDat$k[1], links=T, xlab="", 
                                       ylab="",cex=0.1,lwd=2,pch=0,
                                       col=colors,ylim=c(-15.5,-4.5),
                                       xlim=c(-83,-74), axes=FALSE)
    axis(2,las=1,seq(-16,-6,2),labels = c("16°S","14°S","12°S","10°S","8°S","6°S"))
    map("worldHires", add=T, fill=T, col=8)
    legend("bottomleft", lty=1L, lwd=3, col=colors[1:cclustDat$k[1]],
           legend=sprintf("Group %d",1:cclustDat$k[1]), cex=0.9, bty="n")
    legend(-79.25,-4.5,legend=name.phase[j],bty="n",cex=1.2)
    text(ports[c(7,9,12,22),2]+0.85,ports[c(7,9,12,22),3],labels=ports[c(7,9,12,22),1],cex=0.75)
    text(ports[26,2]+0.5,ports[26,3],labels=ports[26,1],cex=0.75)
    text(ports[3,2]+0.5,ports[3,3],labels=ports[3,1],cex=0.75)
    box()
  }
  dev.off()
}


##Defining number of clusters k=3
for (i in seq_along(anch.out)) {
  name.beta <- paste0(dire,"/New_Data/ENSO_phases/",y.folder,"/BetaDiv/",
                      anch.out[i])
  name.cclust <- paste0(dire,"/New_Data/ENSO_phases/",y.folder,"/SpatialClust/",
                        anch.out[i]) 
  
  #number of clusters
  k <- kclust[[i]]
  
  jpeg(filename = paste0(dire,"/New_Data/ENSO_phases/",y.folder,"/Plots/",
                         anch.out[i],"_spat_clusterFINAL.jpeg"),
       width = 270, height = 235, units = "mm", res = 500)
  par(mar=c(1,3,0.5,0.5),mfrow=c(2,3),font=2)
  for (j in seq_along(phasesDat)) {
    #loading data
    cclustDat <- load(paste0(name.cclust,"_",phasesDat[j],".RData")) %>% get()
    
    #plot
    constr.hclust:::plot.constr.hclust(cclustDat$clust, k=k[j], links=T, xlab="", 
                                       ylab="",cex=0.1,lwd=2,pch=0,
                                       col=colors,ylim=c(-15.3,-4.5),
                                       xlim=c(-82.5,-75), axes=FALSE)
    axis(2,las=1,seq(-16,-6,2),labels = c("16°S","14°S","12°S","10°S","8°S","6°S"))
    maps:::map("worldHires", add=T, fill=T, col=8)
    legend("bottomleft",lty=1L, lwd=3, col=colors[1:3],
           legend=sprintf("Group %d",1:k[j]), cex=1.2, bty="n")
    legend(-79.25,-4.5,legend=name.phase[j],bty="n",cex=1.8)
    text(ports[c(7,9,12,22),2]+0.75,ports[c(7,9,12,22),3],labels=ports[c(7,9,12,22),1],cex=0.8)
    text(ports[26,2]+0.35,ports[26,3],labels=ports[26,1],cex=0.8)
    text(ports[3,2]+0.3,ports[3,3],labels=ports[3,1],cex=0.8)
    box()
  }
  dev.off()
}





sitesLCBDpval <- sites_ENSOphase.1 %>% filter(pLCBD1 <= 0.05)

lcbd.plot2 <- ggplot() +
  geom_point(data = sites_ENSOphase.1,aes(x=lonNot,y=lat,fill=Phase,
                                          size=LCBD*100),shape=21,
             colour="black") +
  scale_fill_manual(values=col1) +
  geom_point(data=sitesLCBDpval,aes(x=lonNot,y=lat,size=LCBD*100),
             shape=21,colour="black",fill="black") +
  geom_polygon(data = peruMap, aes(long, lat, group = group),
               fill = 8, color="black") +
  coord_cartesian(ylim=c(-16,-3),xlim=c(-100,-74)) +
  scale_y_continuous(breaks=seq(-16,-4,4),
                     labels = c("16°S","12°S","8°S","4°S")) +
  theme(panel.background = element_rect(fill = NA,
                                        color="black"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank()) 


###plotting points with number trawls

sites_ENSOphase.1$n_trawl2 <- (sqrt(sites_ENSOphase.1$n_trawl))*0.125

trawls.plot <- ggplot() +
  geom_point(data = sites_ENSOphase.1,aes(x=lonNot,y=lat,fill=Phase,
                                          size=(n_trawl2)),shape=21,
             colour="black",alpha=0.8) +
  scale_fill_manual(values=col1) +
  scale_size_continuous(breaks=seq(0.25,1,0.25),labels=c("4","16","36","64")) +
  geom_polygon(data = peruMap, aes(long, lat, group = group),
               fill = 8, color="black") +
  coord_cartesian(ylim=c(-16,-3),xlim=c(-99.5,-75.5)) +
  scale_y_continuous(breaks=seq(-16,-4,4),
                     labels = c("16°S","12°S","8°S","4°S")) +
  geom_text(data=ports[c(7,9,12,22),],aes(x=lon+0.8,y=lat,label=Puertos),
            size=2.5) +
  geom_text(data=ports[26,],aes(x=lon+0.6,y=lat,label=Puertos),
            size=2.5) +
  geom_text(data=ports[3,],aes(x=lon+0.4,y=lat,label=Puertos),
            size=2.5) +
  guides(color="none",fill="none") +
  labs(size="N trawls",fill="ENSO phase") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white", 
                                        colour = "black"),
        legend.position = c(0.06,0.18)) 
  
  
jpeg(filename = paste0(dire,"/New_Data/ENSO_phases/",y.folder,"/Plots/",
                       "Sites_trawlsALL.jpeg"),
     width = 180, height = 100, units = "mm", res = 500)
trawls.plot
dev.off()


### Plotting points by ENSO phase ##

#66 species selected - high occurrence at least 1 time
sp_sel <- read.csv(paste0(dire,"/table_species1.csv"))

##data from segmented areas by ENSO phase
en_1516_surv <- read.csv(file=paste0(dire,"/",y.folder,"/Dat_surveys/en_1516_surv.csv"))
en_9798_surv <- read.csv(file=paste0(dire,"/",y.folder,"/Dat_surveys/en_9798_surv.csv"))
wk.en_surv <- read.csv(file=paste0(dire,"/New_Data/ENSO_phases/",y.folder,
                                   "/Dat_surveys/wk.en_surv.csv"))
st.ln_surv <- read.csv(file=paste0(dire,"/New_Data/ENSO_phases/",y.folder,
                                   "/Dat_surveys/st.ln_surv.csv"))
summ_surv <- read.csv(file=paste0(dire,"/New_Data/ENSO_phases/",y.folder,
                                  "/Dat_surveys/summ_surv.csv"))
spr_surv <- read.csv(file=paste0(dire,"/New_Data/ENSO_phases/",y.folder,
                                "/Dat_surveys/spr_surv.csv"))

plot(en_1516_surv[,3:2])

