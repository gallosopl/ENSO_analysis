#################
## Temporal Beta diversity index ###

#using segment areas as site point
#Pair of TBI that will be calculated:
#1#T1: Neutral summer - T2: Neutral Spring
#2#T1: EN 9798 - T2: EN 1516
#3#T1: EN 1516 - T2: LN strong


require(adespatial)
require(dplyr)
require(tidyverse)
require(ggplot2)
require(ggpubr)
require(gridExtra)
require(mapdata)

###
rm(list = ls())
#working directory
dire <- "~/GitHub/ENSO_analysis"
setwd(dire)

y.folder <- "C_Isopara95-19New"

phasesDat <- list.files(path=paste0(dire,"/",y.folder,"/Dat_surveys"))
phase.folder <- strsplit(phasesDat,".csv",fixed=TRUE) %>% as.character()

#67 species selected 
sp_sel <- read.csv(paste0(dire,"/table_species1.csv"))

#loading data from segment sites 
#with latitude and dc intervals
site_ENSOphases <- load(paste0(dire,"/",y.folder,
                 "/sites_ENSOphaseAll.RData")) %>% get
site_ENSOphasesL <- site_ENSOphases[[1]] %>% select(lon,lat,Site,latG,dc,
                                                    latGInt,dcInt) %>% 
                                                              distinct(.,) %>%
                                                    as.data.frame()
site_ENSOphasesL$latGInt <- factor(site_ENSOphasesL$latGInt,
                                   levels=c("5-7","8-10","11-13","14-16"))
site_ENSOphasesL$dcInt <- factor(site_ENSOphasesL$dcInt,
                                 levels=c("10-20","30-40","50-60","70-80"))

  
##first comparison: #### 
#Neutral Summer - Neutral Spring

segmSumm <- load(paste0(dire,"/",y.folder,
                         "/Segm_summ_surv.RData")) %>% get
  
segmSpri <- load(paste0(dire,"/",y.folder,
                        "/Segm_spr_surv.RData")) %>% get

sameSites <- as.data.frame(table(c(segmSumm$geog$SiteLD, segmSpri$geog$SiteLD))) %>% 
                 filter(.,Freq %in% 2)
Summ.T1 <- segmSumm$sp_comp[which(segmSumm$geog$SiteLD %in% sameSites$Var1),
                            which(colnames(segmSumm$sp_comp) %in%
                                    sp_sel$T.st)]
Spri.T2 <- segmSpri$sp_comp[which(segmSpri$geog$SiteLD %in% sameSites$Var1),
                            which(colnames(segmSpri$sp_comp) %in%
                                    sp_sel$T.st)]

#noAnchovy
Summ.T1NoAnch <- Summ.T1[,-(which(colnames(Summ.T1) %in% "E.ringens"))]
Spri.T2NoAnch <- Spri.T2[,-(which(colnames(Spri.T2) %in% "E.ringens"))]

#with anchovy
TBI.summ.spri <- TBI(Summ.T1 %>% log1p(),Spri.T2 %>% log1p(),
                     method="%diff",pa.tr=FALSE,
                     BCD=TRUE,test.BC=TRUE,test.t.perm=TRUE,nperm=9999)
site.loss <- TBI.summ.spri$BCD.mat %>% filter(.,Change %in% "-  ")
site.gain <- TBI.summ.spri$BCD.mat %>% filter(.,Change %in% "+  ")

summary(TBI.summ.spri$TBI)

#arranging data for ggplot
dat.TBI.SummSpri <- TBI.summ.spri$BCD.mat %>% mutate(segmSumm$geog %>% filter(SiteLD %in% 
                                    sameSites$Var1) %>% select(SiteLD))
in.site <- match(dat.TBI.SummSpri$SiteLD, site_ENSOphasesL$Site)
dat.TBI.SummSpri1 <- mutate(dat.TBI.SummSpri,site_ENSOphasesL$lon[in.site],
                            site_ENSOphasesL$lat[in.site],site_ENSOphasesL$latG[in.site],
                            site_ENSOphasesL$dc[in.site],site_ENSOphasesL$latGInt[in.site],
                            site_ENSOphasesL$dcInt[in.site])
colnames(dat.TBI.SummSpri1) <- c("B","C","D","Change","SiteLD","lon","lat","latG",
                                 "dc","latGInt","dcInt") 
dat.TBI.SummSpri1$Code_change <- rep("SummSpri",nrow(dat.TBI.SummSpri))

##looking into species biomass change
sp.table1 <- cbind(apply(Summ.T1,2,mean),apply(Spri.T2,2,mean))
sp.tableSummSpri <- cbind(sp.table1,(sp.table1[,2]-sp.table1[,1])) %>% 
                    as.data.frame
sp.tableSummSpri$sp_names <- rownames(sp.tableSummSpri)
sp.tableSummSpri$code_change <- rep("SummSpri",nrow(sp.tableSummSpri))

#x11()
plot(TBI.summ.spri, pch.loss=19, pch.gain=15,
     main="B-C plot",xlim=c(0,1),ylim=c(0,1),diam=TRUE,
     silent = FALSE)


##second comparison: ##### 
#EN 1997-98 - EN 2015-16

segmEN97 <- load(paste0(dire,"/",y.folder,
                        "/Segm_en_9798_surv.RData")) %>% get

segmEN15 <- load(paste0(dire,"/",y.folder,
                        "/Segm_en_1516_surv.RData")) %>% get

sameSites2 <- as.data.frame(table(c(segmEN97$geog$SiteLD, segmEN15$geog$SiteLD))) %>% 
                          filter(.,Freq %in% 2)
EN97.T1 <- segmEN97$sp_comp[which(segmEN97$geog$SiteLD %in% sameSites2$Var1),
                            which(colnames(segmEN97$sp_comp) %in%
                                    sp_sel$T.st)]
EN15.T2 <- segmEN15$sp_comp[which(segmEN15$geog$SiteLD %in% sameSites2$Var1),
                            which(colnames(segmEN15$sp_comp) %in%
                                    sp_sel$T.st)]

#noAnchovy
EN97.T1NoAnch <- EN97.T1[,-(which(colnames(EN97.T1) %in% "E.ringens"))]
EN15.T2NoAnch <- EN15.T2[,-(which(colnames(EN15.T2) %in% "E.ringens"))]


#with anchovy
TBI.en97.15 <- TBI(EN97.T1 %>% log1p(),EN15.T2 %>% log1p(),
                     method="%diff",pa.tr=FALSE,
                     BCD=TRUE,test.BC=TRUE,test.t.perm=TRUE,nperm=9999)
site.loss2 <- TBI.en97.15$BCD.mat %>% filter(.,Change %in% "-  ")
site.gain2 <- TBI.en97.15$BCD.mat %>% filter(.,Change %in% "+  ")

#arranging data for ggplot
dat.TBI.en9715 <- TBI.en97.15$BCD.mat %>% mutate(segmEN97$geog %>% filter(SiteLD %in% 
                                              sameSites2$Var1) %>% select(SiteLD))
in.site2 <- match(dat.TBI.en9715$SiteLD, site_ENSOphasesL$Site)
dat.TBI.en9715.1 <- mutate(dat.TBI.en9715,site_ENSOphasesL$lon[in.site2],
                            site_ENSOphasesL$lat[in.site2],site_ENSOphasesL$latG[in.site2],
                            site_ENSOphasesL$dc[in.site2],site_ENSOphasesL$latGInt[in.site2],
                            site_ENSOphasesL$dcInt[in.site2])
colnames(dat.TBI.en9715.1) <- c("B","C","D","Change","SiteLD","lon","lat","latG",
                                "dc","latGInt","dcInt") 
dat.TBI.en9715.1$Code_change <- rep("EN9715",nrow(dat.TBI.en9715))

##looking into species biomass change
sp.table2 <- cbind(apply(EN97.T1,2,mean),apply(EN15.T2,2,mean))
sp.tableEN9715 <- cbind(sp.table2,(sp.table2[,2]-sp.table2[,1])) %>% 
                    as.data.frame
sp.tableEN9715$sp_names <- rownames(sp.tableEN9715)
sp.tableEN9715$code_change <- rep("EN9715",nrow(sp.tableEN9715))

#x11()
plot(TBI.en97.15, pch.loss=19, pch.gain=15,
     main="B-C plot",xlim=c(0,1),ylim=c(0,1),diam=TRUE,
     silent = FALSE)


##third comparison: ##### 
#EN 15-16- LN strong

segmEN1516 <- load(paste0(dire,"/",y.folder,
                        "/Segm_en_1516_surv.RData")) %>% get

segmLNst <- load(paste0(dire,"/",y.folder,
                        "/Segm_st.ln_surv.RData")) %>% get

sameSites3 <- as.data.frame(table(c(segmEN15$geog$SiteLD, segmLNst$geog$SiteLD))) %>% 
  filter(.,Freq %in% 2)
EN1516.T1 <- segmEN1516$sp_comp[which(segmEN1516$geog$SiteLD %in% sameSites3$Var1),
                            which(colnames(segmEN1516$sp_comp) %in%
                                    sp_sel$T.st)]
LNst.T2 <- segmLNst$sp_comp[which(segmLNst$geog$SiteLD %in% sameSites3$Var1),
                            which(colnames(segmLNst$sp_comp) %in%
                                    sp_sel$T.st)]

#with anchovy
TBI.en15.ln <- TBI(EN1516.T1 %>% log1p(),LNst.T2 %>% log1p(),
                   method="%diff",pa.tr=FALSE,
                   BCD=TRUE,test.BC=TRUE,test.t.perm=TRUE,nperm=9999)
site.loss3 <- TBI.en15.ln$BCD.mat %>% filter(.,Change %in% "-  ")
site.gain3 <- TBI.en15.ln$BCD.mat %>% filter(.,Change %in% "+  ")

#arranging data for ggplot
dat.TBI.en15.ln <- TBI.en15.ln$BCD.mat %>% mutate(segmEN1516$geog %>% filter(SiteLD %in% 
                                          sameSites3$Var1) %>% select(SiteLD))
in.site3 <- match(dat.TBI.en15.ln$SiteLD, site_ENSOphasesL$Site)
dat.TBI.en15.ln.1 <- mutate(dat.TBI.en15.ln,site_ENSOphasesL$lon[in.site3],
                           site_ENSOphasesL$lat[in.site3],site_ENSOphasesL$latG[in.site3],
                           site_ENSOphasesL$dc[in.site3],site_ENSOphasesL$latGInt[in.site3],
                           site_ENSOphasesL$dcInt[in.site3])
colnames(dat.TBI.en15.ln.1) <- c("B","C","D","Change","SiteLD","lon","lat","latG",
                                "dc","latGInt","dcInt") 
dat.TBI.en15.ln.1$Code_change <- rep("EN15.LN",nrow(dat.TBI.en15.ln))

##looking into species biomass change
sp.table3 <- cbind(apply(EN1516.T1,2,mean),apply(LNst.T2,2,mean))
sp.tableEN15LN <- cbind(sp.table3,(sp.table3[,2]-sp.table3[,1])) %>% 
                    as.data.frame
sp.tableEN15LN$sp_names <- rownames(sp.tableEN15LN)
sp.tableEN15LN$code_change <- rep("EN15.LN",nrow(sp.tableEN15LN))

#x11()
plot(TBI.en15.ln, pch.loss=19, pch.gain=15,
     main="B-C plot",xlim=c(0,1),ylim=c(0,1),diam=TRUE,
     silent = FALSE)


###Plotting all results ####

dat.TBI.ALL <- bind_rows(dat.TBI.SummSpri1,dat.TBI.en9715.1,dat.TBI.en15.ln.1)
rownames(dat.TBI.ALL) <- NULL

dat.TBI.ALL$Code_change <- factor(dat.TBI.ALL$Code_change,
                                     levels=c("SummSpri","EN9715","EN15.LN"))
dat.TBI.ALL$SiteLD <- as.factor(dat.TBI.ALL$SiteLD)

dat.TBI.ALL$latGInt <- factor(dat.TBI.ALL$latGInt,
                              levels = c("5-7","8-10","11-13","14-16"))

write.csv(dat.TBI.ALL,file=paste0(dire,"/New_Data/ENSO_phases/",y.folder,
                                 "/TBIval-ALL.csv"),row.names = FALSE)

change_labs <- c("T1: N_Summer  -  T2: N_Spring",
                "T1: EN_9798  -  T2: EN_1516",
                "T1: EN_1516 - T2: LN_strong")
names(change_labs) <- levels(dat.TBI.ALL$Code_change)

points <- c(22,21)

dat.abline <- data.frame(latGInt=c(rep("5-7",3),rep("8-10",3),rep("11-13",3),
                                   rep("14-16",3)),
                        Code_change=rep(unique(dat.TBI.ALL$Code_change),4),
                         hline=rep(c(TBI.summ.spri$t.test_B.C[1] %>% as.numeric,
                                 TBI.en97.15$t.test_B.C[1] %>% as.numeric,
                                 TBI.en15.ln$t.test_B.C[1] %>% as.numeric),4),
                         slope=rep(1,12))

dat.abline$latGInt <- factor(dat.abline$latGInt,
                             levels = c("5-7","8-10","11-13","14-16"))


tbiplot.all <-
  ggplot(data=dat.TBI.ALL,aes(x=B,y=C,size=D,color=Change)) +
  facet_grid(cols = vars(Code_change),#vars(latGInt),
             labeller = labeller(Code_change=change_labs)) +
  geom_point(alpha=0.65,stroke=1.2) + 
  # geom_text(aes(label=SiteLD),size=2,hjust=0,vjust=0.5,
  #           check_overlap = TRUE) +
  scale_color_manual(values=c("black","red"),name="Species\nchange",
                     breaks = c("+  ","-  "),
                     labels=c("Gain","Loss")) +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(color="green",linetype="longdash",linewidth=0.8) +
  geom_vline(xintercept = 0,color="black") +
  geom_hline(yintercept = 0,color="black") + 
  geom_abline(data=dat.abline, aes(intercept=hline,
                                     slope=slope),
              color="blue",linetype="dashed",linewidth=0.8) +
  labs(x="Species losses (B)",y="Species gain (C)",size="TBI value") +
  theme(panel.background = element_rect(fill = NA,
                                        color="black"),
        legend.position = "top",#legend.box = "vertical",
        plot.title.position = "plot",
        #legend.key.size = unit(3, units = "mm"),
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size=5),
        legend.title = element_text(size=6))


jpeg(filename = paste0(dire,"/",y.folder,"/Plots/",
                       "TBI_AllLat.jpeg"),
     width = 170, height = 100, units = "mm", res = 500)
tbiplot.all
dev.off()


##plot in a map - TBI values #### 

col.fill = c("black","red")

require(mapdata)
require(ggspatial)
peruMap = map_data('worldHires', region = 'Peru')
peruMap <- sf::st_as_sf(map("world",plot=FALSE,fill=TRUE,region="Peru"))

#puertos
ports <- read.csv(paste0(dire,"/puertos1.csv"))

spat.tbi <-
  ggplot() +
  geom_sf(data=peruMap) +
  coord_sf(ylim=c(-16,-4),xlim=c(-83,-74)) +
  #north(peruMap) +
  annotation_scale() +
  scale_y_continuous(breaks=seq(-16,-4,4),
                     labels = c("16°S","12°S","8°S","4°S")) +
  scale_x_continuous(breaks=seq(-83,-74,3),
                     labels = c("83°W","80°W","77°W","74°W")) +
  geom_point(data=dat.TBI.ALL,aes(x=lon,y=lat,color=Change,size=D),
             shape=19,alpha=0.7,stroke=1.5) +
  facet_grid(cols = vars(Code_change),
             labeller = labeller(Code_change=change_labs)) +
  scale_color_manual(values=col.fill,name="Species\nChange",
                     breaks = c("+  ","-  "),
                     labels =c("Gain","Loss")) +
  labs(size="TBI value") +
  guides(fill="none")+
  geom_text(data=ports[c(7,9,12,22),],aes(x=lon+0.7,y=lat,label=Puertos),
            size=2) +
  geom_text(data=ports[26,],aes(x=lon+0.4,y=lat,label=Puertos),
            size=2) +
  geom_text(data=ports[3,],aes(x=lon+0.3,y=lat,label=Puertos),
            size=2) +
  theme(panel.background = element_rect(fill = NA,
                                        color="black"),
        panel.border = element_rect(fill=NA,color="black"),
        axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),legend.position = "top")


spat.tbi

jpeg(filename = paste0(dire,"/",y.folder,"/Plots/",
                       "Spat_TBI_2VAR.jpeg"),
    width = 170, height = 100, units = "mm", res = 500)
spat.tbi
dev.off()


## PLOTTING SPECIES BIOMASS CHANGE FROM T1 TO T2 ####

#loading species SCBD
dat.scbd <- list()
pha <- 1
for (i in seq_along(phase.folder)) {
  beta.div <- load(paste0(dire,"/",y.folder,"/BetaDiv/Anch_",
                          phase.folder[i],".RData")) %>% get()
  scbd.sp <- beta.div$SCBD
  scbd.sp$phase <- rep(phase.folder[i],nrow(scbd.sp))
  dat.scbd[[pha]] <- scbd.sp
  
  pha <- pha + 1
}
dat.scbd.1 <- bind_rows(dat.scbd)

write.csv(dat.scbd.1,file=paste0(dire,"/",y.folder,"/SCBDspALL.csv"))

sp_in <- c(unique(dat.scbd.1$sp),"T.paitensis","S.deliciosa","O.regia") 

#for all phases change
tbi_species <- bind_rows(sp.tableSummSpri,sp.tableEN9715,sp.tableEN15LN)
rownames(tbi_species) <- NULL

tbi_species$code_change <- factor(tbi_species$code_change,
                                  levels=c("SummSpri","EN9715","EN15.LN"))
tbi_species$V12 <- -log1p(tbi_species$V1)

tbi_species$V21 <- log1p(tbi_species$V2)

aba <- tbi_species %>% select(V12,code_change,sp_names)
aba1 <- tbi_species %>% select(V21,code_change,sp_names)
colnames(aba) <- c("V2","code_change","sp_names")
colnames(aba1) <- c("V2","code_change","sp_names")

sp_change <- bind_rows(aba, aba1) %>% filter(sp_names %in% sp_in)
#aba2$compar <- rep("Summ_Nino",nrow(aba2))

sp_change1 <- sp_change %>% filter(sp_names %in% sp_in)
sp.out <- which(sp_change1$sp_names %in% c("Ophichthidae","Nomeidae",
                                 "G.peruvianus","F.corneta"))
sp_change2 <- sp_change1[-sp.out,]

brks <- seq(-15,15,5)
lbls = paste0(as.character(c(seq(15,0,-5), seq(5, 15, 5))))

colors <- c("darkolivegreen","orangered","steelblue")

sp.change <- 
  ggplot(data=sp_change2,aes(x=V2,y=sp_names,fill=code_change)) +
  geom_bar(position="stack",stat="identity",width = .6,alpha=0.7) +
  scale_fill_manual(values=colors, name  ="Phases Comparison",
                    breaks=c("SummSpri", "EN9715","EN15.LN"),
                    labels=c("T1:N_Summer - T2:N_Spring", 
                             "T1:EN9798 - T2:EN1516",
                             "T1:EN1516 - T2:LNstrong")) +
  geom_vline(xintercept = 0,color="black",linewidth=0.5) +
  scale_x_continuous(breaks = brks,   # Breaks
                     labels = lbls) +
  labs(x="Log1p(CPUE)",y="Species",
       title="                 T1                                          T2") +
  theme(panel.grid.minor = element_line(colour = "gray83", 
                                        linewidth = 0.1),
        panel.grid.major = element_line(colour = "gray83", 
                                        linewidth = 0.1),
        panel.background = element_rect(fill = NA,
                                        color="black"),
        axis.title.y = element_blank(),legend.position = c(0.79,0.5),
        legend.key.size = unit(0.3, units = "cm"))


sp.change

jpeg(filename = paste0(dire,"/",y.folder,"/Plots/",
                       "Sp_ChangeTBI1.jpeg"),
     height = 85, width=150, units = "mm", res = 500)
sp.change
dev.off()
