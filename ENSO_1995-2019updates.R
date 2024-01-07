#### 
#Updating analysis for ENSO_phases
#Considering surveys: 1995-2019

#Summer: months 11-5
#Spring: months 7-10
#Types of ENSO: Extreme EN 97-98,Strong EN 05-06, 
#               Moderate EN, Strong/Moderate LN


rm(list = ls())

require(dplyr)
require(tidyverse)
require(adespatial)
require(vegan)

#working directory
dire <- "~/GitHub/ENSO_analysis"
setwd(dire)

y.folder <- "C_Isopara95-19New"

#loading data of ALL surveys
#this data has ONLY FISHES
#has MPH surveys
#from 1983 to 2020

surv_pelag1 <- read.csv(paste0(dire,"/",y.folder,
                               "/baseFISHES_1983-2020changeLAST.csv"))
# match(surv_pelag$CRUC_LANCE,a1$CRUC_LANCE)
# surv_pelag$Lat_1 <- as.numeric(surv_pelag$Lat_1)
# surv_pelag$Lon_1 <- as.numeric(surv_pelag$Lon_1)
# surv_pelag$Year <- as.numeric(surv_pelag$Year)

surv_pelag <- surv_pelag1 %>% filter(Year >= 1995) %>% #from 1995
                filter(Lat_1 >=-16 & Lat_1 <= -5) #spatial area

#matrix of sp composition #FOR ALL DATA
surv_pelag$CAPTSP_STAND  %>% class

aba <- surv_pelag %>% group_by(CRUC_LANCE,Label) %>% summarise(mean(CAPTSP_STAND))

aba %>% pivot_wider(id_cols ="Label",names_from="CRUC_LANCE",
                    values_from = "mean(CAPTSP_STAND)")

mat1 <- surv_pelag %>% pivot_wider(id_cols = "Label",
                                    names_from = "CRUC_LANCE",
                                    values_from = "CAPTSP_STAND",
                                   values_fn = mean) 
mat2 <- mat1[,-1] %>% as.matrix %>% as.numeric %>%
  matrix(.,ncol = ncol(mat1[,-1]))
colnames(mat2) <- mat1[,-1] %>% colnames(.) 
row.names(mat2) <- mat1$Label
#transpuesta de mat2
sp_mat <- t(mat2) %>% replace_na(.,0)
sp_mat[is.na(sp_mat)] <- 0
# sp_mat1 <- sp_mat
# for (i in 1:ncol(sp_mat1)) {
#   sp_mat1[,i] <- replace_na(sp_mat1[,i],0)
# }
##Adding lon & lat to sp matrix
in_trawl <- rownames(sp_mat)
lon_lat.in <- list()
sel <- 1
for (i in seq_along(in_trawl)) {
  ab <- surv_pelag %>% filter(CRUC_LANCE %in% in_trawl[i]) 
  lon_lat.in[[sel]] = ab[1,] 
  sel <- sel + 1
}
lon_lat <- lon_lat.in %>% bind_rows %>% select(Year,Lat_1:Month_1,Code_surv,
                                               CRUC_LANCE)
#SP COMPOSITION DATA WITH LON-LAT, YEAR, CODE SURVEY
#This is the full data for surveys 1995-2020
sp_matALL <- mutate(lon_lat,sp_mat %>% as.data.frame)

# write.csv(sp_matALL,file=paste0(dire,"/y.folder/Spcomp_ALL9520.csv"),
#           row.names = FALSE)


############
### STARTING TO SELECT SURVEYS BY ENSO PHASE ####

#working directory
dire <- "~/GitHub/ENSO_analysis"
setwd(dire)

y.folder <- "C_Isopara95-19New"

#reading data
sp_matALL_1 <- read.csv(paste0(dire,"/",y.folder,"/Spcomp_ALL9520.csv"))

#data only from 1995 to 2019
sp_matALL <- sp_matALL_1 %>% filter(Year <= 2019)

#selecting surveys by scenario:
#Extreme EN 1997-1998:
en_9798.in <- c("C9709-10","C9803-05","C9805-06","C9808-09")
en_9798_surv <- sp_matALL %>% filter(Code_surv %in% en_9798.in)
en_9798_surv$Code_surv %>% unique
en_9798_surv[,3:2] %>% plot

en_1 <- en_9798_surv %>% select(Year,Month_1,Code_surv,CRUC_LANCE) %>%
  group_by(Year,Month_1,Code_surv) %>% count(Code_surv) 
en_1$Condition <- rep("EN_9798",nrow(en_1))

# write.csv(en_9798_surv,file=paste0(dire,"/",y.folder,
#                                    "/Dat_surveys/en_9798_surv.csv"),row.names = FALSE)

#Strong EN 2015-2016:
en_1516.in <- c("C1508-10","C1603-04","C1605-06")
en_1516_surv <- sp_matALL %>% filter(Code_surv %in% en_1516.in)
en_1516_surv$Code_surv %>% unique
en_1516_surv[,3:2] %>% plot

en_2 <- en_1516_surv %>% select(Year,Month_1,Code_surv,CRUC_LANCE) %>%
  group_by(Year,Month_1,Code_surv) %>% count(Code_surv) 
en_2$Condition <- rep("EN_1516",nrow(en_2))


# write.csv(en_1516_surv,file=paste0(dire,"/",y.folder,
#                                  "/Dat_surveys/en_1516_surv.csv"),row.names = FALSE)

#Strong-moderate LN:
st.ln_in <- c("C0702-04","C0708-09","C1008-09",
              "C1011-12","C1302-04","C1308-09")
st.ln_surv <- sp_matALL %>% filter(Code_surv %in% st.ln_in)
st.ln_surv$Code_surv %>% unique
st.ln_surv[,3:2] %>% plot

ln_1 <- st.ln_surv %>% select(Year,Month_1,Code_surv,CRUC_LANCE) %>%
  group_by(Year,Month_1,Code_surv) %>% count(Code_surv) 
ln_1$Condition <- rep("LN_SM",nrow(ln_1))

# write.csv(st.ln_surv,file=paste0(dire,"/",y.folder,
#                                  "/Dat_surveys/st.ln_surv.csv"),row.names = FALSE)

#Weak/modarate EN:
wk.en_in <- c("C0202-03","C0209-11","C0608-09","C0611-12","C0802-04",
              "C0808-09","C0908-09","C0908-09","C1202-04","C1408-10",
              "C1408-09","C1809-11","C1902-03")
wk.en_surv <- sp_matALL %>% filter(Code_surv %in% wk.en_in)
wk.en_surv$Code_surv %>% unique
wk.en_surv[,3:2] %>% plot

en_wk <- wk.en_surv  %>% select(Year,Month_1,Code_surv,CRUC_LANCE) %>%
  group_by(Year,Month_1,Code_surv) %>% count(Code_surv)
en_wk$Condition <- rep("EN_Weak",nrow(en_wk))

# write.csv(wk.en_surv,file=paste0(dire,"/",y.folder,
#                                  "/Dat_surveys/wk.en_surv.csv"),row.names = FALSE)

#Surveys to not consider: EN moderate (coastal EN) and LN weak
surv.out <- c("C1703-04","C0108-09","C0110-11","C1802-04")


#Neutral surveys
neut_out <- which(sp_matALL$Code_surv %in% c(en_9798.in,en_1516.in,st.ln_in,
                                             wk.en_in,surv.out))
neut_surv <- sp_matALL[-neut_out,]
neut_surv$Code_surv %>% unique

(neut_surv %>% filter(Month_1 %in% 5:7))$Code_surv %>% unique

#Summer neutral
summ_surv <- neut_surv %>% filter(Month_1 %in% c(1:4,11,12))
summ_surv$Code_surv %>% unique 
summ_surv[,3:2] %>% plot

summ_1 <- summ_surv %>% select(Year,Month_1,Code_surv,CRUC_LANCE) %>%
  group_by(Year,Month_1,Code_surv) %>% count(Code_surv) 
summ_1$Condition <- rep("N_Summ",nrow(summ_1))

# write.csv(summ_surv,file=paste0(dire,"/",y.folder,
#                                  "/Dat_surveys/summ_surv.csv"),row.names = FALSE)

#Spring neutral
spr_surv <- neut_surv %>% filter(Month_1 %in% 7:10)
spr_surv$Code_surv %>% unique 
spr_surv[,3:2] %>% plot

spr_1 <- spr_surv %>% select(Year,Month_1,Code_surv,CRUC_LANCE) %>%
  group_by(Year,Month_1,Code_surv) %>% count(Code_surv) 
spr_1$Condition <- rep("N_Spri",nrow(spr_1))

# write.csv(spr_surv,file=paste0(dire,"/",y.folder,
#                                  "/Dat_surveys/spr_surv.csv"),row.names = FALSE)

table_surv <- bind_rows(en_1,en_2,ln_1,en_wk,summ_1,spr_1)
# write.csv(table_surv,file=paste0(dire,"/",y.folder,
#                                "/table_surv.csv"),row.names = FALSE)


##############################
####
###Arrange data to isoparalittoral areas##########
#loading function
source(paste0(dire,"/toAreaIsoparaSpComp.R")) #for isoparalitoras areas
source(paste0(dire,"/toDoSpatialcluster.R")) #for spatial clustering

require(vegan)
require(constr.hclust)
require(factoextra)
require(adespatial)

phasesDat <- list.files(path=paste0(dire,"/",y.folder,
                                    "/Dat_surveys"))

#67 species selected - high occurrence at least 5 times
sp_sel <- read.csv(paste0(dire,"/table_species1.csv"))

#### For spatial cluster of all data ####
##considering first agrregated data##
for (i in seq_along(phasesDat)) {
  mat_phase <- read.csv(file=paste0(dire,"/",y.folder,
                                 "/Dat_surveys/",phasesDat[i]))
  
  phase.folder <- strsplit(phasesDat[i],".csv",fixed=TRUE) %>% as.character()
  
  Mat_isop <- toAreaIsoparaSpComp(data=mat_phase[,-(1:6)],lon=mat_phase$Lon_1,
                                  lat=mat_phase$Lat_1,dire=dire)
  
  ##segment data
  save(Mat_isop,file=paste0(dire,"/",y.folder,
                               "/Segm_",phase.folder,".RData"))
  
  ##community data with anchovy - abundance##
  spe.1Anch <- Mat_isop$sp_comp[,which(colnames(Mat_isop$sp_comp) %in%
                                         sp_sel$T.st)] 
  spe.1Anch1 <- spe.1Anch[,which(colSums(spe.1Anch)>0)] #species mat
  spe.1Anch.dc <- vegdist(decostand(spe.1Anch1,"hell"),"euc") #distance mat
  
  ##community data without anchovy - abundance##
  sp_sel1 <- sp_sel[-which(sp_sel$T.st %in% c("E.ringens")),5]
  spe.1NoAnch <- Mat_isop$sp_comp[,which(colnames(Mat_isop$sp_comp)
                                   %in% sp_sel1)]
  spe.1NoAnch1 <- spe.1NoAnch[,which(colSums(spe.1NoAnch)>0)] #species mat
  spe.1NoAnch.dc <- vegdist(decostand(spe.1NoAnch1,"hell"),"euc") #distance mat 
  
  ##community data presence-absence - including anchovy##
  spe.paAnch <- decostand(spe.1Anch1,"pa") #species mat
  spe.paAnch.dc <- dist(spe.1Anch1,"binary") #distance mat
  
  ##all data in one:
  spe <- list(spe.1Anch1,spe.1NoAnch1,spe.paAnch) #species mat
  spe.dc <- list(spe.1Anch.dc,spe.1NoAnch.dc,spe.paAnch.dc) #distance mat
  sp.xy <- Mat_isop$coord #lonlat
  
  #saving data of species matrix with anchovy condition
  save(spe,file=paste0(dire,"/",
                            y.folder,"/",phase.folder,"_MatDat.RData"))
  
  name <- c("Anch","NoAnch","paAnch")
  
  ###Spatial cluster###
  
  for (j in 1:3) {
    clust_phaseAnch <- toDoSpatialcluster(spe[[j]],sp.xy[,2:3],
                                          spe.dc[[j]],"ward.D2")
    save(clust_phaseAnch,file=paste0(dire,"/",
                                     y.folder,"/SpatialClust/",name[j],"_",
                                     phase.folder,".RData"))
  }
  
  ###Beta diversity###
  method.beta <- c("hell","hell","jaccard")
  
  for (h in 1:3) {
    #abundance with anchovy --
    sp.beta.div <- beta.div((spe[[h]]),method.beta[[h]], nperm=99999, 
                            adj=TRUE)
    #BD total
    sp.beta.div$beta[2]
    #species with SCBD larger than the mean SCBD
    sp.scbdIn <- which(sp.beta.div$SCBD > mean(sp.beta.div$SCBD))
    
    #species mean biomass and frequency
    sp.meanIN <- spe[[h]] %>% select(names(sp.scbdIn)) 
    mean.spBiomass <- apply(sp.meanIN,2,mean)
    freq.sp <- sp.meanIN %>% decostand(.,"pa") %>% apply(.,2,sum)
    
    SCBD.sp <- data.frame(sp=names(sp.meanIN), 
                          SCBD = sp.beta.div$SCBD[sp.scbdIn] %>% as.numeric,
                          meanBiom=mean.spBiomass %>% as.numeric(),
                          freq = freq.sp %>% as.numeric)
    
    #Holm correction of LCBD p-values 
    p.adj1 <- p.adjust(sp.beta.div$p.LCBD,method = "holm")
    #Sites with significant LCBD value after correction
    sit <- which(p.adj1<=0.05) %>% as.numeric
    #sites significant with no correction
    sit2 <- which(sp.beta.div$p.LCBD <= 0.05) %>% as.numeric() 
    
    LCBD.sig <- data.frame(site=sit, id_area = sp.xy[sit,1],
                           LCBD=sp.beta.div$LCBD[sit] %>% as.numeric,
                           p.valNO=sp.beta.div$p.LCBD[sit] %>% as.numeric,
                           p.valYES=p.adj1[sit] %>% as.numeric,
                           spe[[h]][sit,])
    
    
    beta.div <- list(beta=sp.beta.div$beta,
                     SCBD=SCBD.sp,
                     LCBD=sp.beta.div$LCBD,
                     LCBD.sig = LCBD.sig,
                     pLCBD1=p.adj1, pLCBD=sp.beta.div$p.LCBD,
                     Rich=specnumber(spe[[h]]))
    
    save(beta.div,file=paste0(dire,"/",
                              y.folder,"/BetaDiv/",name[h],"_",
                              phase.folder,".RData"))  
    
  }
  
  print(phase.folder)
  
}


##Plotting Spatial clustering
require(mapdata)

#ordering phases
phasesDat <- c("summ_surv","spr_surv","en_9798_surv","en_1516_surv",
               "wk.en_surv","st.ln_surv")

name.phase <- c("N_Summer","N_Spring","EN_9798","EN_1516",
                "EN_weak","LN_strong")

ports <- read.csv(paste0(dire,"/puertos1.csv"))

colors <- rainbow(8)
jpeg(filename = paste0(dire,"/",y.folder,"/Plots/",
                      "spat_cluster_ALL2.jpeg"),
    width = 170, height = 150, units = "mm", res = 500)
par(mar=c(1,3,1,0.5),mfrow=c(2,3))
for (l in seq_along(phasesDat)) {
  
  #1995-2019
  lab <- strsplit(phasesDat[l],".csv",fixed=TRUE) %>% as.character()
  
  #loading cluster data
  cclust.new <- load(paste0(dire,"/",y.folder,"/SpatialClust/Anch_",lab,
                            ".RData")) %>% get()
  if (is.na(cclust.new$k[2]) == TRUE)
    k = 2
  else k = cclust.new$k[2]
  
  constr.hclust:::plot.constr.hclust(cclust.new$clust, k=k, links=T, xlab="", 
                                     ylab="Latitude",cex=0.1,lwd=2,pch=0,
                                     col=colors[1:k],
                                     xlim=c(-83,-74), axes=FALSE)
  axis(2,las=1,seq(-16,-6,2),labels = c("16°S","14°S","12°S","10°S","8°S","6°S"))
  map("worldHires", add=T, fill=T, col=8)
  legend("bottomleft", lty=1L, lwd=3, col=colors[1:k],
         legend=sprintf("Group %d",1:k), cex=0.9, bty="n")
  legend(-79.25,-4.5,legend=name.phase[l],bty="n",cex=1.2)
  text(ports[c(7,9,12,22),2]+0.85,ports[c(7,9,12,22),3],labels=ports[c(7,9,12,22),1],cex=0.75)
  text(ports[26,2]+0.5,ports[26,3],labels=ports[26,1],cex=0.75)
  text(ports[3,2]+0.5,ports[3,3],labels=ports[3,1],cex=0.75)
  box()
}
dev.off()


# ### Plotting the LCBD values 
# 
# 
# #1995-2019
# lab <- strsplit(phasesDat[l],".csv",fixed=TRUE) %>% as.character()
# 
# #loading cluster data
# cclust.new <- load(paste0(dire,"/",y.folder,"/Spat_clust_",
#                            lab,".RData")) %>% get()
# 
# betadiv.new <- load(paste0(dire,"/",y.folder,"/beta_div_",
#                           lab,".RData")) %>% get()
# 
# 
# 
# 
# ##To do the boostrap of sites
# ##Spatial cluster - beta diversity ##
# for (i in seq_along(phasesDat)) {
#   
#   sp_mat <- read.csv(file=paste0(dire,"/",y.folder,
#                                  "/Dat_surveys/",phasesDat[i]))
#   
#   phase.folder <- strsplit(phasesDat[i],".csv",fixed=TRUE) %>% as.character()
#   
#   Mat_isop <- toAreaIsoparaSpComp(data=sp_mat[,-(1:6)],lon=sp_mat$Lon_1,
#                                   lat=sp_mat$Lat_1,dire=dire)
#   
#   #number of trawls selected by isoparalitoral area
#   #sample <- quantile(Mat_isop$geog$n_trawl,.25) %>% round()
#   sample <- min(Mat_isop$geog$n_trawl)
#   
#   nboot <- 10 #boostrapped number 
#   
#   for (xy in 1:nboot) {
#     code.iso.area <- Mat_isop$dat_ini$IsoCode %>% unique
#     
#     dat1 <- list()
#     in.1 <- 1
#     for (j in seq_along(code.iso.area)) {
#         a1 <- which(Mat_isop$dat_ini$IsoCode %in% code.iso.area[j])
#         sel<- sample(x=a1,size=sample,replace=FALSE)
#         dat <- Mat_isop$dat_ini[sel,1:ncol(Mat_isop$sp_comp)]
#         dat_f <- apply(dat,2,function(x) mean(x,na.rm=T)) %>%
#                     matrix(.,nrow=1,ncol=ncol(dat)) %>% as.data.frame() 
#         dat_f$IsoCode <- code.iso.area[j]
#         dat1[[in.1]] <- dat_f
#         in.1 <- in.1 + 1
#       }
#     spmat1 <- dat1 %>% bind_rows()
#     colnames(spmat1) <- c(colnames(Mat_isop$sp_comp),"IsoCode")
#     #species matrix - selected species
#     sp_matf <- spmat1[,which(colnames(spmat1) %in% sp_sel$x)] 
#     out.anch<- which((sp_matf %>% colnames) %in% "E.ringens") 
#     sp_matf1 <- sp_matf[,-out.anch]
#     #geographic data
#     in.match <- match(spmat1[,ncol(spmat1)],Mat_isop$coord$IsoCode)
#     geog_matf <- Mat_isop$coord[in.match,] %>% mutate(.,Mat_isop$geog[in.match,]) 
#     
#     ##Starting with spatial clustering
#     sp.comp <- sp_matf1[,which(colSums(sp_matf1)>0)]
#     sp.xy <- geog_matf[,2:3]
#     
#     #distance matrix 
#     sp.comp.D <- vegdist(decostand(sp.comp,"hel"),"euc")     
#     ###to do spatial cluster###
#     cclust <- toDoSpatialcluster(sp.comp,sp.xy,sp.comp.D,"ward.D2")
#     
#     #plotting
#     par(mfrow=c(1,2),mar=c(2,2,2,0))
#     constr.hclust:::plot.constr.hclust(cclust$clust, k=cclust$k[1], links=T,
#                                        main=paste(phase.folder,xy))
#     constr.hclust:::plot.constr.hclust(cclust$clust, k=cclust$k[2], links=T)
#     
#     
#     ###Beta diversity###
#     sp.beta.div <- beta.div((sp.comp),"hell", nperm=9999)
#     
#     #data of beta div results - all
#     beta.div <- list(beta=sp.beta.div$beta,
#                      SCBD=which(sp.beta.div$SCBD > mean(sp.beta.div$SCBD)),
#                      LCBD=sp.beta.div$LCBD,
#                      pLCBD1=p.adj1, pLCBD=sp.beta.div$p.LCBD,
#                      Rich=specnumber(sp.comp))
#     
#     
#     ##Results of both spatial clustering and beta diversity values
#     clustPhase <- list(cclust,beta.div)
#     
#     ##saving ALL results by boostrapp to its ENSO phase folder
#     save(clustPhase,file=paste0(dire,"/",y.folder,"/Phases/",
#                                 phase.folder,"/clustPhase",xy,".RData"))  
#     
#   }
#   
# }
# 
# i=4
# par(mfrow=c(2,2),mar=c(1,2,2,0))
# for (xy in 1:4) {
# 
# phase.folder <- strsplit(phasesDat[i],".csv",fixed=TRUE) %>% as.character()
# 
# cclust1 <- load(paste0(dire,"/",y.folder,"/Phases/",
#              phase.folder,"/clustPhase",xy,".RData")) %>% get
# 
# constr.hclust:::plot.constr.hclust(cclust1[[1]]$clust, k=3, links=T,
#                                    main=paste(phase.folder,xy))
# }
