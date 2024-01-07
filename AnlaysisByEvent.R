##For global analysis by event 
#from 1995 to 2019

#Each EN and LN event from 1997 to 2019

rm(list = ls())

require(dplyr)
require(tidyverse)
require(adespatial)
require(vegan)

#working directory
dire <- "~/GitHub/ENSO_analysis"
setwd(dire)

y.folder <- "C_Isopara95-19New"

#reading data
sp_matALL_1 <- read.csv(paste0(dire,"/",y.folder,"/Spcomp_ALL9520.csv"))

#data only from 1995 to 2019
sp_matALL <- sp_matALL_1 %>% filter(Year <= 2019)

sp_matALL$Event <- NA
sp_matALL$Scenario <- NA

##--EL NINO EVENTS
#EN9798 extreme
EN97 <- which(sp_matALL$Code_surv %in% c("C9709-10","C9803-05","C9805-06","C9808-09"))
sp_matALL$Event[EN97] <- rep("EN9798",length(EN97)) 
sp_matALL$Scenario[EN97] <- rep("EN_9798",length(EN97))
#EN2002 weak
EN02 <- which(sp_matALL$Code_surv %in% c("C0202-03","C0209-11"))
sp_matALL$Event[EN02] <- rep("EN02",length(EN02))
sp_matALL$Scenario[EN02] <- rep("EN_weak",length(EN02))
#EN2006 weak
EN06 <- which(sp_matALL$Code_surv %in% c("C0608-09","C0611-12"))
sp_matALL$Event[EN06] <- rep("EN06", length(EN06))
sp_matALL$Scenario[EN06] <- rep("EN_weak",length(EN06))
#EN2008 weak
EN08 <- which(sp_matALL$Code_surv %in% c("C0802-04","C0808-09"))
sp_matALL$Event[EN08] <- rep("EN08", length(EN08))
sp_matALL$Scenario[EN08] <- rep("EN_weak",length(EN08))
#EN2009 weak
EN09 <- which(sp_matALL$Code_surv %in% c("C0908-09"))
sp_matALL$Event[EN09] <- rep("EN09", length(EN09))
sp_matALL$Scenario[EN09] <- rep("EN_weak",length(EN09))
#EN2012 weak
EN12 <- which(sp_matALL$Code_surv %in% c("C1202-04"))
sp_matALL$Event[EN12] <- rep("EN12", length(EN12))
sp_matALL$Scenario[EN12] <- rep("EN_weak",length(EN12))
#EN2014 weak
EN14 <- which(sp_matALL$Code_surv %in% c("C1408-10","C1408-09"))
sp_matALL$Event[EN14] <- rep("EN14", length(EN14))
sp_matALL$Scenario[EN14] <- rep("EN_weak",length(EN14))
#EN1516 strong
EN1516 <- which(sp_matALL$Code_surv %in% c("C1508-10","C1603-04","C1605-06"))
sp_matALL$Event[EN1516] <- rep("EN1516", length(EN1516))
sp_matALL$Scenario[EN1516] <- rep("EN_1516",length(EN1516))
#EN2017 coastal 
EN17 <- which(sp_matALL$Code_surv %in% c("C1703-04"))
sp_matALL$Event[EN17] <- rep("EN17", length(EN17))
sp_matALL$Scenario[EN17] <- rep("EN_coastal",length(EN17))
#EN2018 weak
EN18 <- which(sp_matALL$Code_surv %in% c("C1809-11","C1902-03"))
sp_matALL$Event[EN18] <- rep("EN1819", length(EN18))
sp_matALL$Scenario[EN18] <- rep("EN_weak",length(EN18))
##--LN events
#LN2001 weak
LN01 <- which(sp_matALL$Code_surv %in% c("C0108-09","C0110-11"))
sp_matALL$Event[LN01] <- rep("LN01", length(LN01))
sp_matALL$Scenario[LN01] <- rep("LN_weak",length(LN01))
#LN2007 moderate
LN07 <- which(sp_matALL$Code_surv %in% c("C0702-04","C0708-09"))
sp_matALL$Event[LN07] <- rep("LN07", length(LN07))
sp_matALL$Scenario[LN07] <- rep("LN_strong",length(LN07))
#LN2010 moderate
LN10 <- which(sp_matALL$Code_surv %in% c("C1008-09","C1011-12"))
sp_matALL$Event[LN10] <- rep("LN10", length(LN10))
sp_matALL$Scenario[LN10] <- rep("LN_strong",length(LN10))
#LN2013 strong
LN13 <- which(sp_matALL$Code_surv %in% c("C1302-04","C1308-09"))
sp_matALL$Event[LN13] <- rep("LN13", length(LN13))
sp_matALL$Scenario[LN13] <- rep("LN_strong",length(LN13))
#LN2018 weak
LN18 <- which(sp_matALL$Code_surv %in% c("C1802-04"))
sp_matALL$Event[LN18] <- rep("LN18", length(LN18))
sp_matALL$Scenario[LN18] <- rep("LN_weak",length(LN18))


##Neutral phases -- 
#only adding scenario description
neut.in <- which(is.na(sp_matALL$Event))
#summer months: 11,12,1:4
su.month <- which(sp_matALL$Month_1 %in% c(11:12,1:4))
su.out <- which(neut.in %in% su.month)
sp_matALL$Scenario[neut.in[su.out]] <- rep("N_Summ",length(neut.in[su.out]))
#spring months: 7:10
sp.month <- which(sp_matALL$Month_1 %in% 7:10)
sp.out <- which(neut.in %in% sp.month)
sp_matALL$Scenario[neut.in[sp.out]] <- rep("N_Spri",length(neut.in[sp.out]))

##For plot of alpha diversity 
#only EN and LN events
Eventdat <- sp_matALL[!is.na(sp_matALL$Event),]
q1 <- table(Eventdat$Event) %>% as.data.frame()

#species selected
#66 species selected - high occurrence at least 1 time
sp_sel <- read.csv(paste0(dire,"/table_species1.csv"))

#events
Ev.1 <- Eventdat$Event %>% unique()
nSam <- 50
nboot <- 500
out.dat = NULL
for (i in seq_along(Ev.1)) {
  sel.ev <- Eventdat %>% filter(Event %in% Ev.1[i])
  sp.ev <- sel.ev[,which(colnames(sel.ev)%in% sp_sel$T.st)]
  
  saveAlpha = numeric(nboot)
  saveEven = numeric(nboot)
  for(k in seq_along(1:nboot)){
    sp.ev2 = sp.ev[sample(x = 1:nrow(sp.ev), size = nSam,replace=T), ]
    saveAlpha[k] = mean(specnumber(sp.ev2))
    saveEven[k] = diversity(colMeans(sp.ev2))/log(specnumber(colMeans(sp.ev2)))
    saveBeta[k] = (beta.div(sp.ev2,"hell",nperm=999,adj=TRUE))$beta[2]
  }
  
  sp.ev = data.frame(Event = rep(Ev.1[i],nboot), Alpha = saveAlpha,
                    Even = saveEven, Beta = saveBeta,
                    Scenario = rep(unique(sel.ev$Scenario),
                                                   nboot))
  out.dat = rbind(out.dat, sp.ev)
}


out.dat$Event <- factor(out.dat$Event,levels=Ev.1)
out.dat$Scenario <- factor(out.dat$Scenario, levels=c("EN_9798","EN_1516",
                                                      "EN_weak","EN_coastal",
                                                      "LN_strong","LN_weak"))

col1 <- c("#CC0000","#FF3333","#FF6666","#C1675B","#0033FF","#12BFEE")

event.1 <- ggplot(aes(y = Alpha, x = Event, fill=Scenario), data = out.dat) + 
  geom_boxplot() + 
  scale_fill_manual(values=col1)+
  scale_x_discrete(labels=Ev.1) +
  theme_bw() +
  labs(x = "", y = "Boostrapped Alpha diversity") +
  theme(axis.title = element_text(size=6), axis.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

event.2 <- ggplot(aes(y = Even, x = Event, fill=Scenario), data = out.dat) + 
  geom_boxplot() + 
  scale_fill_manual(values=col1)+
  scale_x_discrete(labels=Ev.1) +
  theme_bw() +
  labs(x = "", y = "Boostrapped Evenness") +
  theme(axis.title = element_text(size=6), axis.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

event.3 <- ggplot(aes(y = Beta, x = Event, fill=Scenario), data = out.dat) + 
  geom_boxplot() + 
  scale_fill_manual(values=col1)+
  scale_x_discrete(labels=Ev.1) +
  theme_bw() +
  labs(x = "", y = "Boostrapped Beta diversity") +
  theme(axis.title = element_text(size=6), axis.text = element_text(size = 5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


jpeg(filename = paste0(dire,"/",y.folder,"/Plots/EventAlphaDiv.jpeg"),
     width = 170, height = 100, units = "mm", res = 500)
event.plot
dev.off()


jpeg(filename = paste0(dire,"/",y.folder,"/Plots/EventDivALL.jpeg"),
     width = 170, height = 170, units = "mm", res = 500)
ggarrange(event.1,event.2,event.3,nrow=3,labels=c("a","b","c"))
dev.off()

jpeg(filename = paste0(dire,"/",y.folder,"/Plots/EventDivALL2.jpeg"),
     width = 170, height = 100, units = "mm", res = 500)
ggarrange(event.1,event.3,nrow=2,labels=c("a","b"),common.legend = TRUE,
          legend="right")
dev.off()



###### Plotting all sites - trawls 
#by ENSO phase

require(mapdata)
require(ggspatial)

#loading plot of segment areas
area_plot <- load(paste0(dire,"/",y.folder,"/SegmPlotVal.RData")) %>% get()

sp_matALL %>% colnames
geog_phases <- sp_matALL %>% select(Year,Lon_1,Lat_1,Month_1,Code_surv,CRUC_LANCE,
                                    Event,Scenario) %>% 
                            filter(Scenario %in% c("N_Summ","N_Spri","EN_9798",
                                                   "EN_1516","EN_weak","LN_strong"))
geog_phases1 <- geog_phases[which(!is.na(geog_phases$Scenario)),] #%>%
                    #filter(Lon_1 <= -74)

geog_phases1$Scenario <- factor(geog_phases1$Scenario,levels=c("N_Summ","N_Spri",
                                               "EN_9798","EN_1516","EN_weak","LN_strong"))

##to remove trawls on the ground 
coast_l <- read.csv("./borde.csv")
in.point <- sp::point.in.polygon(geog_phases1$Lon_1,geog_phases1$Lat_1,
                             coast_l$lon,coast_l$lat)
geog_phases2 <- geog_phases1[which(in.point %in% 1),]

col1 <- c("#336633","#CC9966","#CC0000","#FF3333","#FF6666","#0033FF")

phase.plot <- area_plot +
  geom_point(data=geog_phases2,aes(x=Lon_1,y=Lat_1,color=Scenario,shape=Scenario),
             size=0.3,alpha=0.5) +
  scale_color_manual(values=col1) +
  theme(legend.position = c(0.15,0.3), legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size=5),legend.title = element_text(size=6))

jpeg(filename = paste0(dire,"/",y.folder,"/Plots/Trawl_pointsALL.jpeg"),
     width = 85, height = 100, units = "mm", res = 500)
phase.plot
dev.off()
