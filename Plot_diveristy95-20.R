##Plotting diversity values

#working directory
dire <- "~/GitHub/ENSO_analysis"
setwd(dire)

y.folder <- "C_Isopara95-19New"

require(tidyverse)
require(ggplot2)
require(vegan)
require(adespatial)

phasesDat <- list.files(path=paste0(dire,"/",y.folder,
                                    "/Dat_surveys"))
#67 species selected - high occurrence at least 1 time
sp_sel <- read.csv(paste0(dire,"/table_species1.csv"))

nSam <- 150
nboot <- 500
out2 = NULL
for (i in seq_along(phasesDat)) {
  
  mat_phase <- read.csv(file=paste0(dire,"/",y.folder,
                                    "/Dat_surveys/",phasesDat[i]))
  
  phase.folder <- strsplit(phasesDat[i],".csv",fixed=TRUE) %>% as.character()
  
  #no considering anchovy
  sp_sel1 <- sp_sel[-which(sp_sel$T.st %in% "E.ringens"),5]
  
  
  tmp <- mat_phase[,which(colnames(mat_phase)%in% sp_sel$T.st)]
  
  saveAlpha = numeric(nboot)
  saveEven = numeric(nboot)
  saveBeta = numeric(nboot)
  for(k in seq_along(1:nboot)){
    tmp2 = tmp[sample(x = 1:nrow(tmp), size = nSam,replace=T), ]
    saveAlpha[k] = mean(specnumber(tmp2))
    saveEven[k] = diversity(colMeans(tmp2))/log(specnumber(colMeans(tmp2)))
    saveBeta[k] = (beta.div(tmp2,"hell",nperm=999,adj=TRUE))$beta[2]
  }
  
  tmp3 = data.frame(Scenario = rep(phase.folder,nboot), Alpha = saveAlpha,
                    Even = saveEven, Beta = saveBeta)
  out2 = rbind(out2, tmp3)
  
}

out2$Scenario <- factor(out2$Scenario,levels=c("summ_surv","spr_surv",
                                               "en_9798_surv","en_1516_surv",
                                               "wk.en_surv","st.ln_surv"))

col1 <- c("#336633","#CC9966","#CC0000","#FF3333","#FF6666","#0033FF")

a1 <- ggplot(aes(y = Alpha, x = Scenario), data = out2) + 
  geom_boxplot(fill=col1,col=1) + 
  scale_x_discrete(labels=c("N_Summer","N_Spring","EN_9798","EN_1516",
                              "EN_weak","LN_strong")) +
  theme_bw() +
  labs(x = "", y = "Boostrapped Alpha diversity") +
  theme(axis.title = element_text(size=6), axis.text = element_text(size = 5))#,
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

a2<- ggplot(aes(y = Even, x = Scenario), data = out2) + 
  geom_boxplot(fill=col1,col=1) + 
  scale_x_discrete(labels=c("N_Summer","N_Spring","EN_9798","EN_1516",
                            "EN_weak","LN_strong")) +
  theme_bw() +
  labs(x = "", y = "Boostrapped Evenness") +
  theme(axis.title = element_text(size=6), axis.text = element_text(size =5))#,
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


a3 <- ggplot(aes(y = Beta, x = Scenario), data = out2) + 
  geom_boxplot(fill=col1,col=1) + 
  scale_x_discrete(labels=c("N_Summer","N_Spring","EN_9798","EN_1516",
                            "EN_weak","LN_strong")) +
  theme_bw() +
  labs(x = "", y = "Boostrapped Beta diversity") +
  theme(axis.title = element_text(size=6), axis.text = element_text(size = 5))#,
#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


require(ggpubr)


jpeg(filename = paste0(dire,"/",y.folder,"/Plots/","AlphaDiv_all.jpeg"),
     width = 85, height = 170, units = "mm", res = 500)
a1
dev.off()


jpeg(filename = paste0(dire,"/",y.folder,"/Plots/","ALLDiv_boots.jpeg"),
     width = 85, height = 170, units = "mm", res = 500)
ggarrange(a1,a2,a3,nrow=3,labels=c("a","b","c"))
dev.off()

jpeg(filename = paste0(dire,"/",y.folder,"/Plots/","ALLDiv_boots2.jpeg"),
     width = 85, height = 100, units = "mm", res = 500)
ggarrange(a1,a3,nrow=2,labels=c("a","b"))
dev.off()



