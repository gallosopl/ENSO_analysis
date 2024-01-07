

toDoSpatialcluster <- function(sp.comp,sp.xy,sp.comp.D,
                                 clust.method){
  
  require(spdep)
  require(codep)
  require(rgdal)
  require(rgeos)
  require(magrittr)
  require(indicspecies)
  require(dplyr)
  require(labdsv)
  require(factoextra)
  require(cluster)
  
  #sp.comp : species composition matrix
  #sp.xy : lon,lat 
  #sp.comp.D : distance matrix
  #clust.method : clustering method ("ward.D","ward.D2","complete",etc.)
  
  
  #generating data Spatial points
  coords.dat <- SpatialPointsDataFrame(coords=sp.xy,data=as.data.frame(sp.comp),
                                   proj4string=CRS("+proj=longlat +zone=18 +ellps=WGS84 +datum=WGS84 +no_defs"))
  #List of neighbours
  listW <- nb2listw(tri2nb(as.matrix(sp.xy)), style = "B")
  neighbors <- listW %>% listw2sn %>% {.[,1:2]}
  neighbors %<>% {.[.[,1L]<.[,2L],]}
  
  #To know distance between points
  spTa  <- SpatialPoints(data.frame(sp.xy))
  proj4string(spTa) <- CRS("+proj=longlat")
  spTa.proj1  <- spTransform(spTa, CRS("+proj=utm +zone=18 ellips=WGS84"))
  dist <- gDistance(spTa.proj1,byid=T) #distance 
  nbdst <- numeric(nrow(neighbors))
  for(l in 1L:nrow(neighbors)) {
    ## i=1L
    nbdst[l] <- dist%>%as.matrix%>%{.[neighbors[l,1],neighbors[l,2]]}
  }
  #to select distances more closer
  threshold <- quantile(nbdst,prob=0.95)
  closeNeighbors <- neighbors[nbdst<threshold,]
  
  #######
  #### To start with the spatial clustering
  cclust <- constr.hclust(sp.comp.D, method=clust.method, 
                          links=closeNeighbors[complete.cases(closeNeighbors),], 
                          coords.dat@coords)
  
  hc <- cclust
  #to optain optimal number of cluster
  # Average silhouette widths (Rousseeuw quality index)
  Si <- numeric(nrow(sp.comp))
  for (k in 2:(nrow(sp.comp) - 1))
  { sil <- silhouette(cutree(hc, k = k), sp.comp.D)
    Si[k] <- summary(sil)$avg.width
  }
  #K OPTIMAL NUMBER
  k.best1 <- which.max(Si)


  #optimal number of clusters by IndVal
  IndVal <- numeric(nrow(sp.comp))
  ng <- numeric(nrow(sp.comp))
  for (k in 2:(nrow(sp.comp) - 1))
  {
    iva <- indval(sp.comp, cutree(hc, k = k), numitr = 10000)
    gr <- factor(iva$maxcls[iva$pval <= 0.05])
    ng[k] <- length(levels(gr)) / k
    iv <- iva$indcls[iva$pval <= 0.05]
    IndVal[k] <- sum(iv)
  }
  k.best2 <- which.max(IndVal[ng == 1]) + 1
  col3 <- rep(1, nrow(sp.comp))
  col3[ng == 1] <- 3

  
  OUT_DF <- list(clust=cclust,k=c(k_sil=k.best1,k_IndVal=k.best2))
  return(OUT_DF)
  
  
}