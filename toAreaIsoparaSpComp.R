toAreaIsoparaSpComp <- function(data,lon,lat,dire){


#data: spComposition data (only sp_comp)
#lon/la: position of point of survey
#dire: working directory where is isoparalitoral data

##
#isopara: isoparalitoral data
#id: area number
#lon/lat
require(rgeos)
require(dplyr)
require(tidyverse)

#isoparalitoral areas
isopara <- read.csv(paste0(dire,"/areas-isoparalitorales.csv"))
#information of isoparalitoral area
load(paste0(dire,"/infoAreasIsop.RData"))

infoAreasIsop_in <- infoAreasIsop %>% filter(.,lat >=-16 & lat <= -5) %>% 
  filter(dc <= 80)
#intervals for selection by latitude
infoAreasIsop_in$latInt <- infoAreasIsop_in$lat %>% cut(.,seq(-16,-4,1))
#intervals for selection by distance to the coast
infoAreasIsop_in$dcInt <- infoAreasIsop_in$dc %>% cut(.,seq(0,80,10))

isoLat <- infoAreasIsop_in$latInt %>% unique
lat1 <- 5:15
isoDc <- infoAreasIsop_in$dcInt %>% unique
Dc1 <- seq(10,80,10)

#number of isoparalitoral areas
in_1 <- isopara %>% group_by(id,area) %>% count 
in_2 <- which(in_1$area %in% unique(infoAreasIsop_in$code))  
iso_area <- in_1[in_2,1:2] %>% as.data.frame

select <- list()
dat <- list()
indv <- 1

for (i in seq_along(isoLat)) {
  #isoparalittoral polygons by latitude
  iso_select <- which(infoAreasIsop_in$latInt %in% isoLat[i]) 
  
  for (j in seq_along(isoDc)) {
    #isoparalittoral polygons
    #by distance to coast
    sel_1 <- which(infoAreasIsop_in$dcInt[iso_select] %in% isoDc[j])
    iso_zone1 <- infoAreasIsop_in[iso_select[sel_1],1]
    
    #to have lon-lat of the centroid of area
    xym = isopara[which(isopara$area %in% iso_zone1),3:4]
    p = Polygon(xym)
    ps = Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    proj4string(sps) = CRS("+proj=longlat +zone=18 +ellps=WGS84 +datum=WGS84 +no_defs")
    spdf = SpatialPolygonsDataFrame(sps,data.frame(f=99.9))
    
    centroid = gCentroid(spdf,byid=TRUE)
    
    ##survey points inside the polygons
    in_point <- point.in.polygon(lon,lat, #surveys sites
                                 isopara[which(isopara$area %in% iso_zone1),3],
                                 isopara[which(isopara$area %in% iso_zone1),4])
    sel_point <- data[which(in_point %in% c(1,2,3)),]
    
    #to know what sites correspond to each isop area
    in_iso <- which(in_point %in% c(1,2,3))
    code <- paste0(iso_zone1[1],"-",iso_zone1[2])
    dat[[indv]] <- mutate(data[in_iso,],lon[in_iso],lat[in_iso],
                          rep(code,nrow(data[in_iso,])))
    
    #to have the sp comp matrix final with mean of all aggregated values
    mat_species <- apply(sel_point,2,function(x) mean(x,na.rm=T)) %>%
      matrix(.,nrow=1,ncol=ncol(data))
    
    lon_lat <- mutate(centroid %>% as.data.frame(),lat1[i],Dc1[j],
                      isoLat[i],isoDc[j],paste0("S",lat1[i],"-",Dc1[j]))
    lon_lat$n_trawl <- nrow(sel_point)
    lon_lat$shelf <- infoAreasIsop$shelf[which(infoAreasIsop$code %in% 
                                                 iso_zone1)] %>% mean() 
    
    select[[indv]] <- mat_species %>% as.data.frame %>% 
      mutate(.,code,lon_lat) 
    indv <- indv + 1
    
    }
   }

    sp_list <- select %>% bind_rows
    dat_ini <- dat %>% bind_rows()
    colnames(dat_ini) <- c(colnames(data),"lon","lat","IsoCode")

    #eliminate all NA values  
    Sp_Mat <- sp_list[complete.cases(sp_list),]

    #ordering site points for latitude
    sp_comp <- Sp_Mat[,1:ncol(data)] 
    colnames(sp_comp) <- colnames(data)
    #rownames(sp_comp) <- c(1:nrow(sp_comp))

    coord <- Sp_Mat[,((ncol(data)+1):(ncol(data)+3))]
    colnames(coord) <- c("IsoCode","lon","lat")

    geog <- Sp_Mat[,((ncol(data)+4):ncol(Sp_Mat))]
    colnames(geog) <- c("latG","dc","latInt","dcInt","SiteLD","n_trawl","shelf")
    
    OUT_DF = list(sp_comp=sp_comp,coord=coord, geog=geog,
              dat_ini=dat_ini)
    return(OUT_DF)

}
