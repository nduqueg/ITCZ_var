rm(list=ls())
cat("\014")

library(reshape)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(hydroTSM)
library(ggplot2)
library(metR)
"%>%" = magrittr::`%>%`

#### DAtes ----
Dates.a <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="month")
Volc <- read.csv("./01_Data/Volcanic_erup.csv")[,-1] %>% 
  melt() %>% 
  dplyr::mutate(Date= paste0(value,"-01-01") %>% as.Date(),
                variable = factor(variable, levels =c("Fischer","VEI5")))
seasons <- c("DJF","MAM","JJA","SON")

years.m <- format(Dates.a, format="%Y"); Years <- unique(format(Dates.a, format="%Y"))
Season.y <- paste0(years.m[-1] ,"-",time2season(Dates.a)[-length(Dates.a)]) %>%  c(.,.[length(.)])
Year.s <- unique(Season.y)

## load data -----
basins <- shapefile("../../01_DataSets/South_America/hybas_sa_lev01-12_v1c/hybas_sa_lev03_v1c.shp")
VIMF.e <- brick("./01_Data/10_VIMF/ModE-Sim_set_1420-1_to_3_ensmean_1420-1849_VIMFeast_seasonal.nc",varname="u") %>% 
  crop(., extent(c(-90,-30,-30,30)))
VIMF.n <- brick("./01_Data/10_VIMF/ModE-Sim_set_1420-1_to_3_ensmean_1420-1849_VIMFnorth_seasonal.nc",varname="v") %>% 
  crop(., extent(c(-90,-30,-30,30)))

cord.vimf <<- rasterToPoints(VIMF.e[[1]])[,c(1,2)]
VIMF.e <- rasterToPoints(VIMF.e)[,-c(1,2)]
VIMF.n <- rasterToPoints(VIMF.n)[,-c(1,2)]

MDiv <- brick("./01_Data/10_VIMF/ModE-Sim_set_1420-1_to_3_ensmean_1420-1849_MDiv_seasonal.nc",varname="vimd") %>%
  shift(.,dx=-180) %>% 
  crop(., extent(c(-90,-30,-30,30)))* 86400 * 92 # transform to mm/season, 92 days between june and august

cord.Mdiv <<- rasterToPoints(MDiv[[1]])[,c(1,2)]
MDiv <- rasterToPoints(MDiv)[,-c(1,2)]

dir.ModRA <- "/mnt/climstor/ERC_PALAEO/ModE-RA/outdata/lowres_20mem_Set_1420-3_1850-1/abs/ensstat/by_var/jja/"
Ppt <- brick(paste0(dir.ModRA,"ModE-RA_lowres_20mem_Set_1420-3_1850-1_ensmean_totprec_abs_1421-2008_jja.nc")) %>% 
  crop(., extent(c(-90,-30,-30,30)))* 86400 * 92 # transform to mm/season, 92 days between june and august

cord.ppt <<- rasterToPoints(Ppt[[1]])[,c(1,2)]
Ppt <- rasterToPoints(Ppt)[,-c(1,2)]

## seasonal Separation  ----
VIMF.e.s <- list(); VIMF.n.s <- list()
MDiv.s <- list()
Ppt.s <- list()

print("Seasonal separation")
for ( i in seasons){
  print(i)
  VIMF.e.s[[i]] <- VIMF.e[, substr(Year.s,6,8)== i]
  VIMF.n.s[[i]] <- VIMF.n[, substr(Year.s,6,8)== i]
  MDiv.s[[i]] <- MDiv[, substr(Year.s,6,8)== i]
  
  # Ppt already selected JJA
  # Ppt.s[[i]] <- Ppt[, substr(Year.s,6,8)== i]
}
Ppt.s[["JJA"]] <- Ppt

## composites ----
VIMF.e.c <- list(); VIMF.n.c <- list()
MDiv.c <- list()
Ppt.c <- list()
l.comp <- c(which(Years == 1783),
            which(Years == 1810),
            which(Years == 1816),
            which(Years == 1823),
            which(Years == 1832),
            which(Years == 1836))
eruptions <- c("Laki","unknown","Tambora","Galanggung","Zavaritskii","Cosigüina")
clim <- which(Years >= 1730 & Years <= 1780)

print("Composites")
VIMF.e.c <- lapply(VIMF.e.s,
                   function(x,t.comp) x[,t.comp] %>% rowMeans(),
                   l.comp) %>% as.data.frame()
VIMF.n.c <- lapply(VIMF.n.s,
                   function(x,t.comp) x[,t.comp] %>% rowMeans(),
                   l.comp) %>% as.data.frame()
MDiv.c <- lapply(MDiv.s,
                 function(x,t.comp) x[,t.comp] %>% rowMeans(),
                 l.comp) %>% as.data.frame()
Ppt.c <- lapply(Ppt.s,
                function(x,t.comp) x[,t.comp] %>% rowMeans(),
                l.comp-1) %>% as.data.frame()

print("Climatology")
VIMF.e.clim <- list(); VIMF.n.clim <- list()
MDiv.clim <- list()
Ppt.clim <- list()
VIMF.e.clim <- lapply(VIMF.e.s,
                   function(x,t.comp) x[,t.comp] %>% rowMeans(),
                   clim) %>% as.data.frame()
VIMF.n.clim <- lapply(VIMF.n.s,
                   function(x,t.comp) x[,t.comp] %>% rowMeans(),
                   clim) %>% as.data.frame()
MDiv.clim <- lapply(MDiv.s,
                 function(x,t.comp) x[,t.comp] %>% rowMeans(),
                 clim) %>% as.data.frame()
Ppt.clim <- lapply(Ppt.s,
                   function(x,t.comp) x[,t.comp] %>% rowMeans(),
                   clim-1) %>% as.data.frame()


print("individual eruptions")
VIMF.e.eru <- list(); VIMF.n.eru <- list()
MDiv.eru <- list()
Ppt.eru <- list()
for( i in 1:length(l.comp)){
  VIMF.e.eru[[i]] <- lapply(VIMF.e.s,
                          function(x,t.comp) x[,t.comp],
                          l.comp[i]) %>% as.data.frame()
  VIMF.n.eru[[i]] <- lapply(VIMF.n.s,
                          function(x,t.comp) x[,t.comp],
                          l.comp[i]) %>% as.data.frame()
  MDiv.eru[[i]] <- lapply(MDiv.s,
                        function(x,t.comp) x[,t.comp],
                        l.comp[i]) %>% as.data.frame()
  Ppt.eru[[i]] <- lapply(Ppt.s,
                         function(x,t.comp) x[,t.comp],
                         l.comp[i]-1) %>% as.data.frame()
}

## differences between the individual eruptions and the climatology ----
VIMF.e.d <- list(); VIMF.n.d <- list()
MDiv.d <- list();
Ppt.d <- list();

for( i in 1:length(l.comp)){
  VIMF.e.d[[i]] <- VIMF.e.eru[[i]] - VIMF.e.clim
  VIMF.n.d[[i]] <- VIMF.n.eru[[i]] - VIMF.n.clim
  MDiv.d[[i]] <- MDiv.eru[[i]] - MDiv.clim
  Ppt.d[[i]] <- Ppt.eru[[i]] - Ppt.clim
}
MDiv.d <- MDiv.d %>% lapply(.,cbind,cord.Mdiv) %>% magrittr::set_names(., eruptions)
Ppt.d <- Ppt.d %>% lapply(.,cbind,cord.ppt) %>% magrittr::set_names(., eruptions)

## angle and magnitud of the differences ----
W.Angle <- list(); W.Mag <- list()

angle.ene <- function(north, east){ atan2(dlat(north), dlon(east, cord.vimf[,2]))*180/pi }
for( i in 1:length(VIMF.e.d)){
  W.Angle[[i]] <- mapply(angle.ene, VIMF.n.d[[i]], VIMF.e.d[[i]]) %>% as.data.frame()
  W.Mag[[i]] <- mapply(metR::Mag, VIMF.n.d[[i]], VIMF.e.d[[i]]) %>% as.data.frame()
}

W.Angle <- lapply(W.Angle, cbind, cord.vimf) %>% magrittr::set_names(., eruptions)
W.Mag <- lapply(W.Mag, cbind, cord.vimf) %>% magrittr::set_names(., eruptions)

### data transformation for plotting individual eruptions ----
data.ang <- melt(W.Angle, id=c("x","y")) %>% magrittr::set_colnames(., c("lon","lat","Season","Angle","Period")) %>% within(., Period <- factor(Period, levels = eruptions))
data.mag <- melt(W.Mag, id=c("x","y")) %>% magrittr::set_colnames(., c("lon","lat","Season","Mag","Period")) %>% within(., Period <- factor(Period, levels = eruptions))
data.g <- merge(data.ang, data.mag, by=c("lon","lat","Season","Period")) %>% 
  magrittr::set_colnames(., c("lon","lat","Season","Period", "Angle","Mag")) %>% within(., Period <- factor(Period, levels = eruptions))

data.mdiv <- melt(MDiv.d, id=c("x","y")) %>% magrittr::set_colnames(., c("lon","lat","Season","MDiv","Period")) %>% 
  within(., {data.class <- MDiv > 0; Period <- factor(Period, levels = eruptions)})
data.Ppt <- melt(Ppt.d, id=c("x","y")) %>% magrittr::set_colnames(., c("lon","lat","Season","Ppt","Period"))  %>% within(., Period <- factor(Period, levels = eruptions))

xlim <- c(-85,-31); ylim <- c(-23,18) 
range <- 200
at.m <- seq(-range,range,length.out = 11) %>% round(.,2); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]

ggplot()+
  facet_wrap(.~ Period, ncol=3)+
  geom_raster(data= data.Ppt %>% subset(., Season=="JJA"),
              aes(x=lon,y=lat,fill=Ppt))+
  scale_fill_stepsn(colours=brewer.pal(11,"BrBG"), breaks=at.m,
                    limits=c(min(at.m),max(at.m)), guide=guide_colorsteps(barwidth=unit(10,"cm")),
                    name="Ppt [mm/season]")+
  
  geom_contour(data=data.mdiv %>% subset(., data.class==T & Season =="JJA"), 
               aes(lon,lat, z=MDiv),color="red",binwidth =200,linewidth=0.35)+ #  linetype=1,
  geom_contour(data=data.mdiv %>% subset(., data.class==F & Season =="JJA"), 
               aes(lon,lat, z=MDiv),color="blue",binwidth =200,linewidth=0.35)+ 
  
  geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  geom_vector(data= data.g %>% subset(., Season =="JJA"), 
              aes(x=lon,y=lat,angle=Angle, mag=Mag), pivot=0.5, skip = 1, col="purple4")+
  scale_mag(name="VIMF", limits=c(0,100))+
  coord_fixed(xlim=xlim,ylim=ylim)+
  labs(x="Longitude [°]",y="Latitude [°]",title= "Year after the volcanic eruption respect the period 1730-1780, JJA")+
  theme_bw()+theme(legend.position = "bottom", strip.text = element_text(size=14))
ggsave("10_VIMF_MDiv_ind_eruptions.png",
       dpi=300,width = 1400*3/300,height = 850*3/300, units = "in")


## differences between the Composite of the eruptions and the climatology ----
VIMF.e.d <- list(); VIMF.n.d <- list()
MDiv.d <- list();
Ppt.d <- list();

VIMF.e.d <- VIMF.e.c - VIMF.e.clim
VIMF.n.d <- VIMF.n.c - VIMF.n.clim
MDiv.d <- MDiv.c - MDiv.clim
Ppt.d <- Ppt.c - Ppt.clim

MDiv.d <- cbind(MDiv.d ,cord.Mdiv)
Ppt.d <- cbind(Ppt.d ,cord.ppt)

## angle and magnitud of the differences ----
W.Angle <- list(); W.Mag <- list()

W.Angle <- mapply(angle.ene,VIMF.n.d, VIMF.e.d) %>% as.data.frame() %>% 
  cbind(., cord.vimf)
W.Mag <- mapply(metR::Mag,VIMF.n.d, VIMF.e.d) %>% as.data.frame() %>% 
  cbind(., cord.vimf)

### data transformation for plotting Composite eruptions ----
data.ang <- melt(W.Angle, id=c("x","y")) %>% magrittr::set_colnames(., c("lon","lat","Season","Angle"))
data.mag <- melt(W.Mag, id=c("x","y")) %>% magrittr::set_colnames(., c("lon","lat","Season","Mag"))
data.g <- merge(data.ang, data.mag, by=c("lon","lat","Season")) %>% 
  magrittr::set_colnames(., c("lon","lat","Season", "Angle","Mag"))

data.mdiv <- melt(MDiv.d, id=c("x","y")) %>% magrittr::set_colnames(., c("lon","lat","Season","MDiv")) %>% 
  within(., data.class <- MDiv > 0)
data.Ppt <- melt(Ppt.d, id=c("x","y")) %>% magrittr::set_colnames(., c("lon","lat","Season","Ppt"))

xlim <- c(-85,-31); ylim <- c(-23,18) 
range <- 200
at.m <- seq(-range,range,length.out = 11) %>% round(.,2); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]

ggplot()+
  geom_raster(data= data.Ppt %>% subset(., Season=="JJA"),
              aes(x=lon,y=lat,fill=Ppt))+
  scale_fill_stepsn(colours=brewer.pal(11,"BrBG"), breaks=at.m,
                    limits=c(min(at.m),max(at.m)), guide=guide_colorsteps(barheight=unit(10,"cm")),
                    name="ModE-RA Ppt\n[mm/season]")+
  
  geom_contour(data=data.mdiv %>% subset(., data.class==T & Season =="JJA"), 
               aes(lon,lat, z=MDiv),color="red",binwidth =200,linewidth=0.35)+ #  linetype=1,
  geom_contour(data=data.mdiv %>% subset(., data.class==F & Season =="JJA"), 
               aes(lon,lat, z=MDiv),color="blue",binwidth =200,linewidth=0.35)+ 
  
  geom_polygon(data=basins,aes(x=long,y=lat, group=group), colour="black",fill="NA",size=0.25)+
  geom_vector(data= data.g %>% subset(., Season =="JJA"), 
              aes(x=lon,y=lat,angle=Angle, mag=Mag), pivot=0.5, skip = 1, col="purple4")+
  scale_mag(name="VIMF", limits=c(0,100))+
  coord_fixed(xlim=xlim,ylim=ylim)+
  labs(x="Longitude [°]",y="Latitude [°]",title= "Composite Year after the six volcanic eruptions, JJA", subtitle="Anomalies respect the period 1730-1780")+
  theme_bw()+
ggsave("11_VIMF_MDiv_Comp_eruptions.png",
       dpi=300,width = 800*3/300,height = 850*3/300, units = "in")
