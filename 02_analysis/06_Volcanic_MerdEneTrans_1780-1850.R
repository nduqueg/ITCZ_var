rm(list=ls())
cat("\014")

library(reshape)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(hydroTSM)
library(ggplot2)
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
Volc <- read.csv("./01_Data/Volcanic_erup.csv")[,-1] %>% 
  melt() %>% 
  within(., variable <- factor(variable, levels =c("Fischer","VEI5"))) %>% 
  magrittr::set_colnames(.,c("variable","Date"))

cdo.cmd <- "cdo seasmean -zonmean ./01_Data/09_VIEnerTrans/ModE-Sim_set_1420-1_to_3_ensmean_VITEFnorth_1420-1849_mon.nc ./01_Data/09_VIEnerTrans/ModE-Sim_set_1420-1_to_3_ensmean_VITEFnorth-ZonMean_1420-1849_seasonal.nc"
system(cdo.cmd)

f <- nc_open("./01_Data/09_VIEnerTrans/ModE-Sim_set_1420-1_to_3_ensmean_VITEFnorth-ZonMean_1420-1849_seasonal.nc")
EneTrans.n <- ncvar_get(f, varid = "dp") %>% 
  .[, seq(from= which(Year.s=="1780-DJF"),
          to= which(Year.s=="1849-SON"))]
lat <- f$dim$lat$vals

Year.s <- Year.s[seq(from= which(Year.s=="1780-DJF"),
                     to= which(Year.s=="1849-SON"))]
Years <- substr(Year.s,1,4) %>% unique()

## seasonal Separation  ----
EneTrans.n.s <- list()

for (i in seasons) EneTrans.n.s[[i]] <- EneTrans.n[, substr(Year.s,6,8)== i] %>% magrittr::set_colnames(., Years) %>% as.data.frame()

## arrange data ----
EneTrans.n.s <- lapply(EneTrans.n.s, cbind, lat)
data.g <- melt(EneTrans.n.s, id="lat") %>% magrittr::set_colnames(., c("lat", "Date","value","Season")) %>% 
  within(.,{
    Date <- as.character(Date) %>% as.numeric()
    value <- value /1e9
  } )

Volc <- Volc %>% subset(., Date >=1780 & Date <=1850)

# plotting world zonal mean ----
range <- max( abs(min(data.g$value)), max(data.g$value))
at.m <- seq(-range,range,length.out = 11) %>% round(.,2); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]

ggplot() + 
  facet_wrap(. ~ Season)+
  geom_tile(data= data.g, aes(Date,lat, fill=value))+
  scale_fill_stepsn(colours=brewer.pal(10,"PRGn"), breaks=at.m,
                    limits=c(min(at.m),max(at.m)))+
  theme_bw()

data.g2 <- data.g %>% subset(.,Season=="JJA")
# range <- max( abs(min(data.g$value)), max(data.g$value))
range <- 3
at.m <- seq(-range,range,length.out = 11) %>% round(.,2); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]
ggplot() + 
  facet_wrap(. ~ Season)+
  geom_tile(data= data.g2, aes(Date,lat, fill=value))+
  geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
  # geom_hline(yintercept = c(0,5,10), linewidth=0.2)+
  scale_fill_stepsn(colours=brewer.pal(10,"PRGn"), breaks=at.m,
                    limits=c(min(at.m),max(at.m)), guide=guide_colorsteps(barwidth=unit(10,"cm")),
                    name="VITEFnorth\n[G W/m]")+
  scale_y_continuous(expand = c(0.01,0.01), breaks = seq(-90,90,by=15))+
  scale_x_continuous(expand = c(0.01,0.01), breaks= seq(1780,1850,by=10))+
  labs(title = "Vertical integral Northward total energy flux, world zonal mean ", y="Latitude [°]")+
  theme_bw()+theme(legend.position = "bottom", 
                   axis.title = element_text(size=12), axis.text = element_text(size=12), strip.text = element_text(size=12),axis.title.x = element_blank(),
                   title = element_text(size=12))


### climatology
EneTrans.n.c <- sapply(EneTrans.n.s, rowMeans) %>% cbind(., lat) %>% as.data.frame()
data.c <- EneTrans.n.c %>% melt(., id="lat") %>% within(., value <- value/1e9)

range <- max( abs(min(data.c$value)), max(data.c$value))
at.m <- seq(-range,range,length.out = 11) %>% round(.,2); at.m.v <- (at.m[-length(at.m)] - at.m[-1])/2 + at.m[-1]

ggplot() + 
  geom_tile(data= data.c, aes(variable, lat, fill=value),width=0.7)+
  geom_hline(yintercept = 0, linewidth=0.2)+
  scale_fill_stepsn(colours=brewer.pal(10,"PRGn"), breaks=at.m,
                    limits=c(min(at.m),max(at.m)), guide=guide_colorsteps(barheight=unit(10,"cm")),
                    name="VITEFnorth\n[G W/m]")+
  scale_y_continuous(expand = c(0.01,0.01), breaks = seq(-90,90,by=15))+
  labs(title = "Vertical integral Northward total energy flux, global zonal mean", 
       subtitle = "ModE-Sim Climatology [1740-1780]", y="Latitude [°]")+
  theme_bw()+theme(axis.title = element_text(size=12), axis.text = element_text(size=12), strip.text = element_text(size=12),axis.title.x = element_blank(),
                   title = element_text(size=12))