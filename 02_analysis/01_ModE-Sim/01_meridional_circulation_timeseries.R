rm(list=ls())
cat("\014")

"%>%"=magrittr::`%>%`
library(ncdf4)
library(abind)
library(hydroTSM)
library(ggplot2)
library(metR)
library(RColorBrewer)

source("00_settings.R")
setwd(dir.base)

### Dates ----
Dates <- list()
Dates$all <- seq(as.Date("1420-01-01"), as.Date("2009-12-01"), by="year")
Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
seasons <- c("DJF","JJA")

## load data ----
mastrf <- list(); mastrf$ep1 <- list(); mastrf$ep2 <- list()

f <- nc_open("./01_Data/01_Streamfunction/ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals

mastrf$ep1$DJF <- ncvar_get(f, varid="mastrfu") %>% .[,, -(length(Dates$ep1)+1)] # exclude the last Dec that belongs to the next year in season DJF
  
f <- nc_open("./01_Data/01_Streamfunction/ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_sJJA.nc")
mastrf$ep1$JJA <- ncvar_get(f, varid="mastrfu")

f <- nc_open("./01_Data/01_Streamfunction/ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_sDJF.nc")
mastrf$ep2$DJF <- ncvar_get(f, varid="mastrfu") %>% .[,, -(length(Dates$ep2)+1)] # exclude the last Dec that belongs to the next year in season DJF

f <- nc_open("./01_Data/01_Streamfunction/ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_sJJA.nc")
mastrf$ep2$JJA <- ncvar_get(f, varid="mastrfu")

################################-
# Mean general circulation ----
################################-

mastrf.longTerm <-list()
for( i in c("ep1","ep2")){ # two time periods
  
  mastrf.longTerm[[i]] <- mastrf[[i]] %>% 
    lapply(., function(x, lat){
      y <- apply(x, c(1,2), mean)  # average time
      colnames(y) <- lev/100
      y <- cbind.data.frame(lat, y)
      return(y)
    },  lat)
}

########   plotting general circulation -

mastrf.longTerm.g <- mastrf.longTerm %>% reshape::melt(., id="lat") %>%
  within(., {
    variable <- as.numeric(as.character(variable))
    value <-  value / 1e9
  })

range <- with(mastrf.longTerm.g, max(abs(min(value)),max(value)) %>% round(., 2))

ggplot() + facet_grid( L1 ~ L2, switch="y" )+
  geom_contour_fill(data= mastrf.longTerm.g, aes(x= lat, y= variable, z= value))+ # ,breaks=MakeBreaks(binwidth=3*10^8, exclude=0)
  scale_fill_stepsn(colours=brewer.pal(11,"RdBu"), breaks=seq(-range,range, length.out=12) %>% round(.,2),limits=c(-range,range), name = "Mass Str F [GKg/s]",
                    guide=guide_colorsteps(even.steps = T,barheight=unit(10,"cm")))+
  # geom_contour(data= mastrf.longTerm.g %>% subset(., value>=0), aes(x= lat, y= variable, z= value), col="blue")+
  # geom_contour(data= mastrf.longTerm.g %>% subset(., value<=0), aes(x= lat, y= variable, z= value), col="red")+
  scale_y_continuous(trans= metR::reverselog_trans())+
  scale_x_continuous(breaks = seq(-90,90,by=15))+
  labs(x="Latitude [°]", y="Level [hPa]", title = "Mass Streamfunction - ModE-Sim ensamble mean")+
  theme_bw()

################################-
# Mid tropospheric calculations & ITCZ features ----
################################-


#######   find mid-tropospheric values

MidT.mastrf <-list()
for( i in c("ep1","ep2")){ # two time periods
  
  MidT.mastrf[[i]] <- mastrf[[i]] %>% 
    lapply(., function(x){
    y <- x[, lev >= 30000 & lev <= 70000,] %>% # select mid-troposphere values
      apply(., c(1,3), mean) %>%  # average the mid-troposphere
      t() %>% zoo(., order.by = Dates[[i]])
    return(y)
  }
  )
}
save(MidT.mastrf, file="./01_Data/01_Streamfunction/ModE-Sim_mastrf_300-700hPa_ensmean_1420-2009.RData")


#######   find features of the ITCZ

Mastrf.prop <- list()
trop.fil <- data.frame(DJF=lat>=-30 & lat <= 30, JJA=lat>=-25 & lat <= 40)

for( i in c("ep1","ep2")){ # two time periods
  for (j in seasons){
    
    x <- MidT.mastrf[[i]][[j]] %>% 
      .[, trop.fil[,j]] # filtering the tropics
    
    Mastrf.prop[[i]][[j]] <- data.frame(Year = index(x),
                                        Max = apply(x, 1, max), lat.Max = apply(x, 1, which.max) %>% lat[ trop.fil[,j] ][.],
                                        Min = apply(x, 1, min), lat.Min = apply(x, 1, which.min) %>% lat[ trop.fil[,j] ][.]
    ) %>% 
      within(.,
             { width <- abs( lat.Max - lat.Min )
             Area <- abs( 2* pi *( 6371e3 )^2* (sin(lat.Max * pi/180) - sin(lat.Min * pi/180)) ) # radious of the earth
             strength <- abs( -9.8*(Max - Min)/Area ) # gravity
             })
  }
}

################################-
## plotting features timeseries ----
################################-

Mastrf.prop.g <- Mastrf.prop %>% reshape::melt(., id=c("Year"))

Mastrf.prop.g %>% subset(.,variable %in% c("Max","lat.Max","Min","lat.Min")) %>%
  ggplot(., aes(x= Year, y= value)) +
  facet_grid(variable ~ L2, switch="y", scales="free_y")+
  geom_line()+ 
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y")+
  labs(title="ModE-Sim Ens. Mean Mass Streamfunction")+
  theme_bw()

Mastrf.prop.g %>% subset(.,variable %in% c("strength","Area","width")) %>%
  ggplot(., aes(x= Year, y= value)) +
  facet_grid(variable ~ L2, switch="y", scales="free_y")+
  geom_line()+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y")+
  labs(title="ModE-Sim Ens. Mean Mass Streamfunction")+
  theme_bw()

