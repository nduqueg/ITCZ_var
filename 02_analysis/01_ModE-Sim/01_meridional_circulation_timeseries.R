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
# Mid tropospheric calculations ----
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
for( i in c("ep1","ep2")){ # two time periods
  
  Mastrf.prop[[i]] <- MidT.mastrf[[i]] %>% 
    
    lapply(., function(x){ # over each season ( DJF & JJA )
      
      max.min.df <- data.frame(Year = index(x),
                               Max = apply(x, 1, max), Max.lat = apply(x, 1, which.max) %>% lat[.],
                               Min = apply(x, 1, min), Min.lat = apply(x, 1, which.min) %>% lat[.]
                               ) %>% 
        within(.,
               { width <- abs( Max.lat - Min.lat )
               Area <- abs( 2* pi *( 6371*10^3 )^2* (sin(Max.lat * pi/180) - sin(Min.lat * pi/180)) ) # radious of the earth
               strength <- abs( -9.8*(Max - Min)/Area ) # gravity
               })
      
    })
}

################################-
## plotting features timeseries ----
################################-

Mastrf.prop.g <- Mastrf.prop %>% reshape::melt(., id=c("Year"))

Mastrf.prop.g %>% subset(.,variable %in% c("Max","Max.lat","Min","Min.lat")) %>%
  ggplot(., aes(x= Year, y= value)) +
  facet_grid(variable ~ L2, switch="y", scales="free_y")+
  geom_line()+ 
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y")+
  theme_bw()

Mastrf.prop.g %>% subset(.,variable %in% c("width","Area","strength")) %>%
  ggplot(., aes(x= Year, y= value)) +
  facet_grid(variable ~ L2, switch="y", scales="free_y")+
  geom_line()+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y")+
  theme_bw()

