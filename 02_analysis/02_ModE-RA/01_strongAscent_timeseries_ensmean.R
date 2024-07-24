rm(list=ls())
cat("\014")

"%>%"=magrittr::`%>%`
library(ncdf4)
library(abind)
library(hydroTSM)
library(ggplot2)
library(RColorBrewer)

source("00_settings.R")
setwd(dir.base)

### Dates ----
Dates <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")

seasons <- c("DJF","JJA")

## load data ----
Omega <- list()

f <- nc_open("./01_Data/02_Omega500/ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals

Omega$DJF <- ncvar_get(f, varid="omega") %>% .[, -(length(Dates)+1)] # exclude the last Dec that belongs to the next year in season DJF

f <- nc_open("./01_Data/02_Omega500/ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_sJJA.nc")
Omega$JJA <- ncvar_get(f, varid="omega")

################################-
# Long term Mean ----
################################-
  
Omega.longTerm <- Omega %>% 
  lapply(., function(x, lat){
    y <- apply(x, 1, mean) %>%  # average time
      cbind.data.frame(lat, .)
    return(y)
  },  lat)

########   plotting general circulation -

Omega.longTerm.g <- Omega.longTerm %>% reshape::melt(., id="lat") %>%
  within(., {
    value <-  value * 1e3
  })

range <- with(Omega.longTerm.g, max(abs(min(value)),max(value)) %>% round(., 2))

ggplot() + facet_wrap( .~ L1 , ncol=2 )+
  geom_line(data= Omega.longTerm.g, aes(x= lat, y= value))+ 
  scale_y_continuous(transform = "reverse" )+
  scale_x_continuous(breaks = seq(-90,90,by=15))+
  labs(x="Latitude [°]", y="Omega [mPa/s]", title = "Omega 500 hPa - ModE-RA ensamble mean")+
  theme_bw()

################################-
# Strong ascent ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

strAsc  <- list()
lat.min.model  <- list()
for (j in seasons){
  
  strAsc[[j]] <- smt.min.max(Omega[[j]] %>% t(),"min", lat) # identify the position of the strong Ascent with the splines
  lat.min.model[[j]] <- Omega[[j]] %>% apply(., 2, which.min) %>% lat[.] # identify it without the smoothing
  
}
strAsc  <- lapply(strAsc , function(x,dates) cbind.data.frame(x,dates), Dates )
lat.min.model  <- lapply(lat.min.model , function(x,dates) cbind.data.frame(x,dates), Dates )


################################-
## plotting features timeseries ----
################################-

lat.min.model.g <- lat.min.model %>% reshape::melt(., id=c("dates")) %>% within(., model <- "original")
Omega.prop.g <- strAsc %>% reshape::melt(., id=c("dates")) %>% within(., model <- "Smooth Spline") %>% 
  rbind(., lat.min.model.g)

Omega.prop.g %>%
  ggplot(., aes(x= dates, y= value, col=model)) +
  facet_wrap(. ~ L1,scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = c("red","black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y")+
  labs(title="ModE-RA Ens. Mean Strong Ascent [Omega 500 hPa]", y="Latitude [°]")+
  theme_bw()+theme(legend.position = "bottom")


