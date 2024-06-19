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
Omega <- list(); Omega$ep1 <- list(); Omega$ep2 <- list()

f <- nc_open("./01_Data/02_Omega500/ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals

Omega$ep1$DJF <- ncvar_get(f, varid="omega") %>% .[, -(length(Dates$ep1)+1)] # exclude the last Dec that belongs to the next year in season DJF
  
f <- nc_open("./01_Data/02_Omega500/ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_sJJA.nc")
Omega$ep1$JJA <- ncvar_get(f, varid="omega")

f <- nc_open("./01_Data/02_Omega500/ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_sDJF.nc")
Omega$ep2$DJF <- ncvar_get(f, varid="omega") %>% .[, -(length(Dates$ep2)+1)] # exclude the last Dec that belongs to the next year in season DJF

f <- nc_open("./01_Data/02_Omega500/ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_sJJA.nc")
Omega$ep2$JJA <- ncvar_get(f, varid="omega")

################################-
# Long term Mean ----
################################-

Omega.longTerm <-list()
for( i in c("ep1","ep2")){ # two time periods
  
  Omega.longTerm[[i]] <- Omega[[i]] %>% 
    lapply(., function(x, lat){
      y <- apply(x, 1, mean) %>%  # average time
        cbind.data.frame(lat, .)
      return(y)
    },  lat)
}

########   plotting general circulation -

Omega.longTerm.g <- Omega.longTerm %>% reshape::melt(., id="lat") %>%
  within(., {
    value <-  value * 1e3
  })

range <- with(Omega.longTerm.g, max(abs(min(value)),max(value)) %>% round(., 2))

ggplot() + facet_grid( L1 ~ L2, switch="y" )+
  geom_line(data= Omega.longTerm.g, aes(x= lat, y= value))+ 
  scale_y_continuous(transform = "reverse" )+
  scale_x_continuous(breaks = seq(-90,90,by=15))+
  labs(x="Latitude [°]", y="Omega [mPa/s]", title = "Omega 500 hPa - ModE-Sim ensamble mean")+
  theme_bw()

################################-
# Strong ascent ----
################################-
smt.min.max <- function(x,type){
  # 'x' is the Omega at 500 hPa, which will be smoothed around the minimum with an spline curve and 
  # the position of its minimum will be also identified
  
  aprox <- apply(x,1, paste0("which.",type) %>% get())
  ind.aprox <-cbind(aprox-3,aprox-2,aprox-1,aprox,aprox+1,aprox+2,aprox+3) %>% t() %>% as.data.frame() # indices of values close to the Max StrFunct for the Spline model
  lat.models <- apply(ind.aprox, 2, function(ind,lat) lat[ind], lat) %>%  as.data.frame() # latitudes for the Spline model
  
  data.spModel <- mapply(FUN= function(x,ind) x[ind], # extracts the MAx StrFunct values near the Max for the Spline model
                         x %>% as.matrix() %>% split(.,row(.)), 
                         ind.aprox,
                         SIMPLIFY = F) %>% 
    
    mapply(function(y,lat) cbind.data.frame(lat,y), # paste for each time step the latitude and the values near the Max StrFunc
           .,
           lat.models, 
           SIMPLIFY = F) %>% 
    lapply(., `colnames<-`,c("x.mod","y.mod"))
  
  sp.models <- lapply(data.spModel, function(x) lm(y.mod ~ splines::ns(x.mod,3), data= x))
  pos.mastrf <- numeric()
  for(k in 1:length(sp.models)){ # for each year identify the position of the Strong Ascent
    pos.mastrf[k] <- predict(sp.models[[k]], newdata = data.frame(x.mod= seq(lat.models[7,k],lat.models[1,k],by=0.1))) %>% as.numeric(.) %>% matrix(.,nrow=1) %>% 
      apply(.,1,get(paste0("which.",type))) %>% 
      seq(lat.models[7,k],lat.models[1,k],by=0.1)[.]
  }
  return(pos.mastrf)
}

strAsc <- list()
lat.min.model <- list()

for( i in c("ep1","ep2")){ # two time periods
  strAsc[[i]] <- list()
  lat.min.model[[i]] <- list()
  for (j in seasons){
    
    strAsc[[i]][[j]] <- smt.min.max(Omega[[i]][[j]] %>% t(),"min")
    lat.min.model[[i]][[j]] <- Omega[[i]][[j]] %>% apply(., 2, which.min) %>% lat[.]
    
  }
  strAsc[[i]] <- lapply(strAsc[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[i]])
  lat.min.model[[i]] <- lapply(lat.min.model[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[i]])
}

################################-
## plotting features timeseries ----
################################-

lat.min.model.g <- lat.min.model %>% reshape::melt(., id=c("dates")) %>% within(., model <- "original")
Omega.prop.g <- strAsc %>% reshape::melt(., id=c("dates")) %>% within(., model <- "Smooth Spline") %>% 
  rbind(., lat.min.model.g)

Omega.prop.g %>%
  ggplot(., aes(x= dates, y= value, col=model)) +
  facet_wrap(. ~ L2,scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = c("red","black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y")+
  labs(title="ModE-Sim Ens. Mean Strong Ascent [Omega 500 hPa]", y="Latitude [°]")+
  theme_bw()+theme(legend.position = "bottom")


