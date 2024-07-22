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
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")
seasons <- c("DJF","JJA")


## load data ----
Omega <<- list()

f <- nc_open("./01_Data/02_Omega500/ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals

read_Mod <- function (f.mod, Set, varn="omega"){
  
  Omega[[Set]] <<- list()
  if (substr(Set, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(Set, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  
  f <- paste0(f.mod,"_sDJF.nc") %>% nc_open(.)
  lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals
  
  Omega[[Set]]$DJF <<- ncvar_get(f, varid=varn) %>% .[, -(length(Dates[[epoch]])+1)] # exclude the last Dec that belongs to the next year in season DJF
  
  f <- paste0(f.mod,"_sJJA.nc") %>% nc_open(.)
  Omega[[Set]]$JJA <<- ncvar_get(f, varid=varn)
  
  return(Omega)
}

print("reading subsets")
a <- read_Mod(f.mod= "./01_Data/02_Omega500/set_1420-1/ModE-Sim_set_1420-1_ensmean_Omega500hPa-ZonMean_1420-1849",
              Set="set_1420-1")
a <- read_Mod(f.mod= "./01_Data/02_Omega500/set_1420-2/ModE-Sim_set_1420-2_ensmean_Omega500hPa-ZonMean_1420-1849",
              Set="set_1420-2")
a <- read_Mod(f.mod= "./01_Data/02_Omega500/set_1420-3/ModE-Sim_set_1420-3_ensmean_Omega500hPa-ZonMean_1420-1849",
              Set="set_1420-3")

a <- read_Mod(f.mod= "./01_Data/02_Omega500/set_1850-1/ModE-Sim_set_1850-1_ensmean_Omega500hPa-ZonMean_1850-2009",
              Set="set_1850-1")
a <- read_Mod(f.mod= "./01_Data/02_Omega500/set_1850-2/ModE-Sim_set_1850-2_ensmean_Omega500hPa-ZonMean_1850-2009",
              Set="set_1850-2")

a <- read_Mod(f.mod= "./01_Data/02_Omega500/ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008",
              Set="ModE-RA")

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

print("Calculating position for subset:")
for( i in names(Omega)){
  print(i)
  
  strAsc[[i]] <- list()
  lat.min.model[[i]] <- list()
  
  for (j in seasons) strAsc[[i]][[j]] <- smt.min.max(Omega[[i]][[j]] %>% t(),"min")
  
  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  strAsc[[i]] <- lapply(strAsc[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
}

################################-
## plotting features timeseries ----
################################-


Omega.prop.g <- strAsc %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset"))

palette <- brewer.pal(9, "Set1")[-c(6:8)]
Omega.prop.g %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Position Mean Strong Ascent [Omega 500 hPa] - Subsets Ens. ", y="Latitude [°]")+
  theme_bw()+theme(legend.position = c(0.2,0.1), legend.direction = "horizontal",
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))


