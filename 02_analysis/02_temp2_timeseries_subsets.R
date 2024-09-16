rm(list=ls())
cat("\014")

"%>%"=magrittr::`%>%`
library(ncdf4)
library(abind)
library(hydroTSM)
library(ggplot2)
library(metR)
library(RColorBrewer)
library(reshape)

source("00_settings.R")
setwd(dir.base)

################################-
### Dates ----
################################-

Dates <- list()
Dates$all <- seq(as.Date("1420-01-01"), as.Date("2009-12-01"), by="year")
Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")
seasons <- c("DJF","JJA")

################################-
## load data ----
################################-

Temp <<- list()

f <- nc_open("./01_Data/03_temp2/ModE-Sim_set_1420-1_to_3_temp2-ZonMean_1420-1849_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals

# function to read several subset of the ensamble from the seasonal data

read_Mod <- function (f.mod, Set, varn="temp2"){
  
  Temp[[Set]] <<- list()
  if (substr(Set, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(Set, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  
  f <- paste0(f.mod,"_sDJF.nc") %>% nc_open(.)
  lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals
  
  Temp[[Set]]$DJF <<- ncvar_get(f, varid=varn) %>%
    .[, -(length(Dates[[epoch]])+1)] %>% # exclude the last Dec that belongs to the next year in season DJF
    t()
  
  f <- paste0(f.mod,"_sJJA.nc") %>% nc_open(.)
  Temp[[Set]]$JJA <<- ncvar_get(f, varid=varn) %>% t()
  
  # define the vector of Dates for the monthly zoo and then creating the annual time series
  if (epoch =="ep1"){ 
    Dates.annual <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"),by="month")
  }else if (epoch =="ep2"){
    Dates.annual <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"),by="month")
  }else if(epoch=="ModERA"){
    Dates.annual <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"),by="month")
    }
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  Temp[[Set]]$annual <<- ncvar_get(f, varid=varn) %>%
    t(.) %>% zoo(., order.by = Dates.annual) %>% 
    monthly2annual(., FUN=mean, na.rm=T) %>% as.matrix(.)
  
  return(Temp)
}

print("reading subsets")
a <- read_Mod(f.mod= "./01_Data/03_temp2/set_1420-1/ModE-Sim_set_1420-1_ensmean_temp2-ZonMean_1420-1849",
              Set="set_1420-1")
a <- read_Mod(f.mod= "./01_Data/03_temp2/set_1420-2/ModE-Sim_set_1420-2_ensmean_temp2-ZonMean_1420-1849",
              Set="set_1420-2")
a <- read_Mod(f.mod= "./01_Data/03_temp2/set_1420-3/ModE-Sim_set_1420-3_ensmean_temp2-ZonMean_1420-1849",
              Set="set_1420-3")

a <- read_Mod(f.mod= "./01_Data/03_temp2/set_1850-1/ModE-Sim_set_1850-1_ensmean_temp2-ZonMean_1850-2009",
              Set="set_1850-1")
a <- read_Mod(f.mod= "./01_Data/03_temp2/set_1850-2/ModE-Sim_set_1850-2_ensmean_temp2-ZonMean_1850-2009",
              Set="set_1850-2")

a <- read_Mod(f.mod= "./01_Data/03_temp2/ModE-RA_lowres_20mem_Set_1420-3_1850-1_temp2-ZonMean_1421-2008",
              Set="ModE-RA")

################################-
# temperature at 2 m ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

North <- list()
South <- list()
Int.Hem <- list()
Tropical <- list()

sqrt(cos(3.14159*lat/180))

mean.lat <- function(x, lat.f){
  
  total.weight <- sqrt(cos(3.14159*lat.f/180)) %>% sum()
  
  if( length(lat.f) == dim(x)[2]){
    y <- sweep(x, 2, FUN="*", sqrt(cos(3.14159*lat.f/180))) %>% # applying the weights
      apply(.,1, sum, na.rm=T) %>% 
      magrittr::divide_by(., total.weight ) # normalizing the weights
    
  } else{
    stop("not equal latitude intervals in both dataset and latitude vector")
  } 
    
  return(y)
}

print("Calculating inter-hemispheric difference by subset:")
for( i in names(Temp)){
  print(i)
  
  North[[i]] <- list()
  South[[i]] <- list()
  Int.Hem[[i]] <- list()
  Tropical[[i]] <- list()
  
  # doing the calculation by the seasons or the annual mean
  for (j in c(seasons, "annual")){
    North[[i]][[j]] <- Temp[[i]][[j]][, lat >= 30] %>% mean.lat(., lat.f = lat[lat >= 30])
    
    South[[i]][[j]] <- Temp[[i]][[j]][, lat <= -30] %>% mean.lat(., lat.f = lat[lat <= -30])
    
    Int.Hem[[i]][[j]] <- North[[i]][[j]] - South[[i]][[j]]
    
    Tropical[[i]][[j]] <- Temp[[i]][[j]][, lat >= -30 & lat <= 30] %>% mean.lat(.,  lat.f = lat[lat >= -30 & lat <= 30])
  } 
  
  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  
  North[[i]] <- lapply(North[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
  South[[i]] <- lapply(South[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
  Int.Hem[[i]] <- lapply(Int.Hem[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
  Tropical[[i]] <- lapply(Tropical[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
  
}

################################-
## Re arrange timeseries fro plotting ----
################################-

Data.g <- list()

for( Series in c("North", "South", "Int.Hem", "Tropical"))  
  Data.g[[Series]] <- get(Series) %>% 
  melt(., id="dates") %>% 
  cbind(.,Series)

Data.g <- do.call(rbind,Data.g)[,-2] %>% magrittr::set_colnames(., c("dates","value","Season","Dataset","Series")) %>% 
  within(., Dataset <- factor(Dataset, levels=names(Temp))  )

################################-
## plotting features timeseries ----
################################-

palette <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","black")

# plot of the Inter-hemispheric Temperature gradient
subset(Data.g, Series == "Int.Hem") %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Inter-hemispheric Temperature gradient [extra tropical >= 30°] - Subsets Ens. ", y="Temperature [°K]")+
  theme_bw()+theme(legend.position = c(0.2,0.95), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# plot of individual regions
subset(Data.g, Series != "Int.Hem") %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_grid(Series ~ Season, switch ="y", scales="free_y")+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Temperature - Subsets Ens. ", y="Temperature [°K]")+
  theme_bw()+theme(legend.position = c(0.2,0.95), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# plot of individual regions
subset(Data.g, Series != "Int.Hem" & Series != "Tropical" & Season !="annual") %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_grid(Series ~ Season, switch ="y", scales="free_y")+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Temperature - Subsets Ens. ", y="Temperature [°K]")+
  theme_bw()+theme(legend.position = c(0.2,0.1), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
