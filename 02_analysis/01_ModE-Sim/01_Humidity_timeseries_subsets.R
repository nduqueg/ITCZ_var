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

Hum <<- list()

f <- nc_open("./01_Data/04_q/ModE-Sim_set_1420-1_to_3_q-ZonMean_1420-1849_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals

# function to read several subset of the ensamble from the seasonal data

read_Mod <- function (f.mod, Set, varn="q"){
  
  Hum[[Set]] <<- list()
  if (substr(Set, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(Set, 5,8) =="1850"){ epoch <- "ep2"}
  
  f <- paste0(f.mod,"_sDJF.nc") %>% nc_open(.)
  lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals
  
  Hum[[Set]]$DJF <<- ncvar_get(f, varid=varn) %>% .[, ,-(length(Dates[[epoch]])+1)] # exclude the last Dec that belongs to the next year in season DJF
  
  f <- paste0(f.mod,"_sJJA.nc") %>% nc_open(.)
  Hum[[Set]]$JJA <<- ncvar_get(f, varid=varn)
  
  return(Hum)
}

print("reading subsets")
a <- read_Mod(f.mod= "./01_Data/04_q/set_1420-1/ModE-Sim_set_1420-1_ensmean_q-ZonMean_1420-1849",
              Set="set_1420-1")
a <- read_Mod(f.mod= "./01_Data/04_q/set_1420-2/ModE-Sim_set_1420-2_ensmean_q-ZonMean_1420-1849",
              Set="set_1420-2")
a <- read_Mod(f.mod= "./01_Data/04_q/set_1420-3/ModE-Sim_set_1420-3_ensmean_q-ZonMean_1420-1849",
              Set="set_1420-3")

a <- read_Mod(f.mod= "./01_Data/04_q/set_1850-1/ModE-Sim_set_1850-1_ensmean_q-ZonMean_1850-2009",
              Set="set_1850-1")
a <- read_Mod(f.mod= "./01_Data/04_q/set_1850-2/ModE-Sim_set_1850-2_ensmean_q-ZonMean_1850-2009",
              Set="set_1850-2")

################################-
# Vertically Integrated Humidity ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

VIM <- function(x, lev){
  
  # create the delta pressures for the vertical integral
  d.lev <- (lev - c(lev[ -1], 0))
  d.lev[-1] <- d.lev[-length(d.lev)]/2 + d.lev[-1]/2
  d.lev[1] <- d.lev[1]/2
  
   VIM <- x %>% apply(., c(1,3), "*", d.lev) %>% # multiply each humidity by the corresponding delta pressure
    apply(., c(2,3), sum) %>% # summatory over the vertical pressure axis
    magrittr::divide_by(., 9.8066) # divide by gravity
   
   return(VIM)
}

VIHum <- list()
VIHum.trop <- list()

print("Calculating position for subset:")
for( i in names(Hum)){
  print(i)
  
  VIHum[[i]] <- lapply(Hum[[i]], VIM, lev) # vertical integral of humidity for every subset
  
  VIHum.trop[[i]] <- lapply(VIHum[[i]],
                            function(x) x[lat >= -30 & lat <= 30,] %>% # filter tropical humidity
                              apply(., 2, mean) 
                            )

  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"}
  
  VIHum.trop[[i]] <- lapply(VIHum.trop[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
}

################################-
## plotting features timeseries ----
################################-

palette <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")

# strength
Hum.trop.g <- VIHum.trop %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>%
  within(.,{
    Dataset <- factor(Dataset, levels=names(Hum))
  })

Hum.trop.30y <- Hum.trop.g %>% dplyr::group_by(dates) %>% 
  dplyr::summarise(DJF = mean(value[Season =="DJF"]), JJA= mean(value[Season =="JJA"])) %>% read.zoo(., index.column = 1) %>% 
  rollmean(., k=31) %>% 
  fortify(melt=T) %>% magrittr::set_colnames(.,c("dates","Season","value"))

ggplot() +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_line(data= Hum.trop.g, aes(x= dates, y= value, col=Dataset))+ 
  geom_line(data= Hum.trop.30y, aes(x= dates, y= value), col="black")+
  scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+

  labs(title="Variability of tropical humidity [Vert. Int. Hum. |lat|>=30°] - Subsets Ens. ", y="Hum [kg/m2]")+
  theme_bw()+theme(legend.position = c(0.2,0.45), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

