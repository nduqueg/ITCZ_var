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

Omega <<- list()

f <- nc_open("./01_Data/02_Omega500/ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals

# function to read several subset of the ensamble from the seasonal data

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
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

loc.strAsc <- list()
strAsc <- list()
lat.min.model <- list()

print("Calculating position for subset:")
for( i in names(Omega)){
  print(i)
  
  loc.strAsc[[i]] <- list()
  strAsc[[i]] <- list()
  
  for (j in seasons){
    loc.strAsc[[i]][[j]] <- smt.min.max(Omega[[i]][[j]] %>% t(),"min", lat)
    strAsc[[i]][[j]] <- apply(Omega[[i]][[j]][lat >= -45 & lat <= 45,], 2, min)
  } 
  
  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  loc.strAsc[[i]] <- lapply(loc.strAsc[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
  strAsc[[i]] <- lapply(strAsc[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
}

################################-
## plotting features timeseries ----
################################-

# location
Omega.loc.g <- loc.strAsc %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>% 
  within(., Dataset <- factor(Dataset, levels=names(Omega))  )

palette <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","black")
Omega.loc.g %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Position Max. Strong Ascent [Omega 500 hPa] - Subsets Ens. ", y="Latitude [°]")+
  theme_bw()+theme(legend.position = c(0.2,0.1), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# strength
Omega.str.g <- strAsc %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>%
  within(.,{
    Dataset <- factor(Dataset, levels=names(Omega))
    value <- value * 100
  })

Omega.str.g %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  scale_y_continuous(transform = "reverse")+
  labs(title="Maximum Strong Ascent [Min. Omega 500 hPa] - Subsets Ens. ", y="Omega 500 hPa [hPa/s]")+
  theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

