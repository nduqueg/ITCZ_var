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

Ppt <<- list()

f <- nc_open("./01_Data/05_Ppt/ModE-Sim_set_1420-1_to_3_totprec-ZonMean_1420-1849_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals

# function to read several subset of the ensamble from the seasonal data

read_Mod <- function (f.mod, Set, varn="totprec"){
  
  Ppt[[Set]] <<- list()
  if (substr(Set, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(Set, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  
  f <- paste0(f.mod,"_sDJF.nc") %>% nc_open(.)
  lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals
  
  Ppt[[Set]]$DJF <<- ncvar_get(f, varid=varn) %>% .[, -(length(Dates[[epoch]])+1)] # exclude the last Dec that belongs to the next year in season DJF
  
  f <- paste0(f.mod,"_sJJA.nc") %>% nc_open(.)
  Ppt[[Set]]$JJA <<- ncvar_get(f, varid=varn)
  
  return(Ppt)
}

print("reading subsets")
a <- read_Mod(f.mod= "./01_Data/05_Ppt/set_1420-1/ModE-Sim_set_1420-1_ensmean_totprec-ZonMean_1420-1849",
              Set="set_1420-1")
a <- read_Mod(f.mod= "./01_Data/05_Ppt/set_1420-2/ModE-Sim_set_1420-2_ensmean_totprec-ZonMean_1420-1849",
              Set="set_1420-2")
a <- read_Mod(f.mod= "./01_Data/05_Ppt/set_1420-3/ModE-Sim_set_1420-3_ensmean_totprec-ZonMean_1420-1849",
              Set="set_1420-3")

a <- read_Mod(f.mod= "./01_Data/05_Ppt/set_1850-1/ModE-Sim_set_1850-1_ensmean_totprec-ZonMean_1850-2009",
              Set="set_1850-1")
a <- read_Mod(f.mod= "./01_Data/05_Ppt/set_1850-2/ModE-Sim_set_1850-2_ensmean_totprec-ZonMean_1850-2009",
              Set="set_1850-2")

a <- read_Mod(f.mod= "./01_Data/05_Ppt/ModE-RA_lowres_20mem_Set_1420-3_1850-1_totprec-ZonMean_1421-2008",
              Set="ModE-RA")

Ppt <- lapply(Ppt, function(set.ppt) lapply(set.ppt, function(x){
  y <- x * 86400 * 30 # transform units to mm/month
  return(y)
}))

################################-
# Strong ascent ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

loc.Ppt <- list()
max.Ppt <- list()

Lat.trop.fil <- lat[lat >= -25 & lat <= 25]

print("Calculating position for subset:")
for( i in names(Ppt)){
  print(i)
  
  loc.Ppt[[i]] <- list()
  max.Ppt[[i]] <- list()
  
  for (j in seasons){
    loc.Ppt[[i]][[j]] <- smt.min.max(Ppt[[i]][[j]][ lat >= -25 & lat <= 25, ] %>% t(),"max", Lat.trop.fil)
    max.Ppt[[i]][[j]] <- apply(Ppt[[i]][[j]][lat >= -25 & lat <= 25,], 2, max)
  } 
  
  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  loc.Ppt[[i]] <- lapply(loc.Ppt[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
  max.Ppt[[i]] <- lapply(max.Ppt[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
}

################################-
## plotting features timeseries ----
################################-

# location
Ppt.loc.g <- loc.Ppt %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>% 
  within(., Dataset <- factor(Dataset, levels=names(Ppt))  )

palette <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","black")
Ppt.loc.g %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Position Max. Ppt - Subsets Ens. ", y="Latitude [°]")+
  theme_bw()+theme(legend.position = c(0.2,0.1), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# strength
Ppt.str.g <- max.Ppt %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>%
  within(.,{
    Dataset <- factor(Dataset, levels=names(Ppt))
  })

Ppt.str.g %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Maximum Precipitation - Subsets Ens. ", y="Ppt  [mm/month]")+
  theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

