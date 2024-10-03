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
# Maximum and mean zonal precipitation in the tropics ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

loc.max.Ppt <- max.Ppt <- mean.Ppt <- list()

Lat.trop.fil <- lat[lat >= -25 & lat <= 25]

mean.lat <- function(x, lat.f){
  
  total.weight <- cos(3.14159*lat.f/180) %>% sum()
  
  if( length(lat.f) == dim(x)[2]){
    y <- sweep(x, 2, FUN="*", cos(3.14159*lat.f/180)) %>% # applying the weights
      apply(.,1, sum, na.rm=T) %>% 
      magrittr::divide_by(., total.weight ) # normalizing the weights
    
  } else{
    stop("not equal latitude intervals in both dataset and latitude vector")
  } 
  
  return(y)
}

print("Calculating position for subset:")
for( i in names(Ppt)){
  print(i)
  
  loc.max.Ppt[[i]] <- list()
  max.Ppt[[i]] <- list()
  mean.Ppt[[i]] <- list()
  
  for (j in seasons){
    loc.max.Ppt[[i]][[j]] <- smt.min.max(Ppt[[i]][[j]][ lat >= -25 & lat <= 25, ] %>% t(),"max", Lat.trop.fil)
    max.Ppt[[i]][[j]] <- apply(Ppt[[i]][[j]][lat >= -25 & lat <= 25,], 2, max)
    mean.Ppt[[i]][[j]] <- mean.lat(Ppt[[i]][[j]][lat >= -25 & lat <= 25,] %>% t(), Lat.trop.fil)
  } 
  
  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  loc.max.Ppt[[i]] <- lapply(loc.max.Ppt[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
  max.Ppt[[i]] <- lapply(max.Ppt[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
  mean.Ppt[[i]] <- lapply(mean.Ppt[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
}

################################-
## plotting features timeseries ----
################################-

# location
Ppt.loc.g <- loc.max.Ppt %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>% 
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
  labs(title="Maximum zonal Precipitation - Subsets Ens. ", y="Ppt  [mm/month]")+
  theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# mean Ppt
Ppt.mean.g <- mean.Ppt %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>%
  within(.,{
    Dataset <- factor(Dataset, levels=names(Ppt))
  })

Ppt.mean.g %>%
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Mean zonal tropical Precipitation [25°S , 25°N] - Subsets Ens. ", y="Ppt  [mm/month]")+
  theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))


################################-
# mean precipitation in tropical bands ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

trop.Ppt <- list()

trop.fil <- data.frame(north= lat >= 13 & lat <= 25,
                       deep= lat >= -13 & lat <= 13,
                       south= lat >= -25 & lat <= -13)

print("Calculating mean Ppt in tropical bands:")
for( i in c("north","deep","south" ) ){
  print(i)
  
  trop.Ppt[[i]] <- list()
  
  for (j in names(Ppt)){
    
    trop.Ppt[[i]][[j]] <- lapply( Ppt[[j]], # apply over each season
                                  function(x, Lat.trop.fil){
                                    z <- mean.lat(x[ trop.fil[,i], ] %>% t(), Lat.trop.fil)
                                    return(z)
                                    },
                                  lat[trop.fil[,i]])
    
    if (substr(j, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(j, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
    
    trop.Ppt[[i]][[j]] <- lapply(trop.Ppt[[i]][[j]], # apply over each season
                                 function(x, dates){
                                   z <- cbind.data.frame(x, dates); return(z)
                                 } , Dates[[epoch]])
  }
}

################################-
# plot mean precipitation in tropical bands ----
################################-
trop.Ppt.g <- trop.Ppt %>% 
  reshape::melt(., id="dates") %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset", "Band")) %>% 
  subset(., select= - variable) %>% 
  within(., {
    Dataset <- factor(Dataset, levels=names(Ppt))
    Band <- factor(Band, levels= c("north","deep","south" ))
    }) 

Band.labs <- c("North [13°N,25°N]", "Deep trop. [13°S,13°N]", "South [25°S,13°S]"); names(Band.labs) <- c("north","deep","south")

trop.Ppt.g %>% 
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_grid(Band ~ Season, scales = "free_y", switch = "y", labeller= labeller(Band = Band.labs))+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Mean zonal tropical Precipitation in bands - Subsets Ens. ", y="Ppt  [mm/month]")+
  theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

################################-
# anomalies precipitation in tropical bands ----
################################-

Data.anom <- list()
anom.1500.1850 <- function(x){
  y <- read.zoo(x, index.column = "dates")
  
  if (x$dates[1] < "1850-01-01") p.anom <- seq(as.Date("1500-01-01"), as.Date("1850-12-01"), by="year") else p.anom <- seq(as.Date("1850-01-01"), as.Date("1900-12-01"), by="year")
  
  anom <- (y - mean( y[p.anom] )) %>% 
    fortify.zoo()
  return(anom)
}

for( Series in c("north", "deep", "south")){
  Data.anom[[Series]] <- trop.Ppt[[Series]] %>% 
    lapply(., function(x) lapply(x,anom.1500.1850) ) %>% 
    lapply(., function(x) lapply(x,`colnames<-`,c("dates","value")))
}

Data.anom %>% 
  reshape::melt(., id=c("dates","value")) %>% magrittr::set_colnames(., c("dates","value","Season","Dataset", "Band")) %>% 
  within(., {
    Dataset <- factor(Dataset, levels=names(Ppt)) 
    Band <- factor(Band, levels= c("north","deep","south" ))
    }) %>% 
  ggplot(., aes(x= dates, y= value, col=Dataset)) +
  facet_grid(Band ~ Season, scales = "free_y", switch = "y", labeller= labeller(Band = Band.labs))+
  geom_line()+ scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Anomalies zonal tropical Precipitation in bands - Subsets Ens. ", y="Ppt  [mm/month]")+
  theme_bw()+theme(legend.position = c(0.7,0.05), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
