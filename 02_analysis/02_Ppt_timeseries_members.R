rm(list=ls())
cat("\014")

"%>%"=magrittr::`%>%`
library(ncdf4)
library(abind)
library(hydroTSM)
library(reshape)
library(ggplot2)
library(metR)
library(RColorBrewer)

source("00_settings.R")
setwd(dir.base)

################################-
### Dates ----
################################-

Dates <- list()
Dates$all <- seq(as.Date("1420-01-01"), as.Date("2009-12-01"), by="month")
Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="month")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="month")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="month")
seasons <- c("DJF","JJA")

################################-
## generate names in ModE-Sim data ----
################################-

Memb.aux <- character()
for ( i in 1:100){ if (i>=10 & i <100) Memb.aux[i] <- paste0("m0",i) else if (i < 10) Memb.aux[i] <- paste0("m00",i) else Memb.aux[i] <- paste0("m",i)}

# names of directories
sets <- data.frame( Set= paste0("set_",rep(1420,3),"-",seq(1,3)) %>% rep(.,each=20),
                    Memb= paste0("m00",seq(1,9)) %>% c(.,paste0("m0",seq(10,60))),
                    Epoch= 1,
                    stringsAsFactors = F
) %>% 
  rbind(., data.frame( Set= "set_1850-1" %>% rep(., each=20),
                       Memb= paste0("m00",seq(1,9)) %>% c(.,paste0("m0",seq(10,20))),
                       Epoch= 2)
  ) %>% 
  rbind(., data.frame( Set= "set_1850-2" %>% rep(., each=16),
                       Memb= paste0("m0",seq(21,36)),
                       Epoch= 2)
  ) %>% rbind(., data.frame( Set= "ModE-RA",
                             Memb= paste0("m0",seq(41,60)),
                             Epoch= 2)
  ) %>% rbind(., data.frame( Set= "ModE-RAclim",
                             Memb= Memb.aux,
                             Epoch= 2)
  )

# names of ensembles in subdirectories
sets.ens <- data.frame( Set= paste0("set_",rep(1420,3),"-",seq(1,3)),
                        Memb= rep("ensmean", 3),
                        Epoch= 1, 
                        stringsAsFactors = F
) %>% 
  rbind(., data.frame( Set= paste0("set_",rep(1850,2),"-",seq(1,2)),
                       Memb= rep("ensmean", 2),
                       Epoch= 2)
  ) %>% 
  rbind(., data.frame( Set= c("ModE-RA", "ModE-RAclim"),
                       Memb= rep("ensmean",2),
                       Epoch= rep(2,2))
  )

################################-
## load data ----
################################-
Volc <- read.csv("./01_Data/Volcanic_erup.csv")[,-1] %>% 
  melt() %>% 
  dplyr::mutate(Date= paste0(value,"-01-01") %>% as.Date(),
                variable = factor(variable, levels =c("Fischer","VEI5")))

Ppt <<- list()
for(i in unique(sets$Set)) Ppt[[i]] <- list()

f <- nc_open("./01_Data/05_Ppt/ModE-Sim_set_1420-1_to_3_totprec-ZonMean_1420-1849_mon.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals; lev <- f$dim$plev$vals

# function to read several members of the ensamble

read_Mod <- function (Set, Memb, Epoch, varn="totprec"){
  Dataset <- "ModE-Sim"
  
  Ppt[[Set]][[Memb]] <<- list()
  
  if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
  
  f.mod <- paste0("./01_Data/05_Ppt/",Set,"/",Dataset,"_",Set,"_",Memb,"_totprec-ZonMean_",t.span)
  if (Set == "ModE-RA"){
    Dataset <- "ModE-RA"
    t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/05_Ppt/lowres_20mem_Set_1420-3_1850-1/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_",Memb,"_totprec-ZonMean_",t.span)
  } else if (Set == "ModE-RAclim"){
    Dataset <- "ModE-RAclim"
    t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/05_Ppt/lowres_100mem_Set_1/",Dataset,"_",Memb,"_totprec-ZonMean_",t.span)
  }
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  Ppt[[Set]][[Memb]] <<- ncvar_get(f, varid=varn)
  
  return(Ppt)
}

# ----------------------------- reading members .............
print("reading members")
pb <- txtProgressBar(min=1,max=nrow(sets),style=3)
for (i in 1:nrow(sets)){
  setTxtProgressBar(pb,i)
  
  a <- with(sets, read_Mod(Set= Set[i], Memb= Memb[i], Epoch= Epoch[i]))
}
close(pb)

print("reading ensemble means")

for (i in 1:nrow(sets.ens)){
  
  P.memb <- with(sets.ens[i,],{
    
    Dataset <- "ModE-Sim"
    if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
    f.mod <- paste0("./01_Data/05_Ppt/",Set,"/",Dataset,"_",Set,"_",Memb,"_totprec-ZonMean_",t.span)
    
    if (Set == "ModE-RA"){
      Dataset <- "ModE-RA"
      t.span <- "1421-2008"
      f.mod <- paste0("./01_Data/05_Ppt/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_totprec-ZonMean_",t.span)
    } else if(Set == "ModE-RAclim"){
      Dataset <- "ModE-RAclim"
      t.span <- "1421-2008"
      f.mod <- paste0("./01_Data/05_Ppt/",Dataset,"_lowres_100mem_Set_1_totprec-ZonMean_",t.span)
    }
    
    data <- list(f.mod=f.mod, Set=Set)
    return(data)
  })
  
  Memb <- "ensmean"
  Set <- P.memb$Set;  f.mod <- P.memb$f.mod
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  Ppt[[Set]][[Memb]] <- ncvar_get(f, varid="totprec")
}

################################-
# Maximum and mean zonal precipitation in the tropics ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

# loc.max.Ppt <- max.Ppt <- mean.Ppt <- list()
# 
# Lat.trop.fil <- lat[lat >= -25 & lat <= 25]
# 
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
# 
# print("Calculating position for subset:")
# for( i in names(Ppt)){
#   print(i)
#   
#   loc.max.Ppt[[i]] <- list()
#   max.Ppt[[i]] <- list()
#   mean.Ppt[[i]] <- list()
#   
#   for (j in seasons){
#     loc.max.Ppt[[i]][[j]] <- smt.min.max(Ppt[[i]][[j]][ lat >= -25 & lat <= 25, ] %>% t(),"max", Lat.trop.fil)
#     max.Ppt[[i]][[j]] <- apply(Ppt[[i]][[j]][lat >= -25 & lat <= 25,], 2, max)
#     mean.Ppt[[i]][[j]] <- mean.lat(Ppt[[i]][[j]][lat >= -25 & lat <= 25,] %>% t(), Lat.trop.fil)
#   } 
#   
#   if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
#   loc.max.Ppt[[i]] <- lapply(loc.max.Ppt[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
#   max.Ppt[[i]] <- lapply(max.Ppt[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
#   mean.Ppt[[i]] <- lapply(mean.Ppt[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
# }

################################-
## plotting features timeseries ----
################################-

# location
# Ppt.loc.g <- loc.max.Ppt %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>% 
#   within(., Dataset <- factor(Dataset, levels=names(Ppt))  )
# 
# palette <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","black")
# Ppt.loc.g %>%
#   ggplot(., aes(x= dates, y= value, col=Dataset)) +
#   facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
#   geom_line()+ scale_color_manual(values = palette)+
#   scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
#                             date_labels = "%Y", expand=c(0.01,0.01))+
#   labs(title="Position Max. Ppt - Subsets Ens. ", y="Latitude [°]")+
#   theme_bw()+theme(legend.position = c(0.2,0.1), legend.direction = "horizontal",
#                    panel.grid = element_line(linetype="dashed",color="lightgrey"),
#                    axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
# 
# # strength
# Ppt.str.g <- max.Ppt %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>%
#   within(.,{
#     Dataset <- factor(Dataset, levels=names(Ppt))
#   })
# 
# Ppt.str.g %>%
#   ggplot(., aes(x= dates, y= value, col=Dataset)) +
#   facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
#   geom_line()+ scale_color_manual(values = palette)+
#   scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
#                date_labels = "%Y", expand=c(0.01,0.01))+
#   labs(title="Maximum zonal Precipitation - Subsets Ens. ", y="Ppt  [mm/month]")+
#   theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
#                    panel.grid = element_line(linetype="dashed",color="lightgrey"),
#                    strip.text = element_text(size=12),
#                    axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
# 
# # mean Ppt
# Ppt.mean.g <- mean.Ppt %>% reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset")) %>%
#   within(.,{
#     Dataset <- factor(Dataset, levels=names(Ppt))
#   })
# 
# Ppt.mean.g %>%
#   ggplot(., aes(x= dates, y= value, col=Dataset)) +
#   facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
#   geom_line()+ scale_color_manual(values = palette)+
#   scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
#                date_labels = "%Y", expand=c(0.01,0.01))+
#   labs(title="Mean zonal tropical Precipitation [25°S , 25°N] - Subsets Ens. ", y="Ppt  [mm/month]")+
#   theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
#                    panel.grid = element_line(linetype="dashed",color="lightgrey"),
#                    strip.text = element_text(size=12),
#                    axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))


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
    
    trop.Ppt[[i]][[j]] <- lapply( Ppt[[j]], # apply over each member
                                  function(x, Lat.trop.fil){
                                    z <- mean.lat(x[ trop.fil[,i], ] %>% t(), Lat.trop.fil)
                                    return(z)
                                    },
                                  lat[trop.fil[,i]])
    
    if (substr(j, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(j, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
    
    trop.Ppt[[i]][[j]] <- lapply(trop.Ppt[[i]][[j]], # apply over each member
                                 function(x, dates){
                                   z <- cbind.data.frame(x, dates); return(z)
                                 } , Dates[[epoch]])
  }
}

################################-
# DJF & JJA aggregation and anomalies ----
################################-

P.s <- list()

DJF_JJA <- function( P){
  
  Memb <- list()
  
  Memb$DJF <- P %>% read.zoo(., index.column=2) %>% dm2seasonal(., "DJF", FUN=sum, na.rm=T) %>% fortify.zoo() %>% within(., Index <- as.numeric(Index))
  Memb$JJA <- P %>% read.zoo(., index.column=2) %>% dm2seasonal(., "JJA", FUN=sum, na.rm=T) %>% fortify.zoo() %>% within(., Index <- as.numeric(Index))
  
  
  return(Memb)
}

library(parallel)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

# Seasonal accumulation
print("seasonal accumulation")
for ( i in names(trop.Ppt)){
  print(i)
  
  P.s[[i]]<- list()
  for(j in names(Ppt))  P.s[[i]][[j]] <- mclapply(trop.Ppt[[i]][[j]], DJF_JJA)
  
}
stopCluster(cl)

# 1850-DJF and 1420-DJF in ModE-Sim have a lack of Ppt due to the not inclusion of December previous year
# the we multiply by 3/2
P.s <- lapply(P.s, function(band){ 
  y <- lapply(band, function(set){
    z <- lapply(set, function(memb){
      memb$DJF[1,2] <- memb$DJF[1,2] * 3/2
      return(memb)
    })
    return(z)
  })
  return(y)
})

# calculate anomalies
P.anom <- list()
calc.anom <- function(x){
  if (x$DJF[1,1] < 1850){
    
    x.anom <- lapply(x, function(y) y[,-1] %>% scale(., scale= FALSE) %>% as.data.frame())
  } else {
    
    x.anom <- lapply(x, function(y){( y[,-1] - mean( y[1:51, -1]))  %>% as.data.frame()}) # from 1850 to 1900
  }
  dates <- x$DJF[,1]
  x.anom <- lapply(x.anom, cbind, dates)
  x.anom <- lapply(x.anom, `colnames<-`, c("value","dates"))
  return(x.anom)
}

print("seasonal anomalies")
for(j in names(trop.Ppt)){
  
  print(j)
  P.anom[[j]] <- list()
  for ( i in names(Ppt)[-7] ){
    # anomalies in 1420-1850 and 1850-1900
    # periods without high antropogenic alteration
    
    P.anom[[j]][[i]] <- lapply(P.s[[j]][[i]], calc.anom)
  }
  
  P.anom[[j]][["ModE-RAclim"]] <- P.s[[j]][["ModE-RAclim"]] %>% lapply(., function(x){ 
    x <- lapply(x, function(y){ 
      colnames(y) <- c("dates","value")
      y$dates <- as.numeric(as.character(y$dates))
      return(y)})
    return(x)})
}


################################-
# plot mean precipitation in tropical bands ----
################################-
# trop.Ppt.g <- trop.Ppt %>% 
#   reshape::melt(., id="dates") %>% magrittr::set_colnames(., c("dates","variable","value","Season","Dataset", "Band")) %>% 
#   subset(., select= - variable) %>% 
#   within(., {
#     Dataset <- factor(Dataset, levels=names(Ppt))
#     Band <- factor(Band, levels= c("north","deep","south" ))
#     }) 
# 
# Band.labs <- c("North [13°N,25°N]", "Deep trop. [13°S,13°N]", "South [25°S,13°S]"); names(Band.labs) <- c("north","deep","south")
# 
# trop.Ppt.g %>% 
#   ggplot(., aes(x= dates, y= value, col=Dataset)) +
#   facet_grid(Band ~ Season, scales = "free_y", switch = "y", labeller= labeller(Band = Band.labs))+
#   geom_line()+ scale_color_manual(values = palette)+
#   scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
#                date_labels = "%Y", expand=c(0.01,0.01))+
#   labs(title="Mean zonal tropical Precipitation in bands - Subsets Ens. ", y="Ppt  [mm/month]")+
#   theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
#                    panel.grid = element_line(linetype="dashed",color="lightgrey"),
#                    strip.text = element_text(size=12),
#                    axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
# 
# trop.Ppt.g %>% subset(., Season == "DJF") %>% 
#   
#   ggplot(., aes(x= dates, y= value, col=Dataset)) +
#   facet_grid(Band ~ Season, scales = "free_y", switch = "y", labeller= labeller(Band = Band.labs))+
#   geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
#   geom_line()+ scale_color_manual(values = palette)+
#   scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
#                date_labels = "%Y", expand=c(0.01,0.01))+
#   labs(title="Mean zonal tropical Precipitation in bands - Subsets Ens. ", y="Ppt  [mm/month]")+
#   theme_bw()+theme(legend.position = c(0.2,0.05), legend.direction = "horizontal",
#                    panel.grid = element_line(linetype="dashed",color="00"),
#                    strip.text = element_text(size=12),
#                    axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

################################-
# anomalies precipitation in tropical bands ----
################################-
P.anom.sets <- P.anom
P.anom <- lapply(P.anom, function(x){ names(x)[1:5] <- "ModE-Sim"; return(x)})

P.anom <- P.anom %>% 
  reshape::melt(., id=c("dates","value")) %>% magrittr::set_colnames(., c("dates","value","Season","Memb","Dataset", "Band")) %>% 
  within(., {
    Dataset <- factor(Dataset, levels=c("ModE-Sim","ModE-RAclim","ModE-RA")) 
    Band <- factor(Band, levels= c("north","deep","south" ))
    })

Memb.anom.g <- subset(P.anom, Memb !="ensmean") %>% 
  dplyr::group_by(., Dataset, dates, Season, Band) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                   p95=quantile(value, probs= 0.95, na.rm = T),
                   min= min(value, na.rm = T),
                   max= max(value, na.rm = T)) %>% 
  as.data.frame()

Ens.anom.g <- subset(P.anom, Memb =="ensmean") %>% 
  dplyr::group_by(., Dataset, dates, Season, Band) %>%
  dplyr::summarise(value.m = mean(value, na.rm=T)) %>% 
  as.data.frame()

Volc <- within(Volc, Date <- format(Date,format="%Y") %>% as.numeric())

palette <- c("#ff7f00","#377eb8","#999999")
Band.labs <- c("North [13°N,25°N]", "Deep trop. [13°S,13°N]", "South [25°S,13°S]"); names(Band.labs) <- c("north","deep","south")

ggplot( ) +
  facet_grid(Band ~ Season, scales = "free_y", switch = "y", labeller= labeller(Band = Band.labs))+
  geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
  geom_ribbon(data=Memb.anom.g, aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=Ens.anom.g,aes(x= dates, y= value.m, col=Dataset))+ scale_color_manual(values = c(palette[-3],"black"))+
  scale_x_continuous(breaks = seq(1450,2000,by=50 ),  expand=c(0.01,0.01))+
  labs(title="Anomalies zonal tropical Precipitation in bands ", y="Ppt anomalies [mm/season]")+
  theme_bw()+theme(legend.position = c(0.7,0.05), legend.direction = "horizontal",
                   panel.grid = element_blank(),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

ggplot( ) +
  facet_wrap(.~Band, scales = "free_y", ncol=1, labeller= labeller(Band = Band.labs))+
  # facet_grid(Band ~ Season, scales = "free_y", switch = "y", labeller= labeller(Band = Band.labs))+
  geom_vline(data= Volc,aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
  geom_ribbon(data=Memb.anom.g %>% subset(., Season=="DJF"),
              aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=Ens.anom.g %>% subset(., Season=="DJF"),
            aes(x= dates, y= value.m, col=Dataset))+ scale_color_manual(values = c(palette[-3],"black"))+
  scale_x_continuous(breaks = seq(1450,2000,by=50 ),  expand=c(0.01,0.01))+
  labs(title="Anomalies zonal tropical Precipitation in DJF ", y="Ppt anomalies [mm/season]")+
  theme_bw()+theme(legend.position = c(0.2,0.4), legend.direction = "horizontal",
                   panel.grid = element_blank(),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
ggsave(paste0("/scratch2/nduque/z_2025_PAGES/Ppt_trop_Bands.png"),
       dpi=300,width = 1400*3/300,height = 700*3/300, units = "in")
