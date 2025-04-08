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
library(parallel)
library(doParallel)

source("00_settings.R")
setwd(dir.base)

################################-
## generate names in ModE-Sim data ----
################################-

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
                       Epoch= 2))

# names of ensembles in subdirectories
sets.ens <- data.frame( Set= paste0("set_",rep(1420,3),"-",seq(1,3)),
                        Memb= rep("ensmean", 3),
                        Epoch= 1, 
                        stringsAsFactors = F
) %>% 
  rbind(., data.frame( Set= paste0("set_",rep(1850,2),"-",seq(1,2)),
                       Memb= rep("ensmean", 2),
                       Epoch= 2)
  )

################################-
### Dates ----
################################-

Dates <- list()
Dates$all <- seq(as.Date("1420-01-01"), as.Date("2009-12-01"), by="year")
Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="month")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="month")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="month")
seasons <- c("DJF","JJA")

################################-
## load data ----
################################-
Volc <- read.csv("./01_Data/Volcanic_erup.csv")[,-1] %>% 
  melt() %>% 
  dplyr::mutate(Date= paste0(value,"-01-01") %>% as.Date(),
                variable = factor(variable, levels =c("Fischer","VEI5")))

Hum <<- list()
for(i in unique(sets$Set)) Hum[[i]] <- list()

Hum.ens <- list()

f <- nc_open("./01_Data/04_q/ModE-Sim_set_1420-1_to_3_q-ZonMean_1420-1849_sDJF.nc")
lat <<- f$dim$lat$vals; lon <<- f$dim$lon$vals; lev <<- f$dim$plev$vals

# function to read several subset of the ensamble from the seasonal data

read_Mod <- function (Set, Memb, Epoch, varn="q"){
  
  Dataset <- "ModE-Sim"
  
  if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
  
  f.mod <- paste0("./01_Data/04_q/",Set,"/",Dataset,"_",Set,"_",Memb,"_q-ZonMean_",t.span)
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  
  Hum[[Set]][[Memb]] <<- ncvar_get(f, varid=varn)
  
  # return(SLP)
}

print("reading members")

pb <- txtProgressBar(min=1,max=nrow(sets),style=3)
for (i in 1:nrow(sets)){
  setTxtProgressBar(pb,i)
  
  a <- with(sets, read_Mod(Set= Set[i], Memb= Memb[i], Epoch= Epoch[i]))
}
close(pb)

print("reading ensemble means")
varn <- "q"
for (i in 1:nrow(sets.ens)){
  
  Dataset <- "ModE-Sim"
  f <- with(sets.ens[i,],{
    
    if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
    f.mod <- paste0("./01_Data/04_q/",Set,"/",Dataset,"_",Set,"_",Memb,"_q-ZonMean_",t.span)

    paste0(f.mod,"_mon.nc") %>% nc_open(.)
  })
  
  Hum.ens[[sets.ens$Set[i]]] <- ncvar_get(f, varid=varn)
}

################################-
# Vertically Integrated Humidity ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

# create the delta pressures for the vertical integral
d.lev <- (lev - c(lev[ -1], 0))
d.lev[-1] <- d.lev[-length(d.lev)]/2 + d.lev[-1]/2
d.lev[1] <- d.lev[1]/2

VIM <- function(x, d.lev){
  
   VIM <- x %>% apply(., c(1,3), "*", d.lev) %>% # multiply each humidity by the corresponding delta pressure
    apply(., c(2,3), sum) %>% # summatory over the vertical pressure axis
    magrittr::divide_by(., 9.8066) # divide by gravity
   
   return(VIM)
}

VIHum <- list()
VIHum.trop <- list()

print("Calculating Vertical integrated humidity and filtering to tropics")
for( i in names(Hum)){
  print(i)
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  VIHum[[i]] <- mclapply(Hum[[i]], VIM, d.lev) # vertical integral of humidity for every subset
  
  VIHum.trop[[i]] <- mclapply(VIHum[[i]],
                            function(x) x[lat >= -30 & lat <= 30,] %>% # filter tropical humidity
                              apply(., 2, mean) 
                            )
  stopCluster(cl)
  
  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"}
  
  VIHum.trop[[i]] <- lapply(VIHum.trop[[i]], function(x,dates) cbind.data.frame(x,dates), Dates[[epoch]])
}

print("Same for the ensemble")
for( i in names(Hum)){
  print(i)
  
  VIHum[[i]][["ensmean"]] <- VIM(Hum.ens[[i]], d.lev)
  
  VIHum.trop[[i]][["ensmean"]] <- VIHum[[i]][["ensmean"]] %>% 
    .[lat >= -30 & lat <= 30,] %>% # filter tropical humidity
    apply(., 2, mean) 
  
  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"}
  VIHum.trop[[i]][["ensmean"]] <- cbind.data.frame(VIHum.trop[[i]][["ensmean"]],Dates[[epoch]])
  }

# save(VIHum,list=("VIHum"), file="./01_Data/04_q/ModE-Sim_allSets_VIHum.RData")
# save(VIHum.trop,list=c("VIHum.trop"), file="./01_Data/04_q/ModE-Sim_allSets_VIHum-trop.RData")

################################-
# DJF & JJA aggregation ----
################################-
VIHum.s <- list()

Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")

DJF_JJA <- function( VIHum, Dates.memb){
  
  Memb <- list()
  
  Memb$DJF <- VIHum %>% read.zoo(., index.column = 2) %>% dm2seasonal(., "DJF", FUN=mean, na.rm=T)
  Memb$JJA <- VIHum %>% read.zoo(., index.column = 2) %>% dm2seasonal(., "JJA", FUN=mean, na.rm=T)
  
  Memb <- lapply(Memb, `index<-`, Dates.memb) %>% lapply(., fortify.zoo)
  
  return(Memb)
}

# cl <- makeCluster(10)
# registerDoParallel(cl)

print("seasonal accumulation")
for ( i in names(VIHum.trop)){
  print(i)
  
  Epoch <- with (sets, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)
  
  VIHum.s[[i]] <- lapply(VIHum.trop[[i]], DJF_JJA, Dates[[Epoch]])
  
}
# stopCluster(cl)

save(VIHum.s,list="VIHum.s",file="./01_Data/04_q/_ModE-Sim_allMemb_VIHum-trop.RData")
################################-
## plotting features timeseries ----
################################-
load("./01_Data/04_q/_ModE-Sim_allMemb_VIHum-trop.RData")

palette <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")

# strength
Hum.trop <- VIHum.s %>% reshape::melt(., id=c("Index")) %>% magrittr::set_colnames(., c("dates","variable","value","Season","Memb","Dataset")) %>%
  within(.,{
    Dataset <- factor(Dataset, levels=names(VIHum.s))
  }) 

Hum.trop.g <- subset(Hum.trop, Memb !="ensmean") %>% 
  dplyr::group_by(., Dataset, dates, Season) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                   p95=quantile(value, probs= 0.95, na.rm = T),
                   min= min(value, na.rm = T),
                   max= max(value, na.rm = T)) %>% 
  as.data.frame()

Hum.trop.ens.g <- subset(Hum.trop, Memb =="ensmean")

# Hum.trop.30y <- Hum.trop.g %>% dplyr::group_by(dates) %>% 
#   dplyr::summarise(DJF = mean(value[Season =="DJF"]), JJA= mean(value[Season =="JJA"])) %>% read.zoo(., index.column = 1) %>% 
#   rollmean(., k=31) %>% 
#   fortify(melt=T) %>% magrittr::set_colnames(.,c("dates","Season","value"))

ggplot() +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
  geom_ribbon(data= Hum.trop.g, aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data= Hum.trop.ens.g, aes(x= dates, y= value, col=Dataset))+
  scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+

  labs(title="Variability of tropical humidity [Vert. Int. Hum. |lat|>=30°] - Subsets Ens. ", y="Hum [kg/m2]")+
  theme_bw()+theme(legend.position = c(0.2,0.5), legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="00"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))


ggplot() +
  geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
  geom_ribbon(data= subset(Hum.trop.g, Season=="DJF"), aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.2)+ scale_fill_manual(values = palette)+
  geom_line(data= subset(Hum.trop.ens.g, Season=="DJF"), aes(x= dates, y= value, col=Dataset))+
  scale_color_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  coord_cartesian(ylim=c(33,36.5))+
  
  labs(title="Variability of tropical humidity in DJF |lat|>=30° - Subsets Ens. ", y="Vert. Int. Hum. [kg/m2]")+
  theme_bw()+theme(legend.position = c(0.3,0.84), legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="00"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

ggsave(paste0("/scratch2/nduque/z_2025_PAGES/Trop_Hum_Zon-mean.png"),
       dpi=300,width = 1200*3/300,height = 350*3/300, units = "in")
