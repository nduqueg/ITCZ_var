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
  rbind(., data.frame( Set= c("ModE-RA"),
                       Memb= c("ensmean"),
                       Epoch= c(2))
  )

################################-
## load data ----
################################-
Volc <- read.csv("./01_Data/Volcanic_erup.csv")[,-1] %>% 
  melt() %>% 
  within(., variable <- factor(variable, levels =c("Fischer","VEI5"))) %>% 
  magrittr::set_colnames(.,c("variable","Date"))

Temp <<- list()

f <- nc_open("./01_Data/03_temp2/ModE-Sim_set_1420-1_to_3_temp2-ZonMean_1420-1849_sDJF.nc")
lat <- f$dim$lat$vals; lon <- f$dim$lon$vals

# function to read several subset of the ensamble from the seasonal data

read_Mod <- function (Set, Memb, varn="temp2"){
  Dataset <- "ModE-Sim"
  
  Temp[[Set]][[Memb]] <<- list()
  if (substr(Set, 5,8) =="1420"){
    epoch <- "ep1"; t.span <- "1420-1849"
    f.mod <- paste0("./01_Data/03_temp2/",Set,"/",Dataset,"_",Set,"_",Memb,"_temp2-ZonMean_",t.span)
    
  } else if(substr(Set, 5,8) =="1850"){
    epoch <- "ep2"; t.span <- "1850-2009"
    f.mod <- paste0("./01_Data/03_temp2/",Set,"/",Dataset,"_",Set,"_",Memb,"_temp2-ZonMean_",t.span)
    
  } else if(Set == "ModE-RA"){ 
    Dataset <- "ModE-RA"; epoch <- "ModERA"; t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/03_temp2/lowres_20mem_Set_1420-3_1850-1/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_",Memb,"_temp2-ZonMean_",t.span)
    
  }else if(Set == "ModE-RAclim"){
    Dataset <- "ModE-RAclim"; epoch <- "ModERA"; t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/03_temp2/lowres_20mem_Set_1420-3_1850-1/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_",Memb,"_temp2-ZonMean_",t.span)
    
  }
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  Temp[[Set]][[Memb]] <<- ncvar_get(f, varid=varn)
  
  return(f)
}

# ----------------------------- reading members .............
print("reading members")
pb <- txtProgressBar(min=1,max=nrow(sets),style=3)
for (i in 1:nrow(sets)){
  setTxtProgressBar(pb,i)
  
  a <- with(sets, read_Mod(Set= Set[i], Memb= Memb[i]))
}
close(pb)


print("reading ensemble means")

for (i in 1:nrow(sets.ens)){
  
  T.memb <- with(sets.ens[i,],{
    
    Dataset <- "ModE-Sim"
    if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
    f.mod <- paste0("./01_Data/03_temp2/",Set,"/",Dataset,"_",Set,"_",Memb,"_temp2-ZonMean_",t.span)
    
    if (Set == "ModE-RA"){
      Dataset <- "ModE-RA"
      t.span <- "1421-2008"
      f.mod <- paste0("./01_Data/03_temp2/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_temp2-ZonMean_",t.span)
    } else if(Set == "ModE-RAclim"){
      Dataset <- "ModE-RAclim"
      t.span <- "1421-2008"
      f.mod <- paste0("./01_Data/03_temp2/",Dataset,"_lowres_100mem_Set_1_temp2-ZonMean_",t.span)
    }
    
    data <- list(f.mod=f.mod, Set=Set)
    return(data)
  })
  
  Memb <- "ensmean"
  Set <- T.memb$Set;  f.mod <- T.memb$f.mod
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  Temp[[Set]][[Memb]] <- ncvar_get(f, varid="temp2")
}


################################-
# temperature at 2 m, interhemispheric ----
################################-
# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

North <- list()
South <- list()
Int.Hem <- list()
Tropical <- list()

# function for average considering the latitude weighting
mean.lat <- function(x, lat.f){
  
  x <- x[match(lat.f,lat),]
  total.weight <- cos(3.14159*lat.f/180) %>% sum()
  
  if( length(lat.f) == dim(x)[1]){
    y <- sweep(x, 1, FUN="*", cos(3.14159*lat.f/180)) %>% # applying the weights
      apply(.,2, sum, na.rm=T) %>% 
      magrittr::divide_by(., total.weight ) # normalizing the weights
    
  } else{
    stop("not equal latitude intervals in both dataset and latitude vector")
  } 
    
  return(y)
}

print("Calculating inter-hemispheric difference by subset:")
for( i in names(Temp)){
  print(i)
  
  North[[i]] <- list(); South[[i]] <- list(); Int.Hem[[i]] <- list(); Tropical[[i]] <- list()
  
  # doing the calculation by the seasons or the annual mean
  North[[i]] <- lapply(Temp[[i]], mean.lat, lat.f = lat[lat >= 30]) %>% as.data.frame()
  
  South[[i]] <- lapply(Temp[[i]], mean.lat, lat.f = lat[lat <= -30]) %>% as.data.frame()
  
  Int.Hem[[i]] <- mapply(North[[i]], South[[i]], FUN=function(x,y){ z <- x - y; return(z)}) %>% as.data.frame()
  
  Tropical[[i]] <- lapply(Temp[[i]], mean.lat, lat.f = lat[lat <= 30 & lat>= -30])  %>% as.data.frame()
  
  if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
  North[[i]] <- cbind.data.frame(North[[i]], Dates[[epoch]])
  South[[i]] <- cbind.data.frame(South[[i]], Dates[[epoch]])
  Int.Hem[[i]] <- cbind.data.frame(Int.Hem[[i]], Dates[[epoch]])
  Tropical[[i]] <- cbind.data.frame(Tropical[[i]], Dates[[epoch]])
}

################################-
# Seasonal aggregation ----
################################-

DJF_JJA <- function( Tmp.m){
  
  Memb <- list()
  
  Memb$DJF <- Tmp.m %>% read.zoo(., index.column="Dates[[epoch]]") %>% dm2seasonal(., "DJF", FUN=mean, na.rm=T) %>% fortify.zoo() %>% within(., Index <- as.numeric(Index))
  Memb$JJA <- Tmp.m %>% read.zoo(., index.column="Dates[[epoch]]") %>% dm2seasonal(., "JJA", FUN=mean, na.rm=T) %>% fortify.zoo() %>% within(., Index <- as.numeric(Index))
  
  return(Memb)
}

# library(parallel)
# library(doParallel)
# cl <- makeCluster(20)
# registerDoParallel(cl)

North.s <- list(); South.s <- list(); Int.Hem.s <- list(); Tropical.s <- list()
# Seasonal accumulation
print("seasonal accumulation")
for ( i in names(Temp)){
  print(i)
  
  North.s[[i]] <- DJF_JJA(North[[i]])
  South.s[[i]] <- DJF_JJA(South[[i]])
  Int.Hem.s[[i]] <- DJF_JJA(Int.Hem[[i]])
  Tropical.s[[i]] <- DJF_JJA(Tropical[[i]])
  
}
# stopCluster(cl)

# Dates <- list()
# Dates$all <- seq(as.Date("1420-01-01"), as.Date("2009-12-01"), by="year")
# Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
# Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
# Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")
# 

################################-
## Re arrange timeseries fro plotting ----
################################-

Data.g <- list()

for( Series in c("North.s", "South.s", "Int.Hem.s", "Tropical.s")){
  Data.g[[Series]] <- get(Series) %>% 
    melt(., id="Index") %>% 
    cbind(.,Series)
}

Data.memb <- do.call(rbind,Data.g) %>% magrittr::set_colnames(., c("dates","Memb","value","Season","Dataset","Series")) %>% 
  subset(., Memb!="ensmean") %>% 
  within(., Dataset <- factor(Dataset, levels=names(Temp))  ) %>% 
  dplyr::group_by(., Dataset, dates, Season, Series) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                   p95=quantile(value, probs= 0.95, na.rm = T),
                   min= min(value, na.rm = T),
                   max= max(value, na.rm = T)) %>% 
  as.data.frame()

Data.ens <- do.call(rbind,Data.g) %>% magrittr::set_colnames(., c("dates","Memb","value","Season","Dataset","Series")) %>% 
  subset(., Memb=="ensmean") %>% 
  within(., Dataset <- factor(Dataset, levels=names(Temp))  )

################################-
## plotting features timeseries ----
################################-

palette <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","black")

# plot of the Inter-hemispheric Temperature gradient
ggplot( ) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_ribbon(data=Data.memb %>% subset(., Series=="Int.Hem.s"), aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=Data.ens %>% subset(., Series=="Int.Hem.s"), aes(x= dates, y= value, col=Dataset))+ scale_color_manual(values = palette)+
  scale_x_continuous(breaks = seq(1450,2000,by=50), expand=c(0.01,0.01))+
  labs(title="Inter-hemispheric Temperature gradient [extra tropical >= 30°] - Subsets Ens. ", y="Temperature [°K]")+
  theme_bw()+theme(legend.position = c(0.2,0.95), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
ggsave("pre-Figures/Interhemis_gradient.png",
       dpi=300,width = 1400*3/300,height = 850*3/300, units = "in")


# plot JJA
Data.memb1 <- Data.memb %>% subset(., Series=="Int.Hem.s" & Season=="JJA" & dates>=1780 & dates <=1850)
Data.ens1 <- Data.ens %>% subset(., Series=="Int.Hem.s" & Season=="JJA" & dates>=1780 & dates <=1850)
Volc <- Volc %>% subset(., Date>1780 & Date <1850)
ggplot( ) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_ribbon(data=Data.memb1, aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=Data.ens1, aes(x= dates, y= value, col=Dataset))+ scale_color_manual(values = palette)+
  geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
  scale_x_continuous(breaks = seq(1780,1850,by=10), expand=c(0.01,0.01))+
  labs(title="Inter-hemispheric Temperature gradient [extra tropical >= 30°]", y="Temperature [°K]")+
  theme_bw()+theme(legend.position = c(0.8,0.2), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
ggsave("/scratch2/nduque/z_2025_PAGES/T_Interhemis_gradient_1780-1850.png",
       dpi=300,width = 1400*3/300,height = 450*3/300, units = "in")

# # plot of individual regions
# subset(Data.g, Series != "Int.Hem" & Series != "Tropical" & Season !="annual") %>%
#   ggplot(., aes(x= dates, y= value, col=Dataset)) +
#   facet_grid(Series ~ Season, switch ="y", scales="free_y")+
#   geom_line()+ scale_color_manual(values = palette)+
#   scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
#                date_labels = "%Y", expand=c(0.01,0.01))+
#   labs(title="Temperature - Subsets Ens. ", y="Temperature [°K]")+
#   theme_bw()+theme(legend.position = c(0.2,0.1), legend.direction = "horizontal",
#                    panel.grid = element_line(linetype="dashed",color="lightgrey"),
#                    strip.text = element_text(size=12),
#                    axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
# 
# 
# ################################-
# # temperature at 2 m,  whole world ----
# ################################-
# World <- list()
# 
# for( i in names(Temp)){
#   print(i)
#   
#   World[[i]] <- Temp[[i]][["annual"]] %>% mean.lat(., lat.f = lat)
#   
#   if (substr(i, 5,8) =="1420"){ epoch <- "ep1" } else if(substr(i, 5,8) =="1850"){ epoch <- "ep2"} else { epoch <- "ModERA" }
#   
#   dates <- Dates[[epoch]]
#   World[[i]] <- cbind.data.frame(World[[i]], dates)
# }
# 
# data.g <- World %>% 
#   melt(., id="dates") %>% 
#   magrittr::set_colnames(., c("dates","variable","value","Dataset")) %>% 
#   within(., Dataset <- factor(Dataset, levels=names(Temp))  )
# 
# data.g %>% 
#   ggplot(., aes(x= dates, y= value, col=Dataset)) +
#   geom_line()+ scale_color_manual(values = palette)+
#   scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
#                date_labels = "%Y", expand=c(0.01,0.01))+
#   labs(title="Temperature - Subsets Ens. ", y="Temperature [°K]")+
#   theme_bw()+theme(legend.position = c(0.2,0.95), legend.direction = "horizontal",
#                    panel.grid = element_line(linetype="dashed",color="lightgrey"),
#                    strip.text = element_text(size=12),
#                    axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
# 
# ################################-
# # anomalies hemispheres ----
# ################################-
# 
# Data.anom <- list()
# anom.1500.1850 <- function(x){
#   y <- read.zoo(x, index.column = "dates")
#   
#   if (x$dates[1] < "1850-01-01") p.anom <- seq(as.Date("1500-01-01"), as.Date("1850-12-01"), by="year") else p.anom <- seq(as.Date("1850-01-01"), as.Date("1900-12-01"), by="year")
#   
#     anom <- (y - mean( y[p.anom] )) %>% 
#     fortify.zoo()
#   return(anom)
# }
# 
# for( Series in c("North", "South", "Tropical")){
#   Data.anom[[Series]] <- get(Series) %>% 
#     lapply(., function(x) lapply(x,anom.1500.1850) ) %>% 
#     lapply(., function(x) lapply(x,`colnames<-`,c("dates","value")))
# }
# 
# Data.anom %>% 
#   melt(., id=c("dates","value")) %>% magrittr::set_colnames(., c("dates","value","Season","Dataset", "Series")) %>% 
#   within(., Dataset <- factor(Dataset, levels=names(Temp)) ) %>% 
#   ggplot(., aes(x= dates, y= value, col=Dataset))+
#   facet_grid(Series ~ Season, switch = "y")+
#   geom_line()+ scale_color_manual(values = palette)+
#   scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
#                date_labels = "%Y", expand=c(0.01,0.01))+
#   labs(title="Temperature anomaly- Subsets Ens. ", y="Temperature [°K]")+
#   theme_bw()+theme(legend.position = c(0.2,0.3), legend.direction = "horizontal",
#                    panel.grid = element_line(linetype="dashed",color="lightgrey"),
#                    strip.text = element_text(size=12),
#                    axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
# 
