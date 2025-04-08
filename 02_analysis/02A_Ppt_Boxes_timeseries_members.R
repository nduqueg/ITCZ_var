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
Dates$all <- seq(as.Date("1420-01-01"), as.Date("2009-12-01"), by="year")
Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")
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
  within(., variable <- factor(variable, levels =c("Fischer","VEI5"))) %>% 
  magrittr::set_colnames(.,c("variable","Date"))

P <<- list()

for(i in unique(sets$Set)) P[[i]] <- list()

read_Mod <- function (Set, Memb, Epoch){
  
  Dataset <- "ModE-Sim"
  
  if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
  
  f.mod <- paste0("./01_Data/05A_Ppt_boxes/",Set,"/",Dataset,"_",Set,"_",Memb,"_totprec-Boxes_",t.span)
  
  if (Set == "ModE-RA"){
    Dataset <- "ModE-RA"
    t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/05A_Ppt_boxes/lowres_20mem_Set_1420-3_1850-1/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_",Memb,"_totprec-Boxes_",t.span)
  } else if (Set == "ModE-RAclim"){
    Dataset <- "ModE-RAclim"
    t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/05A_Ppt_boxes/lowres_100mem_Set_1/",Dataset,"_",Memb,"_totprec-Boxes_anom_",t.span)
  }
  
  load(paste0(f.mod,"_mon.RData"))
  P[[Set]][[Memb]] <<- P.memb
  
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
P.ens <<- list()
for (i in 1:nrow(sets.ens)){
  
  Dataset <- "ModE-Sim"
  P.memb <- with(sets.ens[i,],{
    
    if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
    f.mod <- paste0("./01_Data/05A_Ppt_boxes/",Set,"/",Dataset,"_",Set,"_",Memb,"_totprec-Boxes_",t.span)
    
    if (Set == "ModE-RA"){
      Dataset <- "ModE-RA"
      t.span <- "1421-2008"
      f.mod <- paste0("./01_Data/05A_Ppt_boxes/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_totprec-Boxes_",t.span)
    } else if(Set == "ModE-RAclim"){
      Dataset <- "ModE-RAclim"
      t.span <- "1421-2008"
      f.mod <- paste0("./01_Data/05A_Ppt_boxes/",Dataset,"_lowres_100mem_Set_1_totprec-Boxes_anom_",t.span)
    }
    
    load(paste0(f.mod,"_mon.RData"))
    return(P.memb)
  })
  P.ens[[sets.ens$Set[i]]] <- P.memb
}
rm(P.memb)

################################-
# DJF & JJA aggregation and anomalies ----
################################-

P.s <- list()
P.ens.s <- list()

DJF_JJA <- function( P){
  
  Memb <- list()
  
  Memb$DJF <- P %>% read.zoo(., index.column=1) %>% dm2seasonal(., "DJF", FUN=sum, na.rm=T) %>% fortify.zoo()
  Memb$JJA <- P %>% read.zoo(., index.column=1) %>% dm2seasonal(., "JJA", FUN=sum, na.rm=T) %>% fortify.zoo()
  
  return(Memb)
}

library(parallel)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)

# Seasonal accumulation
print("seasonal accumulation")
for ( i in names(P)){
  print(i)
  
  P.s[[i]] <- mclapply(P[[i]], DJF_JJA)
  
}
stopCluster(cl)

print("seasonal accumulation Ensemble Means")
for ( i in names(P.ens)){
  print(i)
  
  # Epoch <- with (sets.ens, Epoch[i==Set][1])
  # if (i=="ModE-RA") Epoch <- "ModERA" else if (i=="ModE-RAclim") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)
  
  P.ens.s[[i]] <- DJF_JJA(P.ens[[i]])
}

# calculate anomalies
P.anom <- list()
P.ens.anom <- list()
calc.anom <- function(x){
  if (x$DJF[1,1] < 1850){
    
    x.anom <- lapply(x, function(y) y[,-1] %>% scale(., scale= FALSE) %>% as.data.frame())
  } else {
    
    x.anom <- lapply(x, function(y){( t(y[,-1]) - colMeans( y[1:51, -1])) %>% t() %>% as.data.frame()}) # from 1850 to 1900
  }
  dates <- x$DJF[,1]
  x.anom <- lapply(x.anom, cbind, dates)
  return(x.anom)
}

print("seasonal anomalies")
for ( i in names(P.s)[-7] ){
  print(i)
  
  # anomalies in 1420-1850 and 1850-1900
  # periods without high antropogenic alteration
    
  P.anom[[i]] <- lapply(P.s[[i]], calc.anom)
}
P.anom[["ModE-RAclim"]] <- P.s[["ModE-RAclim"]] %>% lapply(., function(x){ 
  x <- lapply(x, function(y){ 
    colnames(y)[1] <- "dates"
    y$dates <- as.numeric(as.character(y$dates))
    return(y)})
  return(x)})

print("seasonal anomlaies Ensemble Means")
for ( i in names(P.ens)[-7]){
  print(i)
  
  P.ens.anom[[i]] <- calc.anom(P.ens.s[[i]])
}
P.ens.anom[["ModE-RAclim"]] <- P.ens.s[["ModE-RAclim"]] %>% lapply(., function(y){ 
  colnames(y)[1] <- "dates"
  y$dates <- as.numeric(as.character(y$dates))
  return(y)})


# Plotting ModE-RAclim time series ----
################################-

# small modification in names of subdatasets

names(P.anom)[1:5] <- names(P.ens.anom)[1:5] <- "ModE-Sim"

# reorganize data and calculate bands
P.anom.g <- P.anom %>% melt(., id="dates") %>% magrittr::set_colnames(.,c("dates","Region","value","Season","memb","Set")) %>% 
  dplyr::group_by(., Set, dates, Season, Region) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                   p95=quantile(value, probs= 0.95, na.rm = T),
                   min= min(value, na.rm = T),
                   max= max(value, na.rm = T)) %>% 
  as.data.frame() %>% within(., Set <- factor(Set, levels=c("ModE-Sim","ModE-RAclim","ModE-RA")))

P.ens.anom.g <- P.ens.anom %>% melt(., id="dates") %>% magrittr::set_colnames(.,c("dates","Region","value.ens","Season","Set")) %>%
  dplyr::group_by(., Set, dates, Season, Region) %>%
  dplyr::summarise(value=mean(value.ens, na.rm = T)) %>% 
  as.data.frame() %>% within(., Set <- factor(Set, levels=c("ModE-Sim","ModE-RAclim","ModE-RA")))

palette <- c("#ff7f00","#377eb8","#999999")

ggplot( ) +
  facet_grid(Region ~ Season, scales = "free_y",switch="y")+
  geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
  geom_ribbon(data=P.anom.g, aes(x= dates, fill=Set,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=P.ens.anom.g, aes(x= dates, y=value, color=Set))+ scale_color_manual(values = c(palette[-3],"black"))+
  scale_x_continuous(breaks = seq(1450,2000,by=50), expand=c(0.01,0.01))+
  labs(title="Precipitation anomalies in proxies regions Ens. Memb.", y="Ppt [mm/season]")+
  theme_bw()+theme(legend.position = "bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="00"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
# 1400 x 700

for(i in levels(P.ens.anom.g$Region)){
  p <- ggplot( ) +
    facet_grid(Season ~ Region, scales = "free_y",switch="y")+
    geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
    geom_ribbon(data= subset(P.anom.g, Region==i), aes(x= dates, fill=Set,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
    geom_line(data= subset(P.ens.anom.g, Region==i), aes(x= dates, y=value, color=Set))+ scale_color_manual(values = c( palette[-3],"black"))+
    scale_x_continuous(breaks = seq(1450,2000,by=50), expand=c(0.01,0.01))+
    labs(title=paste0("Precipitation anomalies in proxies - ",i," - Ens. Memb."), y="Ppt [mm/season]")+
    theme_bw()+theme(legend.position = "bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                     panel.grid = element_line(linetype="dashed",color="00"),
                     axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
  ggsave(paste0("./pre-Figures/05_Ppt_prox_loc/05_ensMemb_Ppt_prox_Loc_",i,".png"),
         plot=p,
         dpi=100,width = 1400/100,height = 700/100, units = "in")
}
# 1400 x 700