rm(list=ls())
cat("\014")

"%>%"=magrittr::`%>%`
library(ncdf4)
library(abind)
library(hydroTSM)
library(ggplot2)
library(reshape)
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
  rbind(., data.frame( Set= "ModE-RA",
                       Memb= "ensmean",
                       Epoch= 2)
        )

### Dates ----
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

SLP <<- list() 
for(i in unique(sets$Set)) SLP[[i]] <- list()

SLP.ens <- list()

f <- nc_open("./01_Data/06A_SLP_ocean/ModE-Sim_set_1420-1_to_3_slp-ZonMean_1420-1849_sDJF.nc")
lat <<- f$dim$lat$vals; lon <<- f$dim$lon$vals; lev <<- f$dim$plev$vals


read_Mod <- function (Set, Memb, Epoch, varn="slp"){
  
  Dataset <- "ModE-Sim"
  
  if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
  
  f.mod <- paste0("./01_Data/06A_SLP_ocean/",Set,"/",Dataset,"_",Set,"_",Memb,"_slp-ZonMean_",t.span)
  
  if (Set == "ModE-RA"){
    Dataset <- "ModE-RA"
    t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/06A_SLP_ocean/lowres_20mem_Set_1420-3_1850-1/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_",Memb,"_slp-ZonMean_",t.span)
  }
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  
  SLP[[Set]][[Memb]] <<- ncvar_get(f, varid=varn)

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
varn <- "slp"
for (i in 1:nrow(sets.ens)){
  
  Dataset <- "ModE-Sim"
  f <- with(sets.ens[i,],{
    
    if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
    f.mod <- paste0("./01_Data/06A_SLP_ocean/",Set,"/",Dataset,"_",Set,"_",Memb,"_slp-ZonMean_",t.span)
    if (Set == "ModE-RA"){
      Dataset <- "ModE-RA"
      t.span <- "1421-2008"
      f.mod <- paste0("./01_Data/06A_SLP_ocean/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_slp-ZonMean_",t.span)
    }
    
    paste0(f.mod,"_mon.nc") %>% nc_open(.)
  })
  SLP.ens[[sets.ens$Set[i]]] <- ncvar_get(f, varid=varn)
}

# save(SLP, file="./01_Data/06A_SLP_ocean/_allMemb_slp.Rdata")
  
################################-
# DJF & JJA aggregation ----
################################-
SLP.s <- list()
SLP.ens.s <- list()

DJF_JJA <- function( SLP, Dates.memb){
  
  Memb <- list()

  Memb$DJF <- t(SLP) %>% zoo(., order.by= Dates.memb) %>% dm2seasonal(., "DJF", FUN=mean, na.rm=T) %>% as.matrix()
  Memb$JJA <- t(SLP) %>% zoo(., order.by= Dates.memb) %>% dm2seasonal(., "JJA", FUN=mean, na.rm=T) %>% as.matrix()
  
  return(Memb)
}

cl <- makeCluster(10)
registerDoParallel(cl)

print("seasonal accumulation")
for ( i in names(SLP)){
  print(i)
  
  Epoch <- with (sets, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)
  
  SLP.s[[i]] <- mclapply(SLP[[i]], DJF_JJA, Dates[[Epoch]])

}
stopCluster(cl)

print("seasonal accumulation Ensemble Means")
for ( i in names(SLP.ens)){
  print(i)
  
  Epoch <- with (sets.ens, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)
  
  SLP.ens.s[[i]] <- DJF_JJA(SLP.ens[[i]], Dates[[Epoch]])
}
# save(SLP.s, file="./01_Data/06A_SLP_ocean/_allMemb_slp_DJFJJA.Rdata")

################################-
# Strong ascent ----
################################-
# load("./01_Data/06A_SLP_ocean/_allMemb_slp_DJFJJA.Rdata")

# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

# function for calculating the feature for both seasons with the smoothing spline
min.max.memb <- function(y, Lat, feature ="loc", Trop.fil){
  
  # filter to tropics
  for (j in seasons) y[[j]] <- y[[j]][, Trop.fil]
  Lat.trop.fil <- Lat[Trop.fil]
  
  Memb <- list()
  if (feature =="loc"){
    Memb[["DJF"]] <- smt.min.max(y[["DJF"]], "max", Lat.trop.fil)
    Memb[["JJA"]] <- smt.min.max(y[["JJA"]], "max", Lat.trop.fil)  
  } else if(feature == "Str"){
    Memb[["DJF"]] <- apply(y[["DJF"]], 1, max)
    Memb[["JJA"]] <- apply(y[["JJA"]], 1, max)
  }
  
  
  return(Memb)
}

loc.subThigh <- list()
subThigh <- list()

trop.fil <- lat >= 2 & lat <= 60

Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")

print("identifying position")
for( i in names(SLP.s)){ # two time periods
  print(i)
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  loc.subThigh[[i]] <- mclapply(SLP.s[[i]], min.max.memb, lat, "loc", trop.fil)
  subThigh[[i]] <- mclapply(SLP.s[[i]], min.max.memb, lat, "Str", trop.fil)
  stopCluster(cl)
  
  Epoch <- with (sets, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)  # defines the Epoch for attaching the dates to the time series of the members
  
  loc.subThigh[[i]] <- lapply(loc.subThigh[[i]], 
                        function(x,dates){ y <- as.data.frame(x); z <- cbind.data.frame(dates,y); return(z)},
                        Dates[[Epoch]])
  subThigh[[i]] <- lapply(subThigh[[i]], 
                            function(x,dates){ y <- as.data.frame(x); z <- cbind.data.frame(dates,y); return(z)},
                            Dates[[Epoch]])
}

print("identifying position Ensemble means")
for( i in names(SLP.ens.s)){ # two time periods
  print(i)

  Epoch <- with (sets.ens, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)  # defines the Epoch for attaching the dates to the time series of the members
  dates <- Dates[[Epoch]]
  
  loc.subThigh[[i]][["ensmean"]] <- min.max.memb(SLP.ens.s[[i]], lat, "loc", trop.fil) %>% cbind.data.frame(dates, .)
  subThigh[[i]][["ensmean"]] <- min.max.memb(SLP.ens.s[[i]], lat, "Str", trop.fil) %>% cbind.data.frame(dates, .)
  
}

save(loc.subThigh,file="./01_Data/06A_SLP_ocean/_allMemb_LocsubThigh.Rdata")
save(subThigh,file="./01_Data/06A_SLP_ocean/_allMemb_subThigh.Rdata")

################################-
## plotting features timeseries ----
################################-
load("./01_Data/06A_SLP_ocean/_allMemb_LocsubThigh.Rdata")
load("./01_Data/06A_SLP_ocean/_allMemb_subThigh.Rdata")

loc.subThigh <- loc.subThigh %>% 
  reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","Season","value","Memb","Dataset"))

SLP.loc.g <- subset(loc.subThigh, Memb !="ensmean") %>% 
  dplyr::group_by(., Dataset, dates, Season) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                  p95=quantile(value, probs= 0.95, na.rm = T),
                  min= min(value, na.rm = T),
                  max= max(value, na.rm = T)) %>% 
  as.data.frame() %>% within(., Dataset <- factor(Dataset, levels=unique(sets$Set)))

SLP.loc.ens.g <- subset(loc.subThigh, Memb =="ensmean") %>% within(., Dataset <- factor(Dataset, levels=unique(sets$Set)))

palette <- brewer.pal(9, "Set1")[-c(6:8)]

# plot of LOCATION of the 90% uncertainty bands for each subset
ggplot( ) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.4)+
  geom_ribbon(data=SLP.loc.g, aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=SLP.loc.ens.g, aes(x= dates, y=value, color=Dataset))+ scale_color_manual(values = c(palette[-6],"black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Position North Subtropical High [Max. SLP] Ens. Memb. - only oceans", y="Latitude [°]")+
  theme_bw()+theme(legend.position = c(0.2,0.5), legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="00"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

ggplot( ) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_vline(data= Volc %>% subset(., Date>"1750-01-01"), aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.5)+
  geom_ribbon(data=SLP.loc.g %>% subset(., dates>"1750-01-01"), aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=SLP.loc.ens.g %>% subset(., dates>"1750-01-01"), aes(x= dates, y=value, color=Dataset))+ scale_color_manual(values = c(palette[-c(6)],"black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="20 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Position North Subtropical High [Max. SLP] Ens. Memb. - only oceans", y="Latitude [°]")+
  theme_bw()+theme(legend.position = c(0.2,0.5), legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="00"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))


subThigh <- subThigh %>% 
  reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","Season","value","Memb","Dataset"))

SLP.str.g <-  subThigh %>% 
  within(., {  
    Dataset <- factor(Dataset, levels=unique(sets$Set))
    value <- value/100
    }) %>% 
  dplyr::group_by(., Dataset, dates, Season) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                   p95=quantile(value, probs= 0.95, na.rm = T),
                   min= min(value, na.rm = T),
                   max= max(value, na.rm = T)) %>% 
  as.data.frame()

SLP.str.ens.g <- subset(subThigh, Memb =="ensmean") %>% within(.,{
  value <- value/100
  Dataset <- factor(Dataset, levels=unique(sets$Set))
})

# plot of STRENGHT of the 90% uncertainty band for each subset
ggplot() +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_vline(data= Volc, aes(xintercept=Date, linetype=variable), col="#a65628", show.legend = FALSE, alpha=0.2)+
  geom_ribbon(data= SLP.str.g, aes(x= dates, fill=Dataset, ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=SLP.str.ens.g, aes(x= dates, y=value, color=Dataset))+ scale_color_manual(values = c(palette[-6],"black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
    labs(title="North Subtropical High [Max. SLP] Ens. Memb. - only oceans", y="SLP [hPa]")+
  theme_bw()+theme(legend.position = c(0.2,0.5), legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="00"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

################################-
## probability distribution of each sub-ensemble ----
################################-

SLP.str.g <- subThigh %>% 
  reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","Season","value","Memb","Dataset")) %>% 
  within(., {  
    Dataset <- factor(Dataset, levels=unique(sets$Set))
    value <- value*100
  }) %>% dplyr::group_by(., Dataset, dates, Season) %>%
  dplyr::reframe(SLP=density(value, na.rm=T)$x, pdf= density(value, na.rm=T)$y, 
                 normal= shapiro.test(value)$p.value >0.05)

spec.Dataset <- "set_1420-1"

# plot of the PDF with shades of red for each time step
SLP.str.g %>%
  subset(.,Season=="DJF" & Dataset== spec.Dataset) %>% 
  ggplot(., aes(x= dates, y= SLP, color=pdf)) +
  facet_wrap(. ~ Dataset, scales = "free_y",ncol=1)+
  
  geom_point(size=0.5)+
  scale_color_gradient(low="#ffffcc",high="#800026")+
  
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  scale_y_continuous(transform = "reverse")+
  labs(title="Strong Ascent [SLP 500 hPa] Ens. Memb. - only oceans", y="SLP 500 hPa [hPa/s]")+
  theme_bw()+theme(legend.position ="bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# plot of test of normal distribution for each subset
SLP.str.g %>%
  ggplot(.,aes(x=dates,y=Dataset,fill=normal))+
  facet_wrap(.~Season, ncol=1)+
  geom_tile()+
  scale_fill_manual(values=c("red","blue"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Normal dist. in Strong Ascent [SLP 500 hPa] Ens. Memb. - only oceans", y="Subset")+
  theme_bw()+theme(legend.position ="bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

################################-
## probability distribution of all the ensemble ----
################################-

ModE_Sim.normality <- subThigh %>% 
  reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","Season","value","Memb","Dataset")) %>% 
  within(., {  
    Dataset <- factor(Dataset, levels=unique(sets$Set))
    value <- value*100
  }) %>% subset(., Dataset!="ModE-RA") %>% 
  dplyr::group_by(., dates, Season) %>%
  dplyr::reframe(SLP=density(value, na.rm=T)$x, pdf= density(value, na.rm=T)$y, 
                 normal= shapiro.test(value)$p.value >0.05)

# plot of normality test
ModE_Sim.normality %>% 
  ggplot(.,aes(x=dates,y=Season,fill=normal))+
  geom_tile()+
  scale_fill_manual(values=c("red","blue"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  scale_y_discrete(limits=rev(levels(ModE_Sim.normality$Season)))+
  labs(title="Normal dist. in Strong Ascent [Min. SLP 500 hPa] All Ens. Memb. - only oceans", y="Subset")+
  theme_bw()+theme(legend.position ="bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# plot for the PDF in shades of red for each day considering the whole ModE-Sim ensemble
ModE_Sim.normality %>%
  ggplot(., aes(x= dates, y= SLP, color=pdf)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_point(size=0.1)+
  scale_color_gradient(low="#ffffcc",high="#800026")+
  
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  scale_y_continuous(transform = "reverse")+
  labs(title="Strong Ascent [SLP 500 hPa] All Ens. Memb.", y="SLP 500 hPa [hPa/s]")+
  theme_bw()+theme(legend.position ="bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
