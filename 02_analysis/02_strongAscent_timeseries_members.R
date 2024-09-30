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
Omega <<- list() 
for(i in unique(sets$Set)) Omega[[i]] <- list()

Omega.ens <- list()

f <- nc_open("./01_Data/02_Omega500/ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_sDJF.nc")
lat <<- f$dim$lat$vals; lon <<- f$dim$lon$vals; lev <<- f$dim$plev$vals


read_Mod <- function (Set, Memb, Epoch, varn="omega"){
  
  Dataset <- "ModE-Sim"
  
  if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
  
  f.mod <- paste0("./01_Data/02_Omega500/",Set,"/",Dataset,"_",Set,"_",Memb,"_Omega500hPa-ZonMean_",t.span)
  
  if (Set == "ModE-RA"){
    Dataset <- "ModE-RA"
    t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/02_Omega500/lowres_20mem_Set_1420-3_1850-1/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_",Memb,"_Omega500hPa-ZonMean_",t.span)
  }
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  
  Omega[[Set]][[Memb]] <<- ncvar_get(f, varid=varn)

  # return(Omega)
}

print("reading members")

pb <- txtProgressBar(min=1,max=nrow(sets),style=3)
for (i in 1:nrow(sets)){
  setTxtProgressBar(pb,i)
  
  a <- with(sets, read_Mod(Set= Set[i], Memb= Memb[i], Epoch= Epoch[i]))
}
close(pb)

print("reading ensemble means")
varn <- "omega"
for (i in 1:nrow(sets.ens)){
  
  Dataset <- "ModE-Sim"
  f <- with(sets.ens[i,],{
    
    if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
    f.mod <- paste0("./01_Data/02_Omega500/",Set,"/",Dataset,"_",Set,"_",Memb,"_Omega500hPa-ZonMean_",t.span)
    if (Set == "ModE-RA"){
      Dataset <- "ModE-RA"
      t.span <- "1421-2008"
      f.mod <- paste0("./01_Data/02_Omega500/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_",t.span)
    }
    
    paste0(f.mod,"_mon.nc") %>% nc_open(.)
  })
  Omega.ens[[sets.ens$Set[i]]] <- ncvar_get(f, varid=varn)
}

# save(Omega, file="./01_Data/02_Omega500/_allMemb_Omega500hPa.Rdata")

################################-
# DJF & JJA aggregation ----
################################-
Omega.s <- list()
Omega.ens.s <- list()

DJF_JJA <- function( Omega, Dates.memb){
  
  Memb <- list()

  Memb$DJF <- t(Omega) %>% zoo(., order.by= Dates.memb) %>% dm2seasonal(., "DJF", FUN=mean, na.rm=T) %>% as.matrix()
  Memb$JJA <- t(Omega) %>% zoo(., order.by= Dates.memb) %>% dm2seasonal(., "JJA", FUN=mean, na.rm=T) %>% as.matrix()
  
  return(Memb)
}

cl <- makeCluster(10)
registerDoParallel(cl)

print("seasonal accumulation")
for ( i in names(Omega)){
  print(i)
  
  Epoch <- with (sets, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)
  
  Omega.s[[i]] <- mclapply(Omega[[i]], DJF_JJA, Dates[[Epoch]])

}
stopCluster(cl)

print("seasonal accumulation Ensemble Means")
for ( i in names(Omega.ens)){
  print(i)
  
  Epoch <- with (sets.ens, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)
  
  Omega.ens.s[[i]] <- DJF_JJA(Omega.ens[[i]], Dates[[Epoch]])
}
# save(Omega.s, file="./01_Data/02_Omega500/_allMemb_Omega500hPa_DJFJJA.Rdata")

################################-
# Strong ascent ----
################################-
# load("./01_Data/02_Omega500/_allMemb_Omega500hPa_DJFJJA.Rdata")

# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

# function for calculating the feature for both seasons with the smoothing spline
min.max.memb <- function(y, Lat, feature ="loc", Trop.fil){
  
  # filter to tropics
  for (j in seasons) y[[j]] <- y[[j]][, Trop.fil[,j]]
  Lat.trop.fil <- Lat[Trop.fil]
  
  Memb <- list()
  if (feature =="loc"){
    Memb[["DJF"]] <- smt.min.max(y[["DJF"]], "min", Lat.trop.fil)
    Memb[["JJA"]] <- smt.min.max(y[["JJA"]], "min", Lat.trop.fil)  
  } else if(feature == "Str"){
    Memb[["DJF"]] <- apply(y[["DJF"]], 1, min)
    Memb[["JJA"]] <- apply(y[["JJA"]], 1, min)
  }
  
  
  return(Memb)
}

loc.strAsc <- list(); loc.strAsc.ens <- list()
strAsc <- list(); strAsc.ens <- list()

trop.fil <- lat >= -45 & lat <= 45

Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")

print("identifying position")
for( i in names(Omega.s)){ # two time periods
  print(i)
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  loc.strAsc[[i]] <- mclapply(Omega.s[[i]], min.max.memb, lat, "loc", trop.fil)
  strAsc[[i]] <- mclapply(Omega.s[[i]], min.max.memb, lat, "Str", trop.fil)
  stopCluster(cl)
  
  Epoch <- with (sets, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)  # defines the Epoch for attaching the dates to the time series of the members
  
  loc.strAsc[[i]] <- lapply(loc.strAsc[[i]], 
                        function(x,dates){ y <- as.data.frame(x); z <- cbind.data.frame(dates,y); return(z)},
                        Dates[[Epoch]])
  strAsc[[i]] <- lapply(strAsc[[i]], 
                            function(x,dates){ y <- as.data.frame(x); z <- cbind.data.frame(dates,y); return(z)},
                            Dates[[Epoch]])
}

print("identifying position Ensemble means")
for( i in names(Omega.ens.s)){ # two time periods
  print(i)

  Epoch <- with (sets.ens, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)  # defines the Epoch for attaching the dates to the time series of the members
  dates <- Dates[[Epoch]]
  
  loc.strAsc[[i]][["ensmean"]] <- min.max.memb(Omega.ens.s[[i]], lat, "loc") %>% cbind.data.frame(dates, .)
  strAsc[[i]][["ensmean"]] <- min.max.memb(Omega.ens.s[[i]], lat, "Str") %>% cbind.data.frame(dates, .)
  
}

save(loc.strAsc,file="./01_Data/02_Omega500/_allMemb_LocStrAsc.Rdata")
save(strAsc,file="./01_Data/02_Omega500/_allMemb_StrAsc.Rdata")

################################-
## plotting features timeseries ----
################################-
load("./01_Data/02_Omega500/_allMemb_LocStrAsc.Rdata")
load("./01_Data/02_Omega500/_allMemb_StrAsc.Rdata")

loc.strAsc <- loc.strAsc %>% 
  reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","Season","value","Memb","Dataset"))

Omega.loc.g <- subset(loc.strAsc, Memb !="ensmean") %>% 
  dplyr::group_by(., Dataset, dates, Season) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                  p95=quantile(value, probs= 0.95, na.rm = T),
                  min= min(value, na.rm = T),
                  max= max(value, na.rm = T)) %>% 
  as.data.frame() %>% within(., Dataset <- factor(Dataset, levels=unique(sets$Set)))

Omega.loc.ens.g <- subset(loc.strAsc, Memb =="ensmean") %>% within(., Dataset <- factor(Dataset, levels=unique(sets$Set)))

palette <- brewer.pal(9, "Set1")[-c(6:8)]

# plot of LOCATION of the 90% uncertainty bands for each subset
ggplot( ) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_ribbon(data=Omega.loc.g, aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=Omega.loc.ens.g, aes(x= dates, y=value, color=Dataset))+ scale_color_manual(values = c(palette[-6],"black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Position Strong Ascent [Min. Omega 500 hPa] Ens. Memb.", y="Latitude [°]")+
  theme_bw()+theme(legend.position = c(0.2,0.1), legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

strAsc <- strAsc %>% 
  reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","Season","value","Memb","Dataset"))

Omega.str.g <-  strAsc %>% 
  within(., {  
    Dataset <- factor(Dataset, levels=unique(sets$Set))
    value <- value*100
    }) %>% 
  dplyr::group_by(., Dataset, dates, Season) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                   p95=quantile(value, probs= 0.95, na.rm = T),
                   min= min(value, na.rm = T),
                   max= max(value, na.rm = T)) %>% 
  as.data.frame()

Omega.str.ens.g <- subset(strAsc, Memb =="ensmean") %>% within(.,{
  value <- value*100
  Dataset <- factor(Dataset, levels=unique(sets$Set))
})

# plot of STRENGHT of the 90% uncertainty band for each subset
ggplot() +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_ribbon(data= Omega.str.g, aes(x= dates, fill=Dataset, ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=Omega.str.ens.g, aes(x= dates, y=value, color=Dataset))+ scale_color_manual(values = c(palette[-6],"black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  scale_y_continuous(transform = "reverse")+
  labs(title="Strong Ascent [Min. Omega 500 hPa] Ens. Memb.", y="Omega 500 hPa [hPa/s]")+
  theme_bw()+theme(legend.position = c(0.2,0.5), legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

################################-
## probability distribution of each sub-ensemble ----
################################-

Omega.str.g <- strAsc %>% 
  reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","Season","value","Memb","Dataset")) %>% 
  within(., {  
    Dataset <- factor(Dataset, levels=unique(sets$Set))
    value <- value*100
  }) %>% dplyr::group_by(., Dataset, dates, Season) %>%
  dplyr::reframe(omega=density(value, na.rm=T)$x, pdf= density(value, na.rm=T)$y, 
                 normal= shapiro.test(value)$p.value >0.05)

spec.Dataset <- "set_1420-1"

# plot of the PDF with shades of red for each time step
Omega.str.g %>%
  subset(.,Season=="DJF" & Dataset== spec.Dataset) %>% 
  ggplot(., aes(x= dates, y= omega, color=pdf)) +
  facet_wrap(. ~ Dataset, scales = "free_y",ncol=1)+
  
  geom_point(size=0.5)+
  scale_color_gradient(low="#ffffcc",high="#800026")+
  
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  scale_y_continuous(transform = "reverse")+
  labs(title="Strong Ascent [Omega 500 hPa] Ens. Memb.", y="Omega 500 hPa [hPa/s]")+
  theme_bw()+theme(legend.position ="bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# plot of test of normal distribution for each subset
Omega.str.g %>%
  ggplot(.,aes(x=dates,y=Dataset,fill=normal))+
  facet_wrap(.~Season, ncol=1)+
  geom_tile()+
  scale_fill_manual(values=c("red","blue"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Normal dist. in Strong Ascent [Omega 500 hPa] Ens. Memb.", y="Subset")+
  theme_bw()+theme(legend.position ="bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

################################-
## probability distribution of all the ensemble ----
################################-

ModE_Sim.normality <- strAsc %>% 
  reshape::melt(., id=c("dates")) %>% magrittr::set_colnames(., c("dates","Season","value","Memb","Dataset")) %>% 
  within(., {  
    Dataset <- factor(Dataset, levels=unique(sets$Set))
    value <- value*100
  }) %>% subset(., Dataset!="ModE-RA") %>% 
  dplyr::group_by(., dates, Season) %>%
  dplyr::reframe(omega=density(value, na.rm=T)$x, pdf= density(value, na.rm=T)$y, 
                 normal= shapiro.test(value)$p.value >0.05)

# plot of normality test
ModE_Sim.normality %>% 
  ggplot(.,aes(x=dates,y=Season,fill=normal))+
  geom_tile()+
  scale_fill_manual(values=c("red","blue"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  scale_y_discrete(limits=rev(levels(ModE_Sim.normality$Season)))+
  labs(title="Normal dist. in Strong Ascent [Min. Omega 500 hPa] All Ens. Memb.", y="Subset")+
  theme_bw()+theme(legend.position ="bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

# plot for the PDF in shades of red for each day considering the whole ModE-Sim ensemble
ModE_Sim.normality %>%
  ggplot(., aes(x= dates, y= omega, color=pdf)) +
  facet_wrap(. ~ Season, scales = "free_y",ncol=1)+
  geom_point(size=0.1)+
  scale_color_gradient(low="#ffffcc",high="#800026")+
  
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  scale_y_continuous(transform = "reverse")+
  labs(title="Strong Ascent [Omega 500 hPa] All Ens. Memb.", y="Omega 500 hPa [hPa/s]")+
  theme_bw()+theme(legend.position ="bottom", legend.direction = "horizontal",legend.background = element_rect(color = "black"),
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
