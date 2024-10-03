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
        ) 

### Dates ----
Dates <- list()
Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="month")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="month")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="month")
Dates$all <- seq(as.Date("1420-01-01"), as.Date("2009-12-01"), by="year")
seasons <- c("DJF","JJA")

################################-
## load data ----
################################-
mastrf <<- list()
for (i in unique(sets$Set)) mastrf[[i]] <- list()

f <- nc_open("./01_Data/01_Streamfunction/ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_sDJF.nc")
lat <<- f$dim$lat$vals; lon <<- f$dim$lon$vals; lev <<- f$dim$plev$vals

read_Mod <- function (Set, Memb, Epoch, varn="mastrfu"){
  
  Dataset <- "ModE-Sim"

  if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
  
  f.mod <- paste0("./01_Data/01_Streamfunction/",Set,"/",Dataset,"_",Set,"_",Memb,"_",t.span,"_mastrf")
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  mastrf[[Set]][[Memb]] <<- ncvar_get(f, varid=varn)
  
  return(mastrf)
}

print("reading members")

pb <- txtProgressBar(min=1,max=nrow(sets),style=3)
for (i in 1:nrow(sets)){
  setTxtProgressBar(pb,i)
  
  a <- with(sets, read_Mod(Set= Set[i], Memb= Memb[i], Epoch= Epoch[i]))
}
close(pb)

print("reading ensemble means")
mastrf.ens <- list()
varn <- "mastrfu"
for (i in 1:nrow(sets.ens)){
  
  Dataset <- "ModE-Sim"
  f.mod <- with(sets.ens[i,],{
    
    if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
    f.mod <- paste0("./01_Data/01_Streamfunction/",Set,"/",Dataset,"_",Set,"_",Memb,"_mastrf_",t.span)
    
    f.mod
  })
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  mastrf.ens[[sets.ens$Set[i]]] <- ncvar_get(f, varid=varn)
  
}

# save(mastrf, file="./01_Data/01_Streamfunction/_allMemb_mastrf500hPa.Rdata")
################################-
# Mid tropospheric calculation ----
################################-
MidT <- function(set.m){
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  midt.strf <- set.m %>% mclapply(., function(x){
    y <- x[, lev >= 30000 & lev <= 70000,] %>% # select mid-troposphere values
      apply(., c(1,3), mean) %>%  # average the mid-troposphere
      t()
    return(y)
  })
  
  stopCluster(cl)
  return(midt.strf)
}


MidT.mastrf <- list()
MidT.mastrf.ens <- MidT(mastrf.ens) # for ensemble means
for(i in names(mastrf)){
  print(i)
  MidT.mastrf[[i]] <- MidT(mastrf[[i]]) # for members
  MidT.mastrf[[i]][["ensmean"]] <- MidT.mastrf.ens[[i]]
}

################################-
# DJF & JJA aggregation ----
################################-
mastrf.s <- list()

DJF_JJA <- function( mastrf, Dates.memb){
  
  Memb <- list()

  Memb$DJF <- mastrf %>% zoo(., order.by= Dates.memb) %>% dm2seasonal(., "DJF", FUN=mean, na.rm=T) %>% as.matrix()
  Memb$JJA <- mastrf %>% zoo(., order.by= Dates.memb) %>% dm2seasonal(., "JJA", FUN=mean, na.rm=T) %>% as.matrix()
  
  return(Memb)
}

cl <- makeCluster(10)
registerDoParallel(cl)

print("seasonal accumulation")
for ( i in names(mastrf)){
  print(i)
  
  Epoch <- with (sets, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)
  
  mastrf.s[[i]] <- mclapply(MidT.mastrf[[i]], DJF_JJA, Dates[[Epoch]])

}
stopCluster(cl)

# save(mastrf.s, file="./01_Data/01_Streamfunction/_allMemb_mastrf_DJFJJA.Rdata")

################################-
# ITCZ features ----
################################-
# load("./01_Data/01_Streamfunction/_allMemb_mastrf500hPa_DJFJJA.Rdata")

# the function for identifying the ITCZ feature' location, based on smoothing spline, is loaded with the settings file (00_settings.R)

# function for calculating the feature for both seasons with the smoothing spline
ITCZ.fea.memb <- function(y, Lat, Trop.fil){
  
  # filter to tropics
  for (j in seasons) y[[j]] <- y[[j]][, Trop.fil[,j]]
  
  Memb <- list()
  
  prop.strFU <- function(z){
    z.prop <- within(z,{ 
      width <- abs( lat.Max - lat.Min )
      Area <- abs( 2* pi *( 6371e3 )^2* (sin(lat.Max * pi/180) - sin(lat.Min * pi/180)) ) # radious of the earth
      strength <- abs( -9.8*(Max - Min)/Area ) # gravity
    })
    return(z.prop)
  }
  Memb[["DJF"]] <- data.frame(lat.Max = smt.min.max(y[["DJF"]], "max", Lat[Trop.fil[,"DJF"]]),
                              lat.Min = smt.min.max(y[["DJF"]], "min", Lat[Trop.fil[,"DJF"]]),
                              Max= apply(y[["DJF"]], 1, "max"),
                              Min= apply(y[["DJF"]], 1, "min")) %>% 
    prop.strFU(.)
  
  Memb[["JJA"]] <- data.frame(lat.Max = smt.min.max(y[["JJA"]], "max", Lat[Trop.fil[,"JJA"]]),
                              lat.Min = smt.min.max(y[["JJA"]], "min", Lat[Trop.fil[,"JJA"]]),
                              Max= apply(y[["JJA"]], 1, "max"),
                              Min= apply(y[["JJA"]], 1, "min")) %>% 
    prop.strFU(.)
  
    return(Memb)
}

Mastrf.prop <- list()
trop.fil <- data.frame(DJF= lat>=-35 & lat <= 30, JJA= lat>=-25 & lat <= 40)

Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")

print("identifying position")
for( i in names(mastrf.s)){ # over each subset
  print(i)
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  
  Mastrf.prop[[i]] <- mclapply(mastrf.s[[i]], ITCZ.fea.memb, lat, trop.fil)
  
  stopCluster(cl)
  
  Epoch <- with (sets, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)  # defines the Epoch for attaching the dates to the time series of the members
  
  Mastrf.prop[[i]] <- lapply(Mastrf.prop[[i]], 
                        function(x,dates.x){ 
                          lapply(x, function(y, dates) {
                          z <- cbind.data.frame(dates, y)
                          return(z)
                          }, dates.x)},
                        Dates[[Epoch]])
}


save(Mastrf.prop,file="./01_Data/01_Streamfunction/_allMemb_Prop_MastrFu.Rdata")

################################-
## plotting features timeseries ----
################################-
load("./01_Data/01_Streamfunction/_allMemb_Prop_MastrFu.Rdata")

Mastrf.prop <- Mastrf.prop %>% 
  reshape::melt(., id=c("dates")) %>% 
  magrittr::set_colnames(., c("dates","variable","value","Season","Memb","Dataset")) %>% 
  within(., variable <- factor(variable, levels = c("lat.Max","lat.Min","Max","Min","strength","width","Area")))

mastrf.loc.g <- Mastrf.prop %>% 
  subset(., Memb !="ensmean" & variable %in% c("strength","width","Area")) %>% 
  dplyr::group_by(., Dataset, dates, Season, variable) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                  p95=quantile(value, probs= 0.95, na.rm = T),
                  min= min(value, na.rm = T),
                  max= max(value, na.rm = T)) %>% 
  as.data.frame() %>% within(., Dataset <- factor(Dataset, levels=unique(sets$Set)))

mastrf.loc.ens.g <- Mastrf.prop %>% 
  subset(., Memb =="ensmean" & variable %in% c("strength","width","Area")) %>%
  within(., Dataset <- factor(Dataset, levels=unique(sets$Set)))

palette <- brewer.pal(9, "Set1")[-c(6:8)]

# plot of LOCATION of the 90% uncertainty bands for each subset
ggplot( ) +
  facet_grid(variable ~ Season, switch= "y", scales = "free_y")+
  geom_ribbon(data=mastrf.loc.g, aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=mastrf.loc.ens.g, aes(x= dates, y=value, color=Dataset))+ scale_color_manual(values = c(palette[-6],"black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="ITCZ properties [mid-tropospheric mass streanfunction] - ModE-Sim Ens. Memb.", y="Feature")+
  theme_bw()+theme(legend.position = "bottom", legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

ggplot( ) +
  facet_wrap(. ~ Season, scales = "free_y", ncol=1)+
  geom_ribbon(data=mastrf.loc.g %>% subset(., variable=="width"), aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=mastrf.loc.ens.g %>% subset(., variable=="width"), aes(x= dates, y=value, color=Dataset))+ scale_color_manual(values = c(palette[-6],"black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="ITCZ width - ModE-Sim Ens. Memb.", y="Width [°]")+
  theme_bw()+theme(legend.position = "bottom", legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))

ggplot( ) +
  facet_wrap(. ~ Season, scales = "free_y", ncol=1)+
  geom_ribbon(data=mastrf.loc.g %>% subset(., variable=="strength"), aes(x= dates, fill=Dataset,ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  geom_line(data=mastrf.loc.ens.g %>% subset(., variable=="strength"), aes(x= dates, y=value, color=Dataset))+ scale_color_manual(values = c(palette[-6],"black"))+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
               date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="ITCZ strength - ModE-Sim Ens. Memb.", y="stregnth [Pa/s]")+
  theme_bw()+theme(legend.position = "bottom", legend.direction = "horizontal",
                   panel.grid = element_line(linetype="dashed",color="lightgrey"),
                   strip.text = element_text(size=12),
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))
