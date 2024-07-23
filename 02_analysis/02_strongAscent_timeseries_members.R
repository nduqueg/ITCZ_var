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
                       Epoch= 2
  )
  ) %>% 
  rbind(., data.frame( Set= "set_1850-2" %>% rep(., each=16),
                       Memb= paste0("m0",seq(21,36)),
                       Epoch= 2
  )
) %>% rbind(., data.frame( Set= "ModE-RA",
                           Memb= paste0("m0",seq(41,60)),
                           Epoch= 2
)
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

f <- nc_open("./01_Data/02_Omega500/ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_sDJF.nc")
lat <<- f$dim$lat$vals; lon <<- f$dim$lon$vals; lev <<- f$dim$plev$vals


read_Mod <- function (Set, Memb, Epoch, varn="omega"){
  
  Dataset <- "ModE-Sim"
  Omega[[Set]] <- list()
  if(Epoch == 1) {t.span <- "1420-1849"} else if(Epoch == 2){ t.span <- "1850-2009"}
  
  f.mod <- paste0("./01_Data/02_Omega500/",Set,"/",Dataset,"_",Set,"_",Memb,"_Omega500hPa-ZonMean_",t.span)
  
  if (Set == "ModE-RA"){
    Dataset <- "ModE-RA"
    t.span <- "1421-2008"
    f.mod <- paste0("./01_Data/02_Omega500/lowres_20mem_Set_1420-3_1850-1/",Dataset,"_lowres_20mem_Set_1420-3_1850-1_",Memb,"_Omega500hPa-ZonMean_",t.span)
  }
  
  f <- paste0(f.mod,"_mon.nc") %>% nc_open(.)
  
  Omega[[Set]][[Memb]] <<- ncvar_get(f, varid=varn)

  return(Omega)
}

print("reading members")

pb <- txtProgressBar(min=1,max=nrow(sets),style=3)
for (i in 1:nrow(sets)){
  setTxtProgressBar(pb,i)
  
  a <- with(sets, read_Mod(Set= Set[i], Memb= Memb[i], Epoch= Epoch[i]))
}
close(pb)

# save(Omega, file="./01_Data/02_Omega500/_allMemb_Omega500hPa.Rdata")

################################-
# DJF & JJA ----
################################-
Omega.s <- list()

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
# save(Omega.s, file="./01_Data/02_Omega500/_allMemb_Omega500hPa_DJFJJA.Rdata")

################################-
# Strong ascent ----
################################-
load("./01_Data/02_Omega500/_allMemb_Omega500hPa_DJFJJA.Rdata")
smt.min.max <- function(x,type, lat){
  # 'x' is the Omega at 500 hPa, which will be smoothed around the minimum with an spline curve and 
  # the position of its minimum will be also identified
  lat.f <- lat[lat >= -45 & lat <= 45]
  
  aprox <- x[,lat >= -45 & lat <= 45] %>% apply(.,1, paste0("which.",type) %>% get()) %>% as.numeric()
  
  ind.aprox <- cbind(aprox-3,aprox-2,aprox-1,aprox,aprox+1,aprox+2,aprox+3) %>% t() %>% as.data.frame() # indices of values close to the Max StrFunct for the Spline model
  lat.models <- apply(ind.aprox, 2, function(ind,lat) lat[ind], lat.f) %>%  as.data.frame() # latitudes for the Spline model
  
  data.spModel <- mapply(FUN= function(x,ind) x[ind], # extracts the MAx StrFunct values near the Max for the Spline model
                         x[,lat >= -45 & lat <= 45] %>% as.matrix() %>% split(.,row(.)), 
                         ind.aprox,
                         SIMPLIFY = F) %>% 
    
    mapply(function(y,lat) cbind.data.frame(lat,y), # paste for each time step the latitude and the values near the Max StrFunc
           .,
           lat.models, 
           SIMPLIFY = F) %>% 
    lapply(., `colnames<-`,c("x.mod","y.mod"))
  
  sp.models <- lapply(data.spModel, function(x) tryCatch( # there are some months that don't have values
    {lm(y.mod ~ splines::ns(x.mod,3), data= x)},
    error= function(cond){NA}
  )
  )
  pos.mastrf <- numeric()
  for(k in 1:length(sp.models)){ # for each year identify the position of the Strong Ascent
    
    if (is.list(sp.models[[k]])){
      
    pos.mastrf[k] <- predict(sp.models[[k]], newdata = data.frame(x.mod= seq(lat.models[7,k],lat.models[1,k],by=0.1))) %>% as.numeric(.) %>% matrix(.,nrow=1) %>% 
      apply(.,1,get(paste0("which.",type))) %>% 
      seq(lat.models[7,k],lat.models[1,k],by=0.1)[.]
    
    } else{
      pos.mastrf[k] <- NA
    }
  }
  return(pos.mastrf)
}

min.max.memb <- function(y, lat){
  
  Memb <- list()
  Memb[["DJF"]] <- smt.min.max(y[["DJF"]], "min", lat)
  Memb[["JJA"]] <- smt.min.max(y[["JJA"]], "min", lat)
  
  return(Memb)
}

strAsc <- list()
Dates$ep1 <- seq(as.Date("1420-01-01"), as.Date("1849-12-01"), by="year")
Dates$ep2 <- seq(as.Date("1850-01-01"), as.Date("2009-12-01"), by="year")
Dates$ModERA <- seq(as.Date("1421-01-01"), as.Date("2008-12-01"), by="year")

print("identifying position")
for( i in names(Omega.s)){ # two time periods
  print(i)
  
  cl <- makeCluster(20)
  registerDoParallel(cl)
  strAsc[[i]] <- mclapply(Omega.s[[i]], min.max.memb, lat)
  stopCluster(cl)
  
  Epoch <- with (sets, Epoch[i==Set][1])
  if (i=="ModE-RA") Epoch <- "ModERA" else Epoch <- paste0("ep",Epoch)
  
  strAsc[[i]] <- lapply(strAsc[[i]], 
                        function(x,dates){ y <- as.data.frame(x); z <- cbind.data.frame(dates,y);return(z)},
                        Dates[[Epoch]])
}

save(strAsc,file="./01_Data/02_Omega500/_allMemb_position.Rdata")

################################-
## plotting features timeseries ----
################################-

Omega.prop.g <- strAsc %>% reshape::melt(., id=c("dates")) %>% 
  dplyr::group_by(., L1, dates, variable) %>%
  dplyr::summarise(p5= quantile(value, probs = 0.05, na.rm = T),
                  p95=quantile(value, probs= 0.95, na.rm = T),
                  min= min(value, na.rm = T),
                  max= max(value, na.rm = T)) %>% 
  as.data.frame()

palette <- brewer.pal(9, "Set1")[-c(6:8)]
Omega.prop.g %>%
  ggplot(., aes(x= dates, fill=L1)) +
  facet_wrap(. ~ variable, scales = "free_y",ncol=1)+
  geom_ribbon(aes(ymin=p5,ymax=p95), alpha=0.3)+ scale_fill_manual(values = palette)+
  scale_x_date(breaks = seq(as.Date("1450-01-01"),as.Date("2000-01-01"),by="50 years"), 
                            date_labels = "%Y", expand=c(0.01,0.01))+
  labs(title="Position Strong Ascent [Omega 500 hPa] Ens. Memb.", y="Latitude [°]")+
  theme_bw()+theme(legend.position = c(0.2,0.1), legend.direction = "horizontal",
                   axis.ticks.length=unit(-4, "pt"), axis.text.x = element_text(margin=margin(2,5,5,5),vjust = -1, size=12), axis.text.y = element_text(margin=margin(0,5,5,0,"pt"),size=12))


