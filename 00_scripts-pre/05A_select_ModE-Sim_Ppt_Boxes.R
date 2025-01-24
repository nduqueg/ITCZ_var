rm(list = ls())
cat("\014")

library(parallel)
library(doSNOW)
library(tcltk)
library(ncdf4)
"%>%"=magrittr::`%>%`

source("00_settings.R")
# dir.base.ModEsim <- "/mnt/climstor/ERC_PALAEO/ModE-Sim/outdata/"

################################-
## generate names in ModE-Sim data and define boxes ----
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
        )

# time periods for each member
Epoch <- list(Epoch1 = seq(1420,1849), Epoch2= seq(1850,2009))

# boxes 
boxes <- data.frame(name=c("NAmMon","SAmMon","Sahel","EAfrMon","SEAsiaMon","NAusMon","NTrPac","Galapagos"),
                    lon.min=c(-100, -70,     -20,     15,       85,        120,      -160,    -95),
                    lon.max=c(-80,  -50,      30,     40,      120,        160,      -120,    -85),
                    lat.min=c(10,   -20,       5,    -25,        5,        -20,       2.5,    -1),
                    lat.max=c(20,    -5,      15,     -5,       25,         -5,        15,     3))

################################-
## fetch Omega500 and transform to zonal mean ----
################################-

# the indiviual members and the ensamble of Omega 500 hPa are already calculated and separated in
# set_1420-1_to_3/abs/ensstat/by_var/mon
# set_1420-1/abs/m001/by_var/mon

# calculate the zonal mean
setwd(paste0(dir.base,"./01_Data/05A_Ppt_boxes/"))
numCores <- 5
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max= nrow(sets), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


# box averages for each one of the members
################################-
results <- foreach ( i = 1:nrow(sets), .combine=rbind, .options.snow=opts, .packages="ncdf4") %dopar% {
  
  if( !dir.exists(sets$Set[i])) dir.create(sets$Set[i])
  
  iEpoch <- sets[i,"Epoch"]
  Per <- Epoch[[iEpoch]] %>% paste0(.,"-",tail(.,1)) %>% .[1]
  
  P.memb <- list()
  for (j in boxes$name){
    boundaries <- subset(boxes, name==j) %>% with(., paste0(lon.min,",",lon.max,",",lat.min,",",lat.max))
    cdo.cmd <- with(sets[i,],
                    paste0("cdo muldpm -fldmean -sellonlatbox,", boundaries," ",
                           dir.data.ModEsim,Set,"/abs/",Memb,"/by_var/mon/ModE-Sim_",Set,"_",Memb,"_totprec_abs_",Per,"_mon.nc ",
                           "./",sets$Set[i],"/ModE-Sim_",Set,"_",Memb,"_totprec-",j,"_",Per,"_mon.nc"))
    system(cdo.cmd)
    
    f <- with(sets[i,],
              nc_open(paste0("./",sets$Set[i],"/ModE-Sim_",Set,"_",Memb,"_totprec-",j,"_",Per,"_mon.nc")))
    
    P.memb[[j]] <- ncvar_get(f, varid="totprec") * 86400 # per seconds to per day
    
    with(sets[i,],
         file.remove(paste0("./",sets$Set[i],"/ModE-Sim_",Set,"_",Memb,"_totprec-",j,"_",Per,"_mon.nc")))
  }
  
  Dates <-  seq(
    substr(Per,1,4) %>%  paste0("-01-01") %>% as.Date(),
    substr(Per,6,9)  %>%  paste0("-12-31") %>% as.Date(),
    by="month")
  
  P.memb <- as.data.frame(P.memb) %>% cbind(Dates,.)
  with(sets[i,],paste0("./",sets$Set[i],"/ModE-Sim_",Set,"_",Memb,"_totprec-Boxes_",Per,"_mon.RData")) %>% 
    save(P.memb, file= .)
}
stopCluster(cl)

################################-
# box average for the general ensemble ----
################################-

set.ens <- c("set_1420-1_to_3","set_1850-1_to_2")
Per <- data.frame(per=c("1420-1849","1850-2009"),row.names = set.ens)
for(i in set.ens){
  
  P.memb <- list()
  for (j in boxes$name){
    boundaries <- subset(boxes, name==j) %>% with(., paste0(lon.min,",",lon.max,",",lat.min,",",lat.max))
    cdo.cmd <- paste0("cdo muldpm -fldmean -sellonlatbox,", boundaries," ",
                      dir.data.ModEsim, i,"/abs/ensstat/by_var/mon/ModE-Sim_", i,"_ensmean_totprec_abs_",Per[i,],"_mon.nc ",
                      "aux.nc")
    system(cdo.cmd)
    
    f <- nc_open("aux.nc")
    P.memb[[j]] <- ncvar_get(f, varid="totprec") * 86400 # per seconds to per day
    
    file.remove("aux.nc")
  }
  
  Dates <- seq(as.Date(paste0(substr(Per[i,],1,4),"-01-01")), as.Date(paste0(substr(Per[i,],6,9),"-12-31")), by="month")
  P.memb <- as.data.frame(P.memb) %>% cbind(Dates,.)
  save(P.memb, file=paste0("ModE-Sim_", i,"_totprec-Boxes_",Per[i,],"_mon.RData"))
}

################################-
# box average for each set ensemble ----
################################-
sets.names <- unique(sets$Set)
epoch.lim <- c(rep("1420-1849",3), rep("1850-2009",2))

for ( i in sets.names){
  epoch.years <- which(i == sets.names) %>% epoch.lim[.]
  
  P.memb <- list()
  for (j in boxes$name){
    boundaries <- subset(boxes, name==j) %>% with(., paste0(lon.min,",",lon.max,",",lat.min,",",lat.max))
    cdo.cmd <- paste0("cdo muldpm -fldmean -sellonlatbox,", boundaries," ",
                      dir.data.ModEsim, i,"/abs/ensstat/by_var/mon/ModE-Sim_",i,"_ensmean_totprec_abs_",epoch.years,"_mon.nc ",
                      "aux.nc")
    system(cdo.cmd)
    
    f <- nc_open("aux.nc")
    P.memb[[j]] <- ncvar_get(f, varid="totprec") * 86400 # per seconds to per day
    
    file.remove("aux.nc")
  }
  
  Dates <- seq(as.Date(paste0(substr(epoch.years,1,4),"-01-01")), as.Date(paste0(substr(epoch.years,6,9),"-12-31")), by="month")
  P.memb <- as.data.frame(P.memb) %>% cbind(Dates,.)
  save(P.memb, file=paste0("./",i,"/ModE-Sim_",i,"_ensmean_totprec-Boxes_",epoch.years,"_mon.RData"))
}


setwd(dir.base)
