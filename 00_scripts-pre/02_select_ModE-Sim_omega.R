rm(list = ls())
cat("\014")

library(parallel)
library(doSNOW)
library(tcltk)
"%>%"=magrittr::`%>%`

source("00_settings.R")
# dir.base.ModEsim <- "/mnt/climstor/ERC_PALAEO/ModE-Sim/outdata/"

## generate names in ModE-Sim data ----

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

## fetch Omega500 and transform to zonal mean ----

# the indiviual members and the ensamble of Omega 500 hPa are already calculated and separated in
# set_1420-1_to_3/abs/ensstat/by_var/mon
# set_1420-1/abs/m001/by_var/mon

# calculate the zonal mean
setwd(paste0(dir.base,"./01_Data/02_Omega500/"))
numCores <- 5
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max= nrow(sets), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

results <- foreach ( i = 1:nrow(sets), .combine=rbind, .options.snow=opts) %dopar% {
  
  if( !dir.exists(sets$Set[i])) dir.create(sets$Set[i])
  
  
  iEpoch <- sets[i,"Epoch"]
  Per <- Epoch[[iEpoch]] %>% paste0(.,"-",tail(.,1)) %>% .[1]
  
  cdo.cmd <- with(sets[i,], paste0("cdo zonmean ",dir.base,Set,"/abs/",Memb,"/by_var/mon/ModE-Sim_",Set,"_",Memb,"_omega_50000_abs_",Per,"_mon.nc ",
                                   "./",sets$Set[i],"/ModE-Sim_",Set,"_",Memb,"_Omega500hPa-ZonMean_",Per,"_mon.nc"))
  system(cdo.cmd)
}
stopCluster(cl)


cdo.cmd <- paste0("cdo zonmean ",dir.base,"set_1420-1_to_3/abs/ensstat/by_var/mon/ModE-Sim_set_1420-1_to_3_ensmean_omega_50000_abs_1420-1849_mon.nc ",
                  "ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_mon.nc")
system(cdo.cmd)

cdo.cmd <- paste0("cdo zonmean ",dir.base,"set_1850-1_to_2/abs/ensstat/by_var/mon/ModE-Sim_set_1850-1_to_2_ensmean_omega_50000_abs_1850-2009_mon.nc ",
                  "ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_mon.nc")
system(cdo.cmd)

cdo.cmd <- "cdo seasmean ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_mon.nc ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_seasonal.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,DJF ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_seasonal.nc ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_sDJF.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,JJA ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_seasonal.nc ModE-Sim_set_1420-1_to_3_Omega500hPa-ZonMean_1420-1849_sJJA.nc"
system(cdo.cmd)

cdo.cmd <- "cdo seasmean ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_mon.nc ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_seasonal.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,DJF ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_seasonal.nc ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_sDJF.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,JJA ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_seasonal.nc ModE-Sim_set_1850-1_to_2_Omega500hPa-ZonMean_1850-2009_sJJA.nc"
system(cdo.cmd)

setwd(dir.base)
