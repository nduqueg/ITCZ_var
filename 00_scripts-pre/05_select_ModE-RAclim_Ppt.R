rm(list = ls())
cat("\014")

library(parallel)
library(doSNOW)
library(tcltk)
"%>%"=magrittr::`%>%`

source("00_settings.R")
# dir.data.ModERA <- "/mnt/climstor/ERC_PALAEO/ModE-RA/outdata/"

################################-
## generate names in ModE-RA data ----
################################-

# names of directories
# the indiviual members and the ensamble of Omega 500 hPa are already calculated and separated in
# lowres_20mem_Set_1420-3_1850-1/abs/ensstat/by_var/mon
# lowres_20mem_Set_1420-3_1850-1/abs/m001/by_var/mon

Set <- "lowres_100mem_Set_1"
abs <- "ModE-RAclim_anom_after_assim_71yr_highpassfiltered_observations"
Memb <- character()
for ( i in 1:100){ if (i>=10 & i <100) Memb[i] <- paste0("m0",i) else if (i < 10) Memb[i] <- paste0("m00",i) else Memb[i] <- paste0("m",i)}

################################-
## fetch Ppt and transform to zonal mean ----
################################-

# calculate the zonal mean
setwd(paste0(dir.base,"./01_Data/05_Ppt/"))
numCores <- 10
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max= length(Memb), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


# zonal mean for each one of the members
################################-
results <- foreach ( i = 1:length(Memb), .combine=rbind, .options.snow=opts) %dopar% {
  
  if( !dir.exists(Set)) dir.create(Set)
  
  cdo.cmd <- paste0("cdo chunit,'kg m-2 s-1','mm/month' -mulc,86400 -muldpm -zonmean ",
                    dir.data.ModERAclim, Set,"/",abs,"/",Memb[i],"/by_var/mon/ModE-RAclim_",Memb[i],"_totprec_anom_1421-2008_mon.nc ",
                    "./",Set,"/ModE-RAclim_",Memb[i],"_totprec-ZonMean_1421-2008_mon.nc")
  system(cdo.cmd)
}
stopCluster(cl)

################################-
# zonal mean for the general ensemble ----
################################-

cdo.cmd <- paste0("cdo chunit,'kg m-2 s-1','mm/month' -mulc,86400 -muldpm -zonmean ",
                  dir.data.ModERAclim,Set,"/",abs,"/ensstat/by_var/mon/ModE-RAclim_ensmean_totprec_anom_1421-2008_mon.nc ",
                  "ModE-RAclim_",Set,"_totprec-ZonMean_1421-2008_mon.nc")
system(cdo.cmd)

setwd("../..")
