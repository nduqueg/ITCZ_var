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

Set <- "lowres_20mem_Set_1420-3_1850-1"
Memb <- paste0("m0",seq(41,60, by=1))

################################-
## fetch Omega500 and transform to zonal mean ----
################################-

# calculate the zonal mean
setwd(paste0(dir.base,"./01_Data/02_Omega500/"))
numCores <- 20
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max= length(Memb), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


# zonal mean for each one of the members
################################-
results <- foreach ( i = 1:length(Memb), .combine=rbind, .options.snow=opts) %dopar% {
  
  if( !dir.exists(Set)) dir.create(Set)
  
  cdo.cmd <- paste0("cdo zonmean ",dir.data.ModERA,Set,"/abs/",Memb[i],"/by_var/mon/ModE-RA_",Set,"_",Memb[i],"_omega_50000_abs_1421-2008_mon.nc ",
                                   "./",Set,"/ModE-RA_",Set,"_",Memb[i],"_Omega500hPa-ZonMean_1421-2008_mon.nc")
  system(cdo.cmd)
}
stopCluster(cl)

################################-
# zonal mean for the general ensemble ----
################################-

cdo.cmd <- paste0("cdo zonmean ",dir.data.ModERA,Set,"/abs/ensstat/by_var/mon/ModE-RA_lowres_20mem_Set_1420-3_1850-1_ensmean_omega_50000_abs_1421-2008_mon.nc ",
                  "ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_mon.nc")
system(cdo.cmd)

# aggregations and selection  
################################-

cdo.cmd <- "cdo seasmean ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_mon.nc ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_seasonal.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,DJF ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_seasonal.nc ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_sDJF.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,JJA ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_seasonal.nc ModE-RA_lowres_20mem_Set_1420-3_1850-1_Omega500hPa-ZonMean_1421-2008_sJJA.nc"
system(cdo.cmd)

setwd("../..")


