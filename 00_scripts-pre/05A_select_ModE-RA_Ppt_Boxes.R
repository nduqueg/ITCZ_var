rm(list = ls())
cat("\014")

library(parallel)
library(doSNOW)
library(tcltk)
library(ncdf4)
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

# boxes 
boxes <- data.frame(name=c("NAmMon","SAmMon","Sahel","EAfrMon","SEAsiaMon","NAusMon","NTrPac","Galapagos"),
                    lon.min=c(-100, -70,     -20,     15,       85,        120,      -160,    -95),
                    lon.max=c(-80,  -50,      30,     40,      120,        160,      -120,    -85),
                    lat.min=c(10,   -20,       5,    -25,        5,        -20,       2.5,    -1),
                    lat.max=c(20,    -5,      15,     -5,       25,         -5,        15,     3))

################################-
## fetch Omega500 and transform to zonal mean ----
################################-

# calculate the zonal mean
setwd(paste0(dir.base,"./01_Data/05A_Ppt_boxes/"))
numCores <- 10
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max= length(Memb), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


# boxes area-average for each one of the members
################################-
results <- foreach ( i = 1:length(Memb), .combine=rbind, .options.snow=opts, .packages="ncdf4") %dopar% {
  
  if( !dir.exists(Set)) dir.create(Set)
  
  P.memb <- list( )
  
  for ( j in boxes$name){
    boundaries <- subset(boxes, name==j) %>% with(., paste0(lon.min,",",lon.max,",",lat.min,",",lat.max))
    cdo.cmd <- paste0("cdo muldpm -fldmean -sellonlatbox,", boundaries," ",
                      dir.data.ModERA,Set,"/abs/",Memb[i],"/by_var/mon/ModE-RA_",Set,"_",Memb[i],"_totprec_abs_1421-2008_mon.nc ",
                      "./",Set,"/ModE-RA_",Set,"_",Memb[i],"_totprec-",j,"_1421-2008_mon.nc")
    system(cdo.cmd)
    
    f <- nc_open(paste0("./",Set,"/ModE-RA_",Set,"_",Memb[i],"_totprec-",j,"_1421-2008_mon.nc"))
    P.memb[[j]] <- ncvar_get(f, varid="totprec") * 86400 # per seconds to per day
    
    file.remove(paste0("./",Set,"/ModE-RA_",Set,"_",Memb[i],"_totprec-",j,"_1421-2008_mon.nc"))
  }
  
  Dates <- seq(as.Date("1421-01-01"), as.Date("2008-12-31"), by="month")
  P.memb <- as.data.frame(P.memb) %>% cbind(Dates,.)
  
  paste0("./",Set,"/ModE-RA_",Set,"_",Memb[i],"_totprec-Boxes_1421-2008_mon.RData") %>% 
    save(P.memb, file= .)
}
stopCluster(cl)

################################-
# box average for the general ensemble ----
################################-

P.memb <- list( )

for ( j in boxes$name){
  boundaries <- subset(boxes, name==j) %>% with(., paste0(lon.min,",",lon.max,",",lat.min,",",lat.max))
  cdo.cmd <- paste0("cdo muldpm -fldmean -sellonlatbox,", boundaries," ",
                    dir.data.ModERA,Set,"/abs/ensstat/by_var/mon/ModE-RA_lowres_20mem_Set_1420-3_1850-1_ensmean_totprec_abs_1421-2008_mon.nc ",
                    "aux.nc")
  system(cdo.cmd)

  f <- nc_open("aux.nc")
  P.memb[[j]] <- ncvar_get(f, varid="totprec") * 86400 # per seconds to per day
  
  file.remove("aux.nc")
}

Dates <- seq(as.Date("1421-01-01"), as.Date("2008-12-31"), by="month")
P.memb <- as.data.frame(P.memb) %>% cbind(Dates,.)
save(P.memb, file="ModE-RA_lowres_20mem_Set_1420-3_1850-1_totprec-Boxes_1421-2008_mon.RData")
