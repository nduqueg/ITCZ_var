rm(list = ls())
cat("\014")

library(parallel)
library(doSNOW)
library(tcltk)
"%>%"=magrittr::`%>%`

dir.base <- "/mnt/climstor/ERC_PALAEO/ModE-Sim/outdata/"

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

## fetch Meridional streamflow ----

l.rerun <- c()
setwd("../01_Data/01_Streamfunction/")
for ( i in 1:nrow(sets)){

  # creating the directories and temporal directories
  if( !dir.exists(sets$Set[i])) dir.create(sets$Set[i])
  setwd(sets$Set[i])
  if (!dir.exists(sets$Memb[i])) dir.create(sets$Memb[i])
  setwd(sets$Memb[i])
  iEpoch <- sets[i,"Epoch"]

  # identifying the files that were "rerun"
  aux <- with(sets[i,], paste0(dir.base,Set,"/abs/",Memb,"/by_year/")) %>% list.files(.,pattern=".meridional.nc")
  name.files.rerun <- aux[aux %>% substr(.,1,5)=="rerun"]; l.rerun <- c(l.rerun,name.files.rerun)
  
  # setting up the cluster for parallel processing
  numCores <- 10
  cl <- makeSOCKcluster(numCores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max= length(Epoch[[iEpoch]]), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  result <- foreach (j = Epoch[[iEpoch]], .combine=rbind, .packages = c(), .options.snow=opts) %dopar% {
    # print(with(sets[i,],paste(Set,"/",Memb,"Year-",j)))

    cdo.cmd <- with(sets[i,], paste0("cdo -selname,mastrfu ",dir.base,Set,"/abs/",Memb,"/by_year/ModE-Sim_",Set,"_",Memb,"_",j,"_mon.meridional.nc ",
                                     "ModE-Sim_",Set,"_",Memb,"_",j,"_mastrf_mon.nc"))
    system(cdo.cmd)
    return(cdo.cmd)
  }
  stopCluster(cl)

  # making the re runs
  cdo.cmd <- with(sets[i,], paste0("cdo -selname,mastrfu ",dir.base,Set,"/abs/",Memb,"/by_year/",name.files.rerun,
                                   " ModE-Sim_",Set,"_",Memb,"_",substr(name.files.rerun,32,35),"_mastrf_mon.nc"))
  for(k in cdo.cmd) system(k)
  
  # pasting the individual yearly files
  setwd("..")
  Per <- Epoch[[iEpoch]] %>% paste0("-",tail(.,1)) %>% .[1]
  cdo.cmd <- paste0("cdo mergetime ./",sets$Memb[i],"/*.nc ","ModE-Sim_",sets$Set[i],"_",sets$Memb[i],"_",Per,"_mastrf_mon.nc")
  system(cdo.cmd)

  f.delete <- list.files(path = sets$Memb[i])
  file.remove(paste0(sets$Memb[i],"/",f.delete)); file.remove(sets$Memb[i])

  setwd("..")
}

cdo.cmd <- "cdo ensmean ./set_1420-1/*.nc ./set_1420-2/*.nc ./set_1420-3/*.nc ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_mon.nc"
system(cdo.cmd)

cdo.cmd <- "cdo ensmean ./set_1850-1/*.nc ./set_1850-2/*.nc ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_mon.nc"
system(cdo.cmd)

cdo.cmd <- "cdo seasmean ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_mon.nc ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_seasonal.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,DJF ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_seasonal.nc ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_sDJF.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,JJA ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_seasonal.nc ModE-Sim_set_1420-1_to_3_ensmean_mastrf_1420-1849_sJJA.nc"
system(cdo.cmd)

cdo.cmd <- "cdo seasmean ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_mon.nc ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_seasonal.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,DJF ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_seasonal.nc ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_sDJF.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,JJA ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_seasonal.nc ModE-Sim_set_1850-1_to_2_ensmean_mastrf_1850-2009_sJJA.nc"
system(cdo.cmd)

