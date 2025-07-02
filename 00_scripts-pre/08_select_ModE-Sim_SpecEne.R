rm(list = ls())
cat("\014")

library(parallel)
library(doSNOW)
library(tcltk)
"%>%"=magrittr::`%>%`

source("00_settings.R")
# dir.base.ModEsim <- "/mnt/climstor/ERC_PALAEO/ModE-Sim/outdata/"

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
  )

# time periods for each member
Epoch <- list(Epoch1 = seq(1420,1850), Epoch2= seq(1850,2009))

################################-
## fetch variables and calculate Total specific energy ----
################################-

l.rerun <- c()
setwd(paste0(dir.base,"/01_Data/08_TotSpecEnergy/"))


for ( i in 1:nrow(sets)){ # run each one of the sets and members
  
  # creating the directories and temporal directories
  if( !dir.exists(sets$Set[i])) dir.create(sets$Set[i])
  setwd(sets$Set[i])
  if (!dir.exists(sets$Memb[i])) dir.create(sets$Memb[i])
  setwd(sets$Memb[i])
  iEpoch <- sets[i,"Epoch"]
  
  # identifying the files that were "rerun"
  aux <- with(sets[i,], paste0(dir.data.ModEsim,Set,"/abs/",Memb,"/by_year/")) %>% list.files(.,pattern=".stdlevs.nc")
  name.files.rerun <- aux[aux %>% substr(.,1,5)=="rerun"]; l.rerun <- c(l.rerun,name.files.rerun)
  
  # setting up the cluster for parallel processing
  numCores <- 10
  cl <- makeSOCKcluster(numCores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max= length(Epoch[[iEpoch]]), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  # processing each year
  
  # Cp= 1005 J/(kg*K);  Lv= 2256400 J/kg
  # operation <- "expr,'E= (1005)*t + (2256.4)*q + 9.80665*geopoth + (u^2+v^2)/2; v=v; u=u' -selname,u,v,t,q,geopoth"
  operation <- "expr,'E= (1005)*t + (2256.4)*q + 9.80665*geopoth + (u^2+v^2)/2;' -selname,u,v,t,q,geopoth"
  result <- foreach (j = Epoch[[iEpoch]], .combine=rbind, .packages = c(), .options.snow=opts) %dopar% {
    
    
    cdo.cmd <- with(sets[i,], paste0("cdo ",operation," ",
                                     dir.data.ModEsim,Set,"/abs/",Memb,"/by_year/ModE-Sim_",Set,"_",Memb,"_",j,"_mon.stdlevs.nc ",
                                     "ModE-Sim_",Set,"_",Memb,"_",j,"_TotSpecEnergy_mon.nc"))
    system(cdo.cmd)
    
    # system(with(sets[i,], paste0("nccopy -k 4 ModE-Sim_",Set,"_",Memb,"_",j,"_TotSpecEnergy_mon.nc ModE-Sim_",Set,"_",Memb,"_",j,"_TotSpecEnergy_mon.nc")))
    return(cdo.cmd)
  }
  stopCluster(cl)
  
  # making the re runs
  cdo.cmd <- with(sets[i,], paste0("cdo ",operation," ",
                                   dir.data.ModEsim,Set,"/abs/",Memb,"/by_year/",name.files.rerun,
                                   " ModE-Sim_",Set,"_",Memb,"_",substr(name.files.rerun,32,35),"_TotSpecEnergy_mon.nc"))
  for(k in cdo.cmd) system(k)
  
  # pasting the individual yearly files
  setwd("..")
  Per <- Epoch[[iEpoch]] %>% paste0("-",tail(.,1)) %>% .[1]
  cdo.cmd <- paste0("cdo -z zip mergetime ./",sets$Memb[i],"/*.nc ","ModE-Sim_",sets$Set[i],"_",sets$Memb[i],"_",Per,"_TotSpecEnergy_mon.nc")
  system(cdo.cmd)
  
  f.delete <- list.files(path = sets$Memb[i])
  file.remove(paste0(sets$Memb[i],"/",f.delete)); file.remove(sets$Memb[i])
  
  setwd("..")
}

################################-
# ensemble mean for the 3 sets in epoch 1 and 2 sets in epoch 2 ----
################################-

cdo.cmd <- "cdo ensmean ./set_1420-1/*.nc ./set_1420-2/*.nc ./set_1420-3/*.nc ModE-Sim_set_1420-1_to_3_ensmean_TotSpecEnergy_1420-1850_mon.nc"
system(cdo.cmd)

cdo.cmd <- "cdo ensmean ./set_1850-1/*.nc ./set_1850-2/*.nc ModE-Sim_set_1850-1_to_2_ensmean_TotSpecEnergy_1850-2009_mon.nc"
system(cdo.cmd)

################################-
# ensemble mean for each sets individually ----
################################-

sets.names <- unique(sets$Set)
epoch.lim <- c(rep("1420-1850",3), rep("1850-2009",2))

for ( i in sets.names){
  epoch.years <- which(i == sets.names) %>% epoch.lim[.]
  
  cdo.cmd <- paste0("cdo ensmean ./",i,"/*.nc ./",i,"/ModE-Sim_",i,"_ensmean_TotSpecEnergy_",epoch.years,"_mon.nc")
  system(cdo.cmd)
}



setwd(dir.base)