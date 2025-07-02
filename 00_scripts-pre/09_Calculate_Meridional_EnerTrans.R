rm(list = ls())
cat("\014")

library(parallel)
library(doSNOW)
library(tcltk)
"%>%"=magrittr::`%>%`


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
Epoch <- list(Epoch1 = seq(1420,1849), Epoch2= seq(1850,2009))

################################-
## fetch variables and calculate Total specific energy ----
################################-

setwd(paste0(dir.base,"01_Data/09_VIEnerTrans/"))
dir.SpecEner <- paste0(dir.base,"01_Data/08_TotSpecEnergy/")

sets.names <- unique(sets$Set)
for ( i in sets.names){
  
  # creating the directories and temporal directories
  if( !dir.exists(i)) dir.create(i)
  setwd(i)
  iEpoch <- subset(sets, Set==i, select=Epoch)[1,]
  Per <- Epoch[[iEpoch]] %>% paste0("-",tail(.,1)) %>% .[1]
  if (iEpoch ==1){
    Per.mod <- "1420-1850"
    tslice <- " -seldate,1420-01-01,1849-12-31 "
  }else{
    Per.mod <- Per
    tslice <- " -seldate,1850-01-01,2009-12-31 "
  } 
  
  # setting up the cluster for parallel processing
  numCores <- 10
  cl <- makeSOCKcluster(numCores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max= nrow(subset(sets,Set==i,select=Memb)), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  result <- foreach (Memb = subset(sets,Set==i,select=Memb)[,1], .combine=rbind, .packages = c(), .options.snow=opts,
                     .export = c("iEpoch","Per","Per.mod","tslice","i","dir.data.ModEsim","dir.SpecEner")) %dopar% {
                       # for ( i in 2:nrow(sets)){ # run each one of the sets and members
                       Set <- i
                       
                       cdo.cmd <- paste0("cdo mul ",
                                         dir.data.ModEsim,Set,"/abs/",Memb,"/by_var/mon/ModE-Sim_",Set,"_",Memb,"_v_abs_",Per,"_mon.nc", # v
                                         tslice,
                                         dir.SpecEner,Set,"/ModE-Sim_",Set,"_",Memb,"_",Per.mod,"_TotSpecEnergy_mon.nc ", # E
                                         Memb,"_dummy.nc")
                       system(cdo.cmd)
                       
                       cdo.cmd <-paste0("cdo vertsum -mul -divc,9.80665",tslice,"../dp.nc ", Memb,"_dummy.nc ", # vertsum -mul -divc,9.80665 dp.nc v*E
                                        "ModE-Sim_",Set,"_",Memb,"_",Per,"_VITEFnorth_mon.nc")
                       system(cdo.cmd)
                       
                       
                     }
  stopCluster(cl)
  f.delete <- list.files(pattern = "_dummy.nc")
  file.remove(f.delete)
  
  setwd("..")
}
################################-
# ensemble mean for the 3 sets in epoch 1 and 2 sets in epoch 2 ----
################################-

cdo.cmd <- "cdo ensmean ./set_1420-1/*.nc ./set_1420-2/*.nc ./set_1420-3/*.nc ModE-Sim_set_1420-1_to_3_ensmean_VITEFnorth_1420-1849_mon.nc"
system(cdo.cmd)

cdo.cmd <- "cdo ensmean ./set_1850-1/*.nc ./set_1850-2/*.nc ModE-Sim_set_1850-1_to_2_ensmean_VITEFnorth_1850-2009_mon.nc"
system(cdo.cmd)

################################-
# ensemble mean for each sets individually ----
################################-

sets.names <- unique(sets$Set)
epoch.lim <- c(rep("1420-1849",3), rep("1850-2009",2))

for ( i in sets.names){
  epoch.years <- which(i == sets.names) %>% epoch.lim[.]
  
  cdo.cmd <- paste0("cdo ensmean ./",i,"/*.nc ./",i,"/ModE-Sim_",i,"_ensmean_VITEFnorth_",epoch.years,"_mon.nc")
  # print(cdo.cmd)
  system(cdo.cmd)
}



setwd(dir.base)