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
Epoch <- list(Epoch1 = seq(1420,1849), Epoch2= seq(1850,2009))

################################-
## fetch members Humidity, join levels in one fileand zonal mean----
################################-

# the indiviual members and the ensamble of the variable are already calculated and separated in
# set_1420-1_to_3/abs/ensstat/by_var/mon
# set_1420-1/abs/m001/by_var/mon

# numCores <- 20
# cl <- makeSOCKcluster(numCores)
# registerDoSNOW(cl)
# 
# pb <- txtProgressBar(min=1, max= nrow(sets), style=3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress=progress)

setwd("./00_scripts-pre/")
plev <- c("100000","85000","70000","50000","30000","20000","10000")
# results <- foreach ( i = 1:nrow(sets), .combine=rbind, .options.snow=opts) %dopar% {
for ( i in 1:nrow(sets)){
  with(sets[i,], paste0(Set,"_",Memb)) %>% print()
  
  iEpoch <- sets[i,"Epoch"]
  Per <- Epoch[[iEpoch]] %>% paste0(.,"-",tail(.,1)) %>% .[1]
  
  z.t <- read.table("myzaxisQ.txt")
  cdo.merge <- "cdo zonmean -merge "
  
  for (j in plev){
    # modify the axis file
    
    z.t[6,3] <- j
    write.table(z.t, "myzaxisQ.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    cdo.cmd <- with(sets[i,], paste0("cdo setzaxis,myzaxisQ.txt ",dir.data.ModEsim,Set,"/abs/",Memb,"/by_var/mon/ModE-Sim_",Set,"_",Memb,"_q_",j,"_abs_",Per,"_mon.nc ",
                                     "./dummy_Q/ModE-Sim_",Set,"_",Memb,"_q_",j,"_abs_",Per,"_mon.nc"))
    system(cdo.cmd)
    
    cdo.merge <- with(sets[i,], paste0(cdo.merge, "./dummy_Q/ModE-Sim_",Set,"_",Memb,"_q_",j,"_abs_",Per,"_mon.nc "))
  }
  
  if( !dir.exists( "../01_Data/04_q/" %>% paste0(.,sets$Set[i]) )) dir.create("../01_Data/04_q/" %>% paste0(.,sets$Set[i]))
  
  cdo.merge <- with(sets[i,], paste0(cdo.merge, "../01_Data/04_q/", sets$Set[i],"/ModE-Sim_",Set,"_",Memb,"_q-ZonMean_",Per,"_mon.nc") )
  system(cdo.merge)
  
  f.rm <- list.files(path = "./dummy_Q/", pattern = with(sets[i,], paste0("ModE-Sim_",Set,"_",Memb)) )
  file.remove("./dummy_Q/" %>% paste0(.,f.rm))
}

# stopCluster(cl)

################################-
# zonal mean for the general ensemble ----
################################-

# first epoch --------------------------------------------------------- s
{
  cdo.merge <- "cdo zonmean -merge "
  Per <- "1420-1849"
  for (j in plev){
    # modify the axis file
    
    z.t[6,3] <- j
    write.table(z.t, "myzaxisQ.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    cdo.cmd <- with(sets[i,], paste0("cdo setzaxis,myzaxisQ.txt ",dir.data.ModEsim,"set_1420-1_to_3/abs/ensstat/by_var/mon/ModE-Sim_set_1420-1_to_3_ensmean_q_",j,"_abs_",Per,"_mon.nc ",
                                     "./dummy_Q/ModE-Sim_set_1420-1_to_3_q_",j,"_abs_",Per,"_mon.nc"))
    system(cdo.cmd)
    
    cdo.merge <- with(sets[i,], paste0(cdo.merge, "./dummy_Q/ModE-Sim_set_1420-1_to_3_q_",j,"_abs_",Per,"_mon.nc "))
  }
  cdo.merge <- with(sets[i,], paste0(cdo.merge, "../01_Data/04_q/ModE-Sim_set_1420-1_to_3_q-ZonMean_",Per,"_mon.nc") )
  system(cdo.merge)
  
  f.rm <- list.files(path = "./dummy_Q/" )
  file.remove("./dummy_Q/" %>% paste0(.,f.rm))
}

# second epoch --------------------------------------------------------- s
{
  cdo.merge <- "cdo zonmean -merge "
  Per <- "1850-2009"
  for (j in plev){
    # modify the axis file
    
    z.t[6,3] <- j
    write.table(z.t, "myzaxisQ.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    cdo.cmd <- with(sets[i,], paste0("cdo setzaxis,myzaxisQ.txt ",dir.data.ModEsim,"set_1850-1_to_2/abs/ensstat/by_var/mon/ModE-Sim_set_1850-1_to_2_ensmean_q_",j,"_abs_",Per,"_mon.nc ",
                                     "./dummy_Q/ModE-Sim_set_1850-1_to_2_q_",j,"_abs_",Per,"_mon.nc"))
    system(cdo.cmd)
    
    cdo.merge <- with(sets[i,], paste0(cdo.merge, "./dummy_Q/ModE-Sim_set_1850-1_to_2_q_",j,"_abs_",Per,"_mon.nc "))
  }
  cdo.merge <- with(sets[i,], paste0(cdo.merge, "../01_Data/04_q/ModE-Sim_set_1850-1_to_2_q-ZonMean_",Per,"_mon.nc") )
  system(cdo.merge)
  
  f.rm <- list.files(path = "./dummy_Q/" )
  file.remove("./dummy_Q/" %>% paste0(.,f.rm))
}


################################-
# ensemble mean zonal mean for each set ensemble ----
################################-
sets.names <- unique(sets$Set)
epoch.lim <- c(rep("1420-1849",3), rep("1850-2009",2))

for ( i in sets.names){
  Per <- which(i == sets.names) %>% epoch.lim[.]
  
  cdo.merge <- "cdo zonmean -merge "
  
  for (j in plev){
    # modify the axis file
    
    z.t[6,3] <- j
    write.table(z.t, "myzaxisQ.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    cdo.cmd <- paste0("cdo setzaxis,myzaxisQ.txt ",dir.data.ModEsim,i,"/abs/ensstat/by_var/mon/ModE-Sim_",i,"_ensmean_q_",j,"_abs_",Per,"_mon.nc ",
                                     "./dummy_Q/ModE-Sim_",i,"_ensmean_q_",j,"_abs_",Per,"_mon.nc")
    system(cdo.cmd)
    
    cdo.merge <- with(sets[i,], paste0(cdo.merge, "./dummy_Q/ModE-Sim_",i,"_ensmean_q_",j,"_abs_",Per,"_mon.nc "))
  }
  
  if( !dir.exists( "../01_Data/04_q/" %>% paste0(.,sets$Set[i]) )) dir.create("../01_Data/04_q/" %>% paste0(.,sets$Set[i]))
  
  cdo.merge <-  paste0(cdo.merge, "../01_Data/04_q/", i,"/ModE-Sim_",i,"_ensmean_q-ZonMean_",Per,"_mon.nc")
  system(cdo.merge)
  
  f.rm <- list.files(path = "./dummy_Q/", pattern = with(sets[i,], paste0("ModE-Sim_",i)) )
  file.remove("./dummy_Q/" %>% paste0(.,f.rm))
  
}



################################-
# aggregations and selection  ----
################################-

setwd("../01_Data/04_q/")

cdo.cmd <- "cdo seasmean ModE-Sim_set_1420-1_to_3_q-ZonMean_1420-1849_mon.nc ModE-Sim_set_1420-1_to_3_q-ZonMean_1420-1849_seasonal.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,DJF ModE-Sim_set_1420-1_to_3_q-ZonMean_1420-1849_seasonal.nc ModE-Sim_set_1420-1_to_3_q-ZonMean_1420-1849_sDJF.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,JJA ModE-Sim_set_1420-1_to_3_q-ZonMean_1420-1849_seasonal.nc ModE-Sim_set_1420-1_to_3_q-ZonMean_1420-1849_sJJA.nc"
system(cdo.cmd)

cdo.cmd <- "cdo seasmean ModE-Sim_set_1850-1_to_2_q-ZonMean_1850-2009_mon.nc ModE-Sim_set_1850-1_to_2_q-ZonMean_1850-2009_seasonal.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,DJF ModE-Sim_set_1850-1_to_2_q-ZonMean_1850-2009_seasonal.nc ModE-Sim_set_1850-1_to_2_q-ZonMean_1850-2009_sDJF.nc"
system(cdo.cmd)
cdo.cmd <- "cdo selseas,JJA ModE-Sim_set_1850-1_to_2_q-ZonMean_1850-2009_seasonal.nc ModE-Sim_set_1850-1_to_2_q-ZonMean_1850-2009_sJJA.nc"
system(cdo.cmd)



# seasonal accumulation and selection

for ( i in sets.names){
  Per <- which(i == sets.names) %>% epoch.lim[.]
  cdo.cmd <- paste0("cdo seasmean ./",i,"/ModE-Sim_",i,"_ensmean_q-ZonMean_",Per,"_mon.nc ./",i,"/ModE-Sim_",i,"_ensmean_q-ZonMean_",Per,"_seasonal.nc")
  system(cdo.cmd)
  
  cdo.cmd <- paste0("cdo selseas,DJF ./",i,"/ModE-Sim_",i,"_ensmean_q-ZonMean_",Per,"_seasonal.nc ./",i,"/ModE-Sim_",i,"_ensmean_q-ZonMean_",Per,"_sDJF.nc")
  system(cdo.cmd)
  
  cdo.cmd <- paste0("cdo selseas,JJA ./",i,"/ModE-Sim_",i,"_ensmean_q-ZonMean_",Per,"_seasonal.nc ./",i,"/ModE-Sim_",i,"_ensmean_q-ZonMean_",Per,"_sJJA.nc")
  system(cdo.cmd)

}


setwd(dir.base)
