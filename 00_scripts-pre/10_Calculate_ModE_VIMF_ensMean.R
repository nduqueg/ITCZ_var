rm(list = ls())
cat("\014")

"%>%"=magrittr::`%>%`

source("00_settings.R")
# dir.base.ModEsim <- "/mnt/climstor/ERC_PALAEO/ModE-Sim/outdata/"

################################-
## generate names in ModE-Sim data ----
################################-

# names of directories
sets <- data.frame( Set= c("set_1420-1_to_3","set_1850-1_to_2"),
                    Epoch= c(1,2),
                    stringsAsFactors = F
)

# time periods for each member
Epoch <- list(Epoch1 = seq(1420,1849), Epoch2= seq(1850,2009))

plev <- c("100000","85000","70000","50000","30000","20000","10000")

################################-
## fetch variables and calculate VIMF ----
################################-

setwd(paste0(dir.base,"01_Data/10_VIMF/"))

for ( i in sets$Set){
  
  # creating the directories and temporal directories

  iEpoch <- subset(sets, Set==i, select=Epoch)[1,]
  Per <- Epoch[[iEpoch]] %>% paste0("-",tail(.,1)) %>% .[1]
  if (iEpoch ==1){
    tslice <- " -seldate,1420-01-01,1849-12-31 "
  }else{
    tslice <- " -seldate,1850-01-01,2009-12-31 "
  }

  # ## join individual levels in q, u and v
  # z.t <- read.table("myzaxisQ.txt")
  # cdo.merge.Q <- cdo.merge.U <- cdo.merge.V <- "cdo merge "
  # if(!dir.exists("dummy_Q")) dir.create("dummy_Q")
  # if(!dir.exists("dummy_U")) dir.create("dummy_U")
  # if(!dir.exists("dummy_V")) dir.create("dummy_V")
  # 
  # for (j in plev){
  #   # modify the axis file
  #   z.t[6,3] <- j
  #   write.table(z.t, "myzaxisQ.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
  # 
  #   cdo.cmd <- paste0("cdo setzaxis,myzaxisQ.txt ",dir.data.ModEsim, i,"/abs/ensstat/by_var/mon/ModE-Sim_",i,"_ensmean_q_",j,"_abs_",Per,"_mon.nc ",
  #                                    "./dummy_Q/ModE-Sim_",i,"_ensmean_q_",j,"_abs_",Per,"_mon.nc")
  #   system(cdo.cmd)
  # 
  #   cdo.merge.Q <- paste0(cdo.merge.Q, "./dummy_Q/ModE-Sim_",i,"_ensmean_q_",j,"_abs_",Per,"_mon.nc ")
  # 
  #   cdo.cmd <- paste0("cdo setzaxis,myzaxisQ.txt ",dir.data.ModEsim, i,"/abs/ensstat/by_var/mon/ModE-Sim_",i,"_ensmean_u_",j,"_abs_",Per,"_mon.nc ",
  #                     "./dummy_U/ModE-Sim_",i,"_ensmean_u_",j,"_abs_",Per,"_mon.nc")
  #   system(cdo.cmd)
  # 
  #   cdo.merge.U <- paste0(cdo.merge.U, "./dummy_U/ModE-Sim_",i,"_ensmean_u_",j,"_abs_",Per,"_mon.nc ")
  # 
  #   cdo.cmd <- paste0("cdo setzaxis,myzaxisQ.txt ",dir.data.ModEsim, i,"/abs/ensstat/by_var/mon/ModE-Sim_",i,"_ensmean_v_",j,"_abs_",Per,"_mon.nc ",
  #                     "./dummy_V/ModE-Sim_",i,"_ensmean_v_",j,"_abs_",Per,"_mon.nc")
  #   system(cdo.cmd)
  # 
  #   cdo.merge.V <- paste0(cdo.merge.V, "./dummy_V/ModE-Sim_",i,"_ensmean_v_",j,"_abs_",Per,"_mon.nc ")
  # }
  # 
  # # merge vertical files in just one
  # cdo.merge.Q <- paste0(cdo.merge.Q, "ModE-Sim_",i,"_ensmean_q_",Per,"_mon.nc")
  # cdo.merge.U <- paste0(cdo.merge.U, "ModE-Sim_",i,"_ensmean_u_",Per,"_mon.nc")
  # cdo.merge.V <- paste0(cdo.merge.V, "ModE-Sim_",i,"_ensmean_v_",Per,"_mon.nc")
  # system(cdo.merge.Q); system(cdo.merge.U); system(cdo.merge.V)
  # 
  # for(k in c("Q","U","V")){
  #   f.rm <- list.files(path = paste0("./dummy_", k,"/")); file.remove(paste0("./dummy_",k) %>% paste0(.,"/",f.rm))
  # }
  # file.remove(c("./dummy_Q/","./dummy_U/","./dummy_V/"))

  ## ------------------------------------------------------------------------------- calculate the VIMF

  cdo.cmd <- paste0("cdo mul ",
                    "ModE-Sim_",i,"_ensmean_u_",Per,"_mon.nc ", # u
                    "ModE-Sim_",i,"_ensmean_q_",Per,"_mon.nc ", # Q
                    "dummy_VIMFu.nc"); system(cdo.cmd)
  cdo.cmd <- paste0("cdo mul ",
                    "ModE-Sim_",i,"_ensmean_v_",Per,"_mon.nc ", # v
                    "ModE-Sim_",i,"_ensmean_q_",Per,"_mon.nc ", # Q
                    "dummy_VIMFv.nc"); system(cdo.cmd)

  cdo.cmd <-paste0("cdo vertsum -mul -divc,-9.80665",tslice,"dp.nc ", "dummy_VIMFu.nc ", # vertsum -mul -divc,9.80665 dp.nc v*q
                   "dummy_VIMFu2.nc")
  system(cdo.cmd)
  cdo.cmd <-paste0("cdo vertsum -mul -divc,-9.80665",tslice,"dp.nc ", "dummy_VIMFv.nc ", # vertsum -mul -divc,9.80665 dp.nc v*q
                   "dummy_VIMFv2.nc")
  system(cdo.cmd)

  cdo.cmd <- paste0("cdo chunit,Pa,kgm**-1s**-1 -chname,dp,u dummy_VIMFu2.nc ","ModE-Sim_",i,"_ensmean_",Per,"_VIMFeast_mon.nc")
  system(cdo.cmd)
  cdo.cmd <- paste0("cdo chunit,Pa,kgm**-1s**-1 -chname,dp,v dummy_VIMFv2.nc ","ModE-Sim_",i,"_ensmean_",Per,"_VIMFnorth_mon.nc")
  system(cdo.cmd)
  
  f.delete <- list.files(pattern = "dummy_")
  file.remove(f.delete)
  
}
setwd(dir.base)
