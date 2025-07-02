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

################################-
## fetch variables and calculate VIMF ----
################################-

setwd(paste0(dir.base,"01_Data/10_VIMF/"))

for ( i in sets$Set){
  
  # creating the directories and temporal directories
  
  iEpoch <- subset(sets, Set==i, select=Epoch)[1,]
  Per <- Epoch[[iEpoch]] %>% paste0("-",tail(.,1)) %>% .[1]
  
  ## ------------------------------------------------------------------------------- calculate the MDiv
  
  cdo.cmd <- paste0("cdo sp2gp -uv2dv -merge ",
                    "ModE-Sim_",i,"_ensmean_",Per,"_VIMFeast_mon.nc ", # u
                    "ModE-Sim_",i,"_ensmean_",Per,"_VIMFnorth_mon.nc ", # v
                    "dummy.nc")
  system(cdo.cmd)
  
  system(paste0("cdo chunit,1/s,kg/m**2/s -chname,sd,vimd -remapbil,r360x180 dummy.nc ModE-Sim_",i,"_ensmean_",Per,"_MDiv_mon.nc"))
  file.remove("dummy.nc")
}

setwd(dir.base)
