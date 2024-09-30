
###############################-
# Insert here the local directory that contains the ModE-Sim data
###############################-

dir.data.ModEsim <- "/mnt/climstor/ERC_PALAEO/ModE-Sim/outdata/"


###############################-
# Insert here the local working directory
###############################-

dir.base <- "/scratch2/nduque/50_ITCZ/ITCZ_var/"


###############################-
# Insert here the local directory that contains the ModE-RA data
###############################-

dir.data.ModERA <- "/mnt/climstor/ERC_PALAEO/ModE-RA/outdata/"


###############################-
# Function for identifying the location of the maximum or minimum value of a ITCZ feature
###############################-
"%>%"=magrittr::`%>%`
smt.min.max <- function(x,type, lat){
  #
  # 'x' is the ITCZ feature and must have the time variable in lines and the latitude in the columns
  # It will be smoothed around the minimum or maximum with an spline curve and the position will be also identified with three values around the models' maximum or minimum
  #
  
  aprox <- x %>% apply(.,1, paste0("which.",type) %>% get()) %>% as.numeric() # position of the model's maximum or minimum
  
  ind.aprox <- cbind(aprox-3,aprox-2,aprox-1,aprox,aprox+1,aprox+2,aprox+3) %>% t() %>% as.data.frame() # indices of values close to the Max StrFunct for the Spline model
  lat.models <- apply(ind.aprox, 2, function(ind,lat) lat[ind], lat) %>%  as.data.frame() # latitudes for the Spline model
  
  data.spModel <- mapply(FUN= function(x,ind) x[ind], # extracts the MAx StrFunct values near the Max for the Spline model
                         x[,lat >= -45 & lat <= 45] %>% as.matrix() %>% split(.,row(.)), 
                         ind.aprox,
                         SIMPLIFY = F) %>% 
    
    mapply(function(y,lat) cbind.data.frame(lat,y), # paste for each time step the latitude and the values near the Max StrFunc
           .,
           lat.models, 
           SIMPLIFY = F) %>% 
    lapply(., `colnames<-`,c("x.mod","y.mod"))
  
  sp.models <- lapply(data.spModel, function(x) tryCatch( # there are some months that don't have values
    {lm(y.mod ~ splines::ns(x.mod,3), data= x)},
    error= function(cond){NA}
  )
  )
  pos.feat <- numeric()
  for(k in 1:length(sp.models)){ # for each year identify the position of the Strong Ascent
    
    if (is.list(sp.models[[k]])){
      
      pos.feat[k] <- predict(sp.models[[k]], newdata = data.frame(x.mod= seq(lat.models[7,k],lat.models[1,k],by=0.1))) %>% as.numeric(.) %>% matrix(.,nrow=1) %>% 
        apply(.,1,get(paste0("which.",type))) %>% 
        seq(lat.models[7,k],lat.models[1,k],by=0.1)[.]
      
    } else{
      pos.feat[k] <- NA
    }
  }
  return(pos.feat)
}