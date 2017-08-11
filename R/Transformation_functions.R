
get.name<-function(x){
  nm <-deparse(substitute(x))
  return(nm)
}



Monomolecular <- function(r, parm){
  parm[3]*(1-exp(-1*r/parm[2]))+1 # Monomolecular
}

Inv.Monomolecular <- function(r, parm){
  if(class(r)=="RasterLayer") {
    R <- parm[3]*(exp(-1*r/parm[2])) 
    (R <- (R-cellStats(R,stat = "min"))+1)  
  } else {
    R <- parm[3]*(exp(-1*r/parm[2])) 
    (R <- (R-min(R))+1)
  }
}

Inv.Rev.Monomolecular <- function(r, parm){
  if(class(r)=="RasterLayer") {
    rev.rast <- SCALE((-1*r),0,10)
    Inv.Monomolecular(rev.rast,parm)
  } else {
    rev.rast <- SCALE.vector((-1*r),0,10)
    Inv.Monomolecular(rev.rast,parm)
  }
}

Rev.Monomolecular <- function(r, parm){
  if(class(r)=="RasterLayer") {
    rev.rast <- SCALE((-1*r),0,10)
    Monomolecular(rev.rast,parm)
  } else {
    rev.rast <- SCALE.vector((-1*r),0,10)
    Monomolecular(rev.rast,parm)
  }
}


Ricker <- function(r,parm){
  parm[3]*r*exp(-1*r/parm[2])+1 # Ricker
}

Inv.Ricker <- function(r,parm){
  if(class(r)=="RasterLayer") {
    R <- (-1*parm[3])*r*exp(-1*r/parm[2])-1 # Ricker
    R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
  } else {    
    R <- (-1*parm[3])*r*exp(-1*r/parm[2])-1 # Ricker
    R <- SCALE.vector(R,MIN=abs(max(R)),MAX=abs(min(R))) # Rescale
  }
}

Inv.Rev.Ricker <- function(r,parm){
  if(class(r)=="RasterLayer") {
    rev.rast <- SCALE((-1*r),0,10)
    Inv.Ricker(rev.rast,parm)
  } else {
    rev.rast <- SCALE.vector((-1*r),0,10)
    Inv.Ricker(rev.rast,parm)
  }
}

Rev.Ricker <- function(r,parm){
  if(class(r)=="RasterLayer") {
    rev.rast <- SCALE((-1*r),0,10)
    Ricker(rev.rast,parm)
  } else {
    rev.rast <- SCALE.vector((-1*r),0,10)
    Ricker(rev.rast,parm)
  }
}
