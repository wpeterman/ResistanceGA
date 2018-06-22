#' Apply transformation to continuous resistance surface
#'
#' Apply one the eight resistance transformations to a continuous resistance surface
#'
#' @param transformation Transformation equation to apply. Can be provided as the name of the transformation or its numeric equivalent (see details)
#' @param shape Value of the shape parameter
#' @param max Value of the maximum value parameter
#' @param scale The standard deviation, in number of raster cells, to use when applying Gaussian kernel smoothing. This is the `sigma` parameter in the `spatstat::blur` function. (Default = NULL)
#' @param r Resistance surface to be transformed. Can be supplied as full path to .asc file or as a raster object
#' @param out Directory to write transformed .asc file. Default is NULL, and will not export .asc file
#' @usage Resistance.tran(transformation, shape, max, scale, r, out)
#' @return R raster object
#' @details Valid arguements for \code{transformation} are:\cr
#' \tabular{ll}{
#'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#'    \tab 2 = "Inverse-Reverse Ricker"\cr
#'    \tab 3 = "Monomolecular"\cr
#'    \tab 4 = "Ricker"\cr
#'    \tab 5 = "Reverse Monomolecular"\cr
#'    \tab 6 = "Reverse Ricker"\cr
#'    \tab 7 = "Inverse Monomolecular"\cr
#'    \tab 8 = "Inverse Ricker"\cr
#'    \tab 9 = "Distance"\cr
#'    }
#'
#' The Distance transformation sets all values equal to one. Because of the flexibility of the Ricker function to take a monomolecular shape (try \code{Plot.trans(PARM=c(10,100), Resistance=c(1,10), transformation="Ricker")} to see this), whenever a shape parameter >6 is selected in combination with a Ricker family transformation, the transformation reverts to a Distance transformation. In general, it seems that using a combination of intermediate Ricker and Monomolecular transformations provides the best, most flexible coverasge of parameter space.
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>

Resistance.tran <- function(transformation,
                            shape,
                            max,
                            scale = NULL,
                            r,
                            out = NULL) {
  if (class(r)[1] != 'RasterLayer') {
    R <- raster(r)
    NAME <- basename(r)
    NAME <- sub("^([^.]*).*", "\\1", NAME)
    names(R) <- NAME
  } else {
    R <- r
    NAME <- r@data@names
  }
  if (is.numeric(transformation)) {
    parm <- c(transformation, shape, max)
  } else {
    parm <- c(get.EQ(transformation), shape, max)
  }
  
  EXPORT.dir <- out
  
  
  
  if (!is.null(scale)) {
    R <- k.smooth(raster = R,
                  sigma = scale,
                  SCALE = TRUE)
    
    # zmat <- as.matrix(R)
    # # Note that sigma is in pixels, NOT map units!
    # x <- spatstat::as.im(zmat)
    # 
    # f <- spatstat::blur(x = x,
    #                     sigma = scale,
    #                     normalise = TRUE,
    #                     bleed = FALSE)
    # R <- R
    # values(R) <- f$v
  }
  
  R <- SCALE(data = R, MIN = 0, MAX = 10)
  
  # Set equation for continuous surface
  equation <- floor(parm[1]) # Parameter can range from 1-9.99
  
  #       # Read in resistance surface to be optimized
  #       SHAPE <-  (parm[2])
  #       Max.SCALE <- (parm[3])
  
  
  
  
  # Apply specified transformation
  if (equation == 1) {
    r <- Inv.Rev.Monomolecular(R, parm)
    EQ <- "Inverse-Reverse Monomolecular"
    
  } else if (equation == 5) {
    r <- Rev.Monomolecular(R, parm)
    EQ <- "Reverse Monomolecular"
    
  } else if (equation == 3) {
    r <- Monomolecular(R, parm)
    EQ <- "Monomolecular"
    
  } else if (equation == 7) {
    r <- Inv.Monomolecular(R, parm)
    EQ <- "Inverse Monomolecular"
    
  } else if (equation == 8) {
    r <- Inv.Ricker(R, parm)
    EQ <- "Inverse Ricker"
    
  } else if (equation == 4) {
    r <- Ricker(R, parm)
    EQ <- "Ricker"
    
  } else if (equation == 6) {
    r <- Rev.Ricker(R, parm)
    EQ <- "Reverse Ricker"
    
  } else if (equation == 2) {
    r <- Inv.Rev.Ricker(R, parm)
    EQ <- "Inverse-Reverse Ricker"
    
  } else {
    r <- (R * 0) + 1 #  Distance
    EQ <- "Distance"
  } # End if-else
  File.name <- NAME
  
  if (cellStats(r, "max") > 1e6)
    r <-
    SCALE(r, 1, 1e6) # Rescale surface in case resistance are too high
  
  if (!is.null(out)) {
    writeRaster(
      x = r,
      filename = paste0(EXPORT.dir, File.name, ".asc"),
      overwrite = TRUE
    )
    return(r)
  } else {
    return(r)
  }
}

#### ORIGINAL ####
#' #' Apply transformation to continuous resistance surface
#' #'
#' #' Apply one the eight resistance transformations to a continuous resistance surface
#' #'
#' #' @param transformation Transformation equation to apply. Can be provided as the name of the transformation or its numeric equivalent (see details)
#' #' @param shape Value of the shape parameter
#' #' @param max Value of the maximum value parameter
#' #' @param r Resistance surface to be transformed. Can be supplied as full path to .asc file or as a raster object
#' #' @param out Directory to write transformed .asc file. Default is NULL, and will not export .asc file
#' #' @usage Resistance.tran(transformation, shape, max, r, out)
#' #' @return R raster object
#' #' @details Valid arguements for \code{transformation} are:\cr
#' #' \tabular{ll}{
#' #'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#' #'    \tab 2 = "Inverse-Reverse Ricker"\cr
#' #'    \tab 3 = "Monomolecular"\cr
#' #'    \tab 4 = "Ricker"\cr
#' #'    \tab 5 = "Reverse Monomolecular"\cr
#' #'    \tab 6 = "Reverse Ricker"\cr
#' #'    \tab 7 = "Inverse Monomolecular"\cr
#' #'    \tab 8 = "Inverse Ricker"\cr
#' #'    \tab 9 = "Distance"\cr
#' #'    }
#' #'
#' #' The Distance transformation sets all values equal to one. Because of the flexibility of the Ricker function to take a monomolecular shape (try \code{Plot.trans(PARM=c(10,100), Resistance=c(1,10), transformation="Ricker")} to see this), whenever a shape parameter >6 is selected in combination with a Ricker family transformation, the transformation reverts to a Distance transformation. In general, it seems that using a combination of intermediate Ricker and Monomolecular transformations provides the best, most flexible coverasge of parameter space.
#' #' @export
#' #' @author Bill Peterman <Bill.Peterman@@gmail.com>
#'
#' Resistance.tran <- function(transformation, shape, max, r, out=NULL){
#'   if(class(r)[1]!='RasterLayer') {
#'     R<-raster(r)
#'     NAME <- basename(r)
#'     NAME<-sub("^([^.]*).*", "\\1", NAME)
#'     names(R)<-NAME
#'   } else {
#'     R<-r
#'     NAME <- r@data@names
#'   }
#'   if(is.numeric(transformation)){
#'     parm<-c(transformation,shape, max)
#'   } else {
#'     parm<-c(get.EQ(transformation),shape, max)
#'   }
#'
#'   EXPORT.dir<-out
#'
#'   ######
#'
#'   r <-SCALE(data=R,MIN=0,MAX=10)
#'
#'   # Set equation for continuous surface
#'   equation <- floor(parm[1]) # Parameter can range from 1-9.99
#'
#'   #       # Read in resistance surface to be optimized
#'   #       SHAPE <-  (parm[2])
#'   #       Max.SCALE <- (parm[3])
#'
#'
#'   ##################
#'
#'   # Apply specified transformation
#'   if(equation==1){
#'     r <- Inv.Rev.Monomolecular(r,parm)
#'     EQ <- "Inverse-Reverse Monomolecular"
#'
#'   } else if(equation==5){
#'     r <- Rev.Monomolecular(r,parm)
#'     EQ <- "Reverse Monomolecular"
#'
#'   } else if(equation==3){
#'     r <- Monomolecular(r,parm)
#'     EQ <- "Monomolecular"
#'
#'   } else if (equation==7) {
#'     r <- Inv.Monomolecular(r,parm)
#'     EQ <- "Inverse Monomolecular"
#'
#'   } else if (equation==8) {
#'     r <- Inv.Ricker(r,parm)
#'     EQ <- "Inverse Ricker"
#'
#'   } else if (equation==4) {
#'     r <- Ricker(r,parm)
#'     EQ <- "Ricker"
#'
#'   } else if (equation==6) {
#'     r <- Rev.Ricker(r,parm)
#'     EQ <- "Reverse Ricker"
#'
#'   } else if (equation==2) {
#'     r <- Inv.Rev.Ricker(r,parm)
#'     EQ <- "Inverse-Reverse Ricker"
#'
#'   } else {
#'     r <- (r*0)+1 #  Distance
#'     EQ <- "Distance"
#'   } # End if-else
#'   File.name <- NAME
#'
#'   if(cellStats(r,"max")>1e6)  r<-SCALE(r,1,1e6) # Rescale surface in case resistance are too high
#'
#'   if(!is.null(out)){
#'     writeRaster(x=r,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
#'     return(r)
#'   } else {
#'     return(r)
#'   }
#' }