#' Kernel smoothing function
#'
#' Apply Gaussian kernel smoothing of specified sigma
#'
#' @param raster A RasterLayer object to be smoothed
#' @param sigma The standard deviation of the Gaussian smoothing parameter (see \code{\link[spatstat]{blur}} documentation in the \code{spatstat} package.)
#' @param SCALE Logical. Should the smoothed raster surface be scaled to range from 0-10 (Default = FALSE)

#' @usage k.smooth (raster, sigma, SCALE)
#' 
#' @details The sigma parameter indicates the standard deviation of the Gaussian smoothing function. Note that sigma is in raster cells, not spatial units.


#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
k.smooth <- function(raster,
                     sigma,
                     SCALE = FALSE) {
  r.mat <- raster::as.matrix(raster)
  x <- spatstat::as.im(r.mat)
  
  x.blur <- spatstat::blur(x = x,
                           sigma = sigma,
                           normalise = TRUE,
                           bleed = FALSE)
  
  
  values(raster) <- spatstat::as.matrix.im(x.blur)
  
  
  # Original using smoothie -------------------------------------------------
  # zmat <- as.matrix(raster)
  # f <- smoothie::kernel2dsmooth(
  #   zmat,
  #   kernel.type = "gauss",
  #   nx = nrow(raster),
  #   ny = ncol(raster),
  #   sigma = sigma
  # )
  # 
  # values(raster) <- f
  
  if(SCALE == TRUE) {
    
    raster <- SCALE(raster, 0, 10)
    return(raster)
    
  } else {
    
    return(raster)
    
  }
}