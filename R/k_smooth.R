#' Kernel smoothing function
#'
#' Apply Gaussian kernel smoothing of specified sigma
#'
#' @param raster A RasterLayer object to be smoothed
#' @param sigma The standard deviation of the Gaussian smoothing parameter (see \code{\link[smoothie]{kernel2dsmooth}} documentation in the \code{smoothie} package.)

#' @usage k.smooth (raster, sigma)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
k.smooth <- function(raster,
                     sigma) {
  zmat <- as.matrix(raster)
  f <- smoothie::kernel2dsmooth(
    zmat,
    kernel.type = "gauss",
    nx = nrow(raster),
    ny = ncol(raster),
    sigma = sigma
  )
  
  values(raster) <- f
  
  raster <- SCALE(raster, 0, 10)
  return(raster)
}