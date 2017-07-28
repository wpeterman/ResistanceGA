#' Get cost distance using gdistance
#'
#' Execute gdistance
#'
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param r Accepts two types of inputs. Provide either the path to the resistance surface file (.asc) or specify an R RasterLayer object
#' @return A costDistance matrix object from gdistance
#' @usage Run_gdistance(gdist.inputs, r)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Run_gdistance <- function(gdist.inputs, r) {
  if (class(r)[1] != 'RasterLayer') {
    r <- raster(r)
  }
  tr <- transition(
    x = r,
    transitionFunction = gdist.inputs$transitionFunction,
    directions = gdist.inputs$directions
  )
  
  if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8) {
    tr <- geoCorrection(tr, "c")
  }
  
  if (gdist.inputs$method == 'costDistance') {
    ## Calculate least cost distance (Default)
    ret <- costDistance(tr, gdist.inputs$samples)
  } else {
    ret <-
      commuteDistance(tr, gdist.inputs$samples) / 1000   ## Calculate commute time = resistance distance
  }
  return(ret)
}