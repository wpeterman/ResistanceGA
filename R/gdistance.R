#' Prepare data for optimization using \code{gdistance}
#' 
#' Creates a necessary input for optimizing resistance surfaces based on pairwise cost distances, implemented using the \code{gdistance} library
#' 
#' @param n.POPS The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param samples Either provide the path to a .txt file containing the xy coordinates, or provide a matrix with x values in column 1 and y values in column 2.
#' @param transitionFunction The function to calculate the gdistance TransitionLayer object. See \code{\link[gdistance]{transition}}. Default = function(x) 1/mean(x)
#' @param directions Directions in which cells are connected (4, 8, 16, or other). Default = 8
#' @param longlat Logical. If true, a \code{\link[gdistance]{geoCorection}} ill be applied to the transition  matrix. Defaault = FALSE
#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage gdist.prep(n.POPS, response, samples, transitionFunction, directions, longlat)
#' 

gdist.prep <- function(n.POPS, response=NULL, samples, transitionFunction=function(x) 1/mean(x), directions=8, longlat=FALSE){
  ID<-To.From.ID(n.POPS)
  ZZ<-ZZ.mat(ID)
  
  (ret<-list(response=response, samples=samples, transitionFunction=transitionFunction, directions=directions, ID=ID, ZZ=ZZ, n.Pops=n.POPS, longlat=longlat))
}


