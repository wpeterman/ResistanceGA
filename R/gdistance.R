#' Prepare data for optimization using \code{gdistance}
#' 
#' Creates a necessary input for optimizing resistance surfaces based on pairwise cost distances, implemented using the \code{gdistance} library
#' 
#' @param n.Pops The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param samples Either provide the path to a .txt file containing the xy coordinates, or provide a matrix with x values in column 1 and y values in column 2. Alternatively, you can provide a \code{\link[sp]{SpatialPoints}} object
#' @param transitionFunction The function to calculate the gdistance TransitionLayer object. See \code{\link[gdistance]{transition}}. Default = function(x) 1/mean(x)
#' @param directions Directions in which cells are connected (4, 8, 16, or other). Default = 8
#' @param longlat Logical. If true, a \code{\link[gdistance]{geoCorrection}} ill be applied to the transition  matrix. Defaault = FALSE
#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage gdist.prep(n.Pops, response, samples, transitionFunction, directions, longlat)
#' 

gdist.prep <- function(n.Pops, response=NULL, samples, transitionFunction=function(x) 1/mean(x), directions=8, longlat=FALSE){
  if(!is.null(response)) {TEST.response <- is.vector(response)
                          if(TEST.response==FALSE) {stop("The object 'response' is not in the form of a single column vector")}}
   
  if(class(samples)[1]=='matrix') {
    if(ncol(samples)>2) {stop("The specified matrix with xy coordinates has too many columns")}
    sp <- SpatialPoints(samples)   
  } else if(class(samples)[1]=='SpatialPoints') { 
    sp <- samples
  } else {
#   grepl(".txt", x = samples)){
    if(!file.exists(samples)) {stop("The path to the specified samples.txt file is incorrect")}
    sp <- SpatialPoints(read.delim(samples,header = F)[,-1])   
  } 
  
  if(n.Pops!=length(sp)) {stop("n.Pops does not equal the number of sample locations")}
  
  ID<-To.From.ID(n.Pops)
  ZZ<-ZZ.mat(ID)
  
  (ret<-list(response=response, samples=sp, transitionFunction=transitionFunction, directions=directions, ID=ID, ZZ=ZZ, n.Pops=n.Pops, longlat=longlat))
}


