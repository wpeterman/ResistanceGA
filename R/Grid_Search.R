#' Conduct grid search of response surface 
#' 
#' Visualize the AICc response surface
#' 
#' @param shape A vector of values for the shape parameter
#' @param max A vector of values for the maximum value parameter
#' @param transformation Transformation to apply. Can be either numeric or character of transformation name
#' @param Resistance An R Raster object, or path to a .asc file
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @usage Grid.Search(shape, max, transformation, Resistance, CS.inputs, gdist.inputs, GA.inputs)
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @return This function will return a \code{filled.contour} plot. Additionally, an object with values that can be plotted with \code{filled.contour} to visualize the response surface
#' @details This function will perform a full factorial grid search of the values provided in the shape and max.scale vectors. Depending on the number of values provided for each, and the time it takes to run each iteration, this process may take a while to complete. \cr Suitable values for transformation:\cr
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



Grid.Search <- function(shape, max, transformation, Resistance, CS.inputs=NULL, gdist.inputs=NULL, GA.inputs) {
  if(class(Resistance)[1]!='RasterLayer') {  
    r <- raster(Resistance)
    r <- SCALE(r,0,10)
  } else {    
    r <-SCALE(Resistance,0,10)
  }
  
  GRID <- expand.grid(shape,max)
  RESULTS <- matrix(nrow=nrow(GRID),ncol=3); colnames(RESULTS)<-c("shape","max","AICc")
  if(!is.numeric(transformation)){
    EQ<-get.EQ(transformation)
  } else {
    EQ <- transformation
  }
  
  if(!is.null(CS.inputs)){
    for(i in 1:nrow(GRID)){
      # Modified 16 September 2015
      AICc <- Resistance.Opt_AICc(PARM = c(EQ,t(GRID[i,])),
                                  Resistance = r,
                                  CS.inputs = CS.inputs, 
                                  Min.Max='min',
                                  GA.inputs = GA.inputs)
      
      # Original function used
      # AICc<-Resistance.Optimization_cont.nlm(PARM=log(c(t(GRID[i,]))),Resistance=r,equation=EQ, get.best=FALSE,CS.inputs,Min.Max='min',write.dir=write.dir)
      
      results<-as.matrix(cbind(GRID[i,],AICc))
      
      RESULTS[i,]<-results  
    }
  } else {
    for(i in 1:nrow(GRID)){
      # Modified 16 September 2015
      AICc <- Resistance.Opt_AICc(PARM = c(EQ,t(GRID[i,])),
                                  Resistance = r,
                                  gdist.inputs = gdist.inputs, 
                                  Min.Max='min',
                                  GA.inputs = GA.inputs)
      
      # Original function used
      # AICc<-Resistance.Optimization_cont.nlm(PARM=log(c(t(GRID[i,]))),Resistance=r,equation=EQ, get.best=FALSE,CS.inputs,Min.Max='min',write.dir=write.dir)
      
      results<-as.matrix(cbind(GRID[i,],AICc))
      
      RESULTS[i,]<-results  
    }
  }
  RESULTS <- data.frame(RESULTS)
  Results.mat <- interp(RESULTS$shape,RESULTS$max,RESULTS$AICc,duplicate='strip')
  filled.contour(Results.mat,col=topo.colors(20),xlab="Shape parameter",ylab="Maximum value parameter")
  
  AICc<-RESULTS
  colnames(AICc)<-c("shape","max","AICc")
  Results.mat<-list(Plot.data=Results.mat,AICc=AICc)
  
  return(Results.mat)
}