#' Prepare data for optimization using \code{gdistance}
#'
#' Creates a necessary input for optimizing resistance surfaces based on pairwise cost distances, implemented using the \code{gdistance} library
#'
#' @param n.Pops The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param samples Either provide the path to a .txt file containing the xy coordinates, or provide a matrix with x values in column 1 and y values in column 2. Alternatively, you can provide a \code{\link[sp]{SpatialPoints}} object
#' @param covariates Data frame of additional covariates that you want included in the MLPE model during opitmization.
#' @param formula If covariates are included in the model, specify the R formula for the fixed effects portion of the MLPE model.
#' @param transitionFunction The function to calculate the gdistance TransitionLayer object. See \code{\link[gdistance]{transition}}. Default = function(x) 1/mean(x)
#' @param directions Directions in which cells are connected (4, 8, 16, or other). Default = 8
#' @param longlat Logical. If true, a \code{\link[gdistance]{geoCorrection}} will be applied to the transition  matrix. Default = FALSE
#' @param method Specify whether pairwise distance should be calulated using the \code{\link[gdistance]{costDistance}} or \code{\link[gdistance]{commuteDistance}} (Default) functions. \code{\link[gdistance]{costDistance}} calculates least cost path distance, \code{\link[gdistance]{commuteDistance}} is equivalent (i.e. nearly perfectly correlated with) resistance distance calculated by CIRCUITSCAPE.
#' @param min.max_dist NOT YET SUPPORTED. Optional. Specify the minimum and maximum distance at which pairwise comparisons will be made(e.g., c(1, 50)). Euclidean distances below and above the minumum and maximum values will be omitted from the analysis. This has potential to reduce analysis time, but also reduces the number of pairwise comparisons.
#' @param keep NOT YET SUPPORTED. An optional vector equal to length \code{response} (i.e. all pairwise observations), with 1 indicating to keep the observation, and 0 to drop the observation. This can be used in conjunction with, or in place of \code{min.max_dist} to select which observations to include in analyses.
#' @return An R object that is a required input into optimization functions
#' @details When specifying a formula, provide it as: \code{response ~ covariate}.
#' the formula \code{response} will use the vector of values specified for the \code{response} parameter. Make sure that covariate names match variable names provided in \code{covariates}

#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>

#' @usage gdist.prep(n.Pops, 
#'                   response = NULL,
#'                   samples,
#'                   covariates = NULL,
#'                   formula = NULL,
#'                   transitionFunction = function(x)  1 / mean(x),
#'                   directions = 8,
#'                   longlat = FALSE,
#'                   method = 'commuteDistance',
#'                   min.max_dist = NULL,
#'                   keep = NULL)
#' 
#' @examples  
#' ## Not run:
#' ## *** TO BE COMPLETED *** ##
#' 
#' ## End (Not run)

gdist.prep <-
  function(n.Pops,
           response = NULL,
           samples,
           covariates = NULL,
           formula = NULL,
           transitionFunction = function(x)
             1 / mean(x),
           directions = 8,
           longlat = FALSE,
           method = 'commuteDistance',
           min.max_dist = NULL,
           keep = NULL) {
    
    if (method != 'commuteDistance') {
      method <- 'costDistance'
    }
    
    if (!is.null(response)) {
      TEST.response <- is.vector(response)
      if (TEST.response == FALSE) {
        stop("The object 'response' is not in the form of a single column vector")
      }
    }
    
    # Checks ------------------------------------------------------------------
    
    # Ensure covariates have same length as response
    if(!is.null(covariates) && !is.data.frame(covariates)) {
      stop("Please provide a data frame when specifying additional covariates")
    } 
    
    if(!is.null(covariates) && nrow(covariates) != length(response)) {
      stop("Response and covariates must have the same number of observations")
    } 
    
    if (class(samples)[1] == 'matrix') {
      if (ncol(samples) > 2) {
        stop("The specified matrix with xy coordinates has too many columns")
      }
      sp <- SpatialPoints(samples)
    } else if (class(samples)[1] == 'SpatialPoints') {
      sp <- samples
    } else {
      #   grepl(".txt", x = samples)){
      if (!file.exists(samples)) {
        stop("The path to the specified samples.txt file is incorrect")
      }
      sp <- SpatialPoints(read.delim(samples, header = F)[,-1])
      
    }
    
    if (n.Pops != length(sp)) {
      stop("n.Pops does not equal the number of sample locations")
    }

    if(!is.null(keep)) {
      ID <- To.From.ID(n.Pops)
      
      if(length(keep) != nrow(ID)) {
        stop("`keep` vector must be same length as the number of pairwise combinations")
      }
      
      drop.obs <- keep 
      
      ZZ <- ZZ.mat(ID, drop.obs)
      # ZZ <- ZZ.mat_select(ID, drop.obs)
      
    } else {
      ID <- To.From.ID(n.Pops)
      ZZ <- ZZ.mat(ID)
    }
    
    df <- NULL
    if(!is.null(response)) {
      if(!is.null(covariates)) {
        df <- data.frame(gd = response,
                         covariates,
                         pop = ID$pop1)
      } else {
        df <- data.frame(gd = response,
                         pop = ID$pop1)
      }
      
      if(!is.null(formula)) {
        formula <- update(formula, gd ~ . + cd + (1 | pop))
      } else {
        formula <- gd ~ cd + (1 | pop)
      }
      
      ## Check slope
      m <- lm(gd ~ c(dist(sp@coords)), data = df)
      if(coef(m)[2] < 0){
        warning('Genetic distance decreases with distance. This is likely to result in a failed optimization.\nCheck your measure carefully and consider subtracting your values from 1 to reverse the relationship.')
      }
    }
    

    
    
    (
      ret <-
        list(
          response = response,
          samples = sp,
          covariates = covariates,
          formula = formula,
          transitionFunction = transitionFunction,
          directions = directions,
          ID = ID,
          ZZ = ZZ,
          keep = keep,
          n.Pops = n.Pops,
          longlat = longlat,
          method = method,
          df = df
        )
    )
    # Min-Max Distance ------------------------------------------------------------
    
    if(!is.null(min.max_dist)) {
      stop("This feature is not yet supported")
      if (length(min.max_dist) != 2) {
        stop("Specify 'min.max_dist' as a 2-element vector: c(min, max)")
      }
      
      toMatch <- min.max_dist
      
      # Function to make list of observations to include
      Min.MAX <- as.matrix(dist(samples@coords))
      Min.MAX[upper.tri(Min.MAX)] <- NA
      
      drop.pairs <- which((Min.MAX < toMatch[1]) | (Min.MAX > toMatch[2]), arr.ind = T)
      
      keep.pairs <- which((Min.MAX > toMatch[1]) & (Min.MAX < toMatch[2]), arr.ind = T)
      
      Min.MAX[cbind(drop.pairs[,1], drop.pairs[,2])] <- NA
      
      drop.obs <- lower(Min.MAX)
      drop.obs <- ifelse(is.na(drop.obs), 0, 1)
      
      ID <- To.From.ID(n.Pops)
      
      if(!is.null(keep)) {
        if(length(keep) != nrow(ID)) {
          stop("`keep` vector must be same length as the number of pairwise combinations")
        }
        
        drop.obs <- (keep * drop.obs) + keep
        drop.obs <- ifelse(drop.obs >= 1, 1, 0)
        
      }
      
      # keep.dat <- lower(Min.MAX)
      # ID <- data.frame(pop1 = c(keep.pairs[,2], drop.pairs[,2]), pop2 = c(keep.pairs[,1], drop.pairs[,1]))
      # drop.row <- which(ID[,1] == ID[,2])
      # ID <- ID[-drop.row,]
      # swap <- which(table(ID$pop1) > 1)
      # n1 <- table(ID$pop1)[[swap[[1]]]]
      # p1 <- ID[n1, 1]
      # p2 <- ID[n1, 2]
      # ID[n1, 1] <- p2
      # ID[n1, 2] <- p1
      # ID.n <- ID
      # ID$pop1 <- factor(ID$pop1)
      # ID$pop2 <- factor(ID$pop2)
      
      # ZZ <- ZZ.mat_select(ID, drop.obs)
      ZZ <- ZZ.mat(ID, drop.obs)
      
      # mx.dist <- as.matrix(dist(samples@coords))
      # mx.dist[mx.dist > max.dist] <- NA
      # mx.dist[!is.na(mx.dist)] <- 1
      # mx.vec <- lower(mx.dist)
      
      #   dat <- data.frame(gd = response,
      #                     ed = lower(Min.MAX),
      #                     pop = ID$pop1)
      #   
      #   dat <- dat[complete.cases(dat),]
      #   
      #   mod <- mlpe_rga(gd ~ ed + (1|pop), data = dat, ZZ = ZZ)
      # }
      
      if(!is.null(covariates)) {
        df <- data.frame(gd = response,
                         covariates,
                         pop = ID$pop1)
      } else {
        df <- NULL
      }
      
      if(!is.null(formula)) {
        formula <- update(formula, gd ~ . + cd + (1 | pop))
      }
      
      (
        ret <-
          list(
            response = response,
            samples = sp,
            covariates = covariates,
            formula = formula,
            transitionFunction = transitionFunction,
            directions = directions,
            ID = ID,
            ZZ = ZZ,
            keep = drop.obs,
            n.Pops = n.Pops,
            longlat = longlat,
            method = method,
            df = df,
            opt.input = 'gdist'
          )
      )
    }
    return(ret)
  }
