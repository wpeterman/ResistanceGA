#' Get cost distance using gdistance
#' Execute gdistance
#'
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function
#' @param r Accepts two types of inputs. Provide either the path to the resistance surface file (.asc) or specify an R RasterLayer object
#' @param scl scale the correction values (default is TRUE). No scaling should be done if the user wants to obtain absolute distance values as output. See \code{\link[gdistance]{geoCorrection}} for details
#' @return A costDistance matrix object from gdistance
#' @usage Run_gdistance(gdist.inputs, r, scl)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Run_gdistance <- function(gdist.inputs, 
                          r, 
                          scl = TRUE) {
  out <- tryCatch(
    {
      
      if (class(r)[1] != 'RasterLayer') {
        r <- raster(r)
      }
      
      tr <- transition(
        x = r,
        transitionFunction = gdist.inputs$transitionFunction,
        directions = gdist.inputs$directions
      )
      if(is.null(gdist.inputs$keep)) {
        if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'costDistance') {
          trC <- geoCorrection(tr, "c", scl = scl)
          ret <- costDistance(trC, gdist.inputs$samples)
          rm(trC, r)
          gc()
          
        } # End costDistance
        
        if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'commuteDistance') {
          
          trR <- geoCorrection(tr, "r", scl = scl)
          ret <- commuteDistance(trR, gdist.inputs$samples) / 1000
          rm(trR, r)
          gc()
        } # End commuteDistance
        
        if(!exists('ret')) {
          if(gdist.inputs$method == 'commuteDistance') {
            trR <- geoCorrection(tr, "r", scl = scl)
            ret <- commuteDistance(trR, gdist.inputs$samples) / 1000
            rm(trR, tr, r)
            
          } else {
            trC <- geoCorrection(tr, "c", scl = scl)
            ret <- costDistance(trC, gdist.inputs$samples)
            rm(trC, tr, r)
          }
          
          gc()
        }
      } else { # Run on select pairs
        if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'costDistance') {
          trC <- geoCorrection(tr, "c", scl = scl)
          
          ret <- vector(mode = 'numeric', length = sum(gdist.inputs$keep))
          kp <- which(gdist.inputs$keep == 1)
          id <- gdist.inputs$ID
          id$pop1 <- as.numeric(id$pop1)
          id$pop2 <- as.numeric(id$pop2)
          id <- as.matrix(id)
          c <- 0
          for(i in kp) {
            c <- c + 1
            pts <- as.vector(id[i,])
            ret[c] <- costDistance(trC, gdist.inputs$samples[pts])
          }
          
          rm(trC, r)
          gc()
          
        } # End costDistance
        
        if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'commuteDistance') {
          
          trR <- geoCorrection(tr, "r", scl = scl)
          
          ret <- vector(mode = 'numeric', length = sum(gdist.inputs$keep))
          kp <- which(gdist.inputs$keep == 1)
          id <- gdist.inputs$ID
          id$pop1 <- as.numeric(id$pop1)
          id$pop2 <- as.numeric(id$pop2)
          id <- as.matrix(id)
          c <- 0
          for(i in kp) {
            c <- c + 1
            pts <- as.vector(id[i,])
            ret[c] <- commuteDistance(trR, gdist.inputs$samples[pts]) / 1000
          }
          
          rm(trR, r)
          gc()
        } # End commuteDistance
      }

      
      return(ret)
    },
    warning = function(warning_condition) {
      # message(warning_condition)
      return(-99999)
    }, 
    error = function(error_condition) {
      # message(error_condition)
      return(-99999)
    }
    
    # tr <- transition(
    #   x = r,
    #   transitionFunction = gdist.inputs$transitionFunction,
    #   directions = gdist.inputs$directions
    # )
    # 
    # if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'costDistance') {
    #   trC <- geoCorrection(tr, "c", scl = scl)
    #   ret <- try(costDistance(trC, gdist.inputs$samples), silent = TRUE)
    #   rm(trC, r)
    #   
    #   if(isTRUE(ret) == 'try-error') {
    #     ret <- -99999
    #     return(ret)
    #   } else {
    #     return(ret)
    #   }
    # } # End costDistance
    # 
    # if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'commuteDistance') {
    #   
    #   trR <- geoCorrection(tr, "r", scl = scl)
    #   ret <- try(commuteDistance(trR, gdist.inputs$samples), silent = TRUE) 
    #   rm(trR, r)
    #   
    #   if(isTRUE(ret) == 'try-error') {
    #     ret <- -99999
    #     return(ret)
    #   } else if (anyNA(ret)) {
    #     ret <- -99999
    #     return(ret)
    #   } else {
    #     ret <- ret / 1000
    #     return(ret)
    #   }
    #   
    # } # End commuteDistance 
  ) # End tryCatch
}