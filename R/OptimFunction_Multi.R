#' Optimize multiple resistance surfaces simultaneously
#'
#' Create composite resistance surface by simultaneously optimizing multiple categoricla and continuous surfaces. This optimization function is designed to be called from GA
#'
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Should be a vector with parameters specified in the order of resistance surfaces. These values are selected during optimization if called within GA function.
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max')
#' @param quiet Logical, if TRUE, AICc and iteration time will not be printed to the screen at the completion of each iteration. Default = FALSE
#' @return AIC value from mixed effect model
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Resistance.Opt_multi <- function(PARM,
                                 CS.inputs = NULL,
                                 gdist.inputs = NULL,
                                 GA.inputs,
                                 Min.Max,
                                 quiet = FALSE) {
  t1 <- proc.time()[3]
  
  method <- GA.inputs$method
  EXPORT.dir <- GA.inputs$Write.dir
  ######
  #   r <- GA.inputs$Resistance.stack
  File.name = "resist_surface"
  # CIRCUITSCAPE ------------------------------------------------------------
  
  if (!is.null(CS.inputs)) {
    Combine_Surfaces(
      PARM = PARM,
      CS.inputs = CS.inputs,
      GA.inputs = GA.inputs,
      out = GA.inputs$Write.dir,
      File.name = File.name,
      rescale = FALSE
    )
    
    r <- raster(paste0(EXPORT.dir, File.name, ".asc"))
    if(cellStats(r, "mean") == 0) { # Skip iteration
      
      obj.func.opt <- -99999
      
    } else { # Continue with iteration
      CS.resist <-
        Run_CS2(
          CS.inputs,
          GA.inputs,
          r = multi_surface,
          EXPORT.dir = GA.inputs$Write.dir,
          File.name = File.name
        )
      
      # Replace NA with 0...a workaround for errors when two points fall within the same cell.
      # CS.resist[is.na(CS.resist)] <- 0
      
      # Run mixed effect model on each Circuitscape effective resistance
      if (method == "AIC") {
        obj.func <- suppressWarnings(AIC(
          MLPE.lmm2(
            resistance = CS.resist,
            response = CS.inputs$response,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ,
            REML = FALSE
          )
        ))
        obj.func.opt <- obj.func * -1
      } else if (method == "R2") {
        obj.func <-
          suppressWarnings(r.squaredGLMM(
            MLPE.lmm2(
              resistance = CS.resist,
              response = CS.inputs$response,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ,
              REML = FALSE
            )
          ))
        obj.func.opt <- obj.func[[1]]
      } else {
        obj.func <- suppressWarnings(logLik(
          MLPE.lmm2(
            resistance = CS.resist,
            response = CS.inputs$response,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ,
            REML = FALSE
          )
        ))
        
        obj.func.opt <- obj.func[[1]]
      }
    }
  }
  
  
  # gdistance ---------------------------------------------------------------
  
  if (!is.null(gdist.inputs)) {
    r <-
      Combine_Surfaces(
        PARM = PARM,
        gdist.inputs = gdist.inputs,
        GA.inputs = GA.inputs,
        out = NULL,
        File.name = File.name,
        rescale = FALSE
      )
    
    if(cellStats(r, "mean") == 0) { # Skip iteration
      
      obj.func.opt <- -99999
      
    } 
    
    cd <- try(Run_gdistance(gdist.inputs, r), TRUE)
    
    if(isTRUE(class(cd) == 'try-error') || isTRUE(exists('obj.func.opt'))) {
      
      obj.func.opt <- -99999
      
    } else { # Continue with iteration
      
      if (method == "AIC") {
        obj.func <- suppressWarnings(AIC(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ,
            REML = FALSE
          )
        ))
        obj.func.opt <- obj.func * -1
      } else if (method == "R2") {
        obj.func <- suppressWarnings(r.squaredGLMM(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ,
            REML = FALSE
          )
        ))
        obj.func.opt <- obj.func[[1]]
      } else {
        obj.func <- suppressWarnings(logLik(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ,
            REML = FALSE
          )
        ))
        
        obj.func.opt <- obj.func[[1]]
      }
    }
  }
  
  rt <- proc.time()[3] - t1
  if (quiet == FALSE) {
    cat(paste0("\t", "Iteration took ", round(rt, digits = 2), " seconds"), "\n")
    cat(paste0("\t", method, " = ", round(obj.func.opt, 4)), "\n", "\n")
  }
  
  return(obj.func.opt)
  # OPTIM.DIRECTION(Min.Max)*(obj.func) # Function to be minimized/maximized
  
}