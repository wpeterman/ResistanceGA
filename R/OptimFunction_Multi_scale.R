#' Optimize multiple resistance surfaces simultaneously with kernel smoothing
#'
#' Create composite resistance surface by simultaneously optimizing multiple binary and continuous surfaces to optimized with a kernel smoothing parameter. This optimization function is designed to be called from GA
#'
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Should be a vector with parameters specified in the order of resistance surfaces. These values are selected during optimization if called within GA function.
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param jl.inputs Object created from running \code{\link[ResistanceGA]{jl.prep}} function. Defined if optimizing using CIRCUITSCAPE run in Julia
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max')
#' @param quiet Logical, if TRUE, AICc and iteration time will not be printed to the screen at the completion of each iteration. Default = FALSE
#' @return Objective function value (either AIC, R2, or LL) from mixed effect model
#' @noRd
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Resistance.Opt_multi.scale <-
  function(PARM,
           CS.inputs = NULL,
           gdist.inputs = NULL,
           jl.inputs = NULL,
           GA.inputs,
           Min.Max,
           quiet = FALSE) {
    if (is.null(GA.inputs$scale)) {
      stop(
        "This function should only be used if you intend to apply kernel smoothing to your resistance surfaces"
      )
    }
    
    t1 <- proc.time()[3]
    
    method <- GA.inputs$method
    EXPORT.dir <- GA.inputs$Write.dir
    ######
    #   r <- GA.inputs$Resistance.stack
    File.name = "resist_surface"
    
    
    # Run CS ------------------------------------------------------------------
    
    
    if (!is.null(CS.inputs)) {
      r <-
        Combine_Surfaces(
          PARM = PARM,
          CS.inputs = CS.inputs,
          GA.inputs = GA.inputs,
          out = NULL,
          File.name = File.name,
          rescale = FALSE
        )
      
      # if(cellStats(r, "mean") == 0) { # Skip iteration
      if(mean(r@data@values, na.rm = TRUE) == 0) { # Skip iteration        
        obj.func.opt <- -99999
        
      }
      
      if(!exists('obj.func.opt')) {
        CS.resist <- try(Run_CS(
          CS.inputs,
          r = r
        ), TRUE)
      }
      
      if(exists('CS.resist') && isTRUE(class(CS.resist) == 'try-error')) {
        obj.func.opt <- -99999
      } 
      
      if(exists('CS.resist') && isTRUE(class(CS.resist) != 'try-error')) { # Continue with iteration
        
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
    
    
    
    # Run gdistance -----------------------------------------------------------
    
    
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
      
      # if(cellStats(r, "mean") == 0) { # Skip iteration
      if(mean(r@data@values, na.rm = TRUE) == 0) { # Skip iteration        
        obj.func.opt <- -99999
        
      } 
      
      if(!exists('obj.func.opt')) {
        cd <- try(Run_gdistance(gdist.inputs, r), TRUE)
      }
      
      if(exists('cd') && isTRUE(class(cd) == 'try-error')) {
        obj.func.opt <- -99999
        rm(cd, r)
        gc()
      } 
      
      if(exists('cd') && isTRUE(class(cd) != 'try-error')) { 
        
        # Continue with iteration
        
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
        rm(cd, r)
        gc()
      }
    }
    
    # Julia -----------------------------------------------------------
    
    
    if (!is.null(jl.inputs)) {
      
      r <-
        Combine_Surfaces(
          PARM = PARM,
          gdist.inputs = gdist.inputs,
          GA.inputs = GA.inputs,
          out = NULL,
          File.name = File.name,
          rescale = FALSE
        )
      
      # if(cellStats(r, "mean") == 0) { # Skip iteration
      if(mean(r@data@values, na.rm = TRUE) == 0) { # Skip iteration        
        obj.func.opt <- -99999
        
      }
      
      if(!exists('obj.func.opt')) {
        cd <- try(Run_CS.jl(jl.inputs, r), TRUE)
      }
      
      if(exists('cd') && isTRUE(class(cd) == 'try-error')) {
        obj.func.opt <- -99999
      } 
      
      if(exists('cd') && isTRUE(class(cd) != 'try-error')) { # Continue with iteration
        
        
        if (method == "AIC") {
          obj.func <- suppressWarnings(AIC(
            MLPE.lmm2(
              resistance = cd,
              response = jl.inputs$response,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ,
              REML = FALSE
            )
          ))
          obj.func.opt <- obj.func * -1
        } else if (method == "R2") {
          obj.func <- suppressWarnings(r.squaredGLMM(
            MLPE.lmm2(
              resistance = cd,
              response =
                jl.inputs$response,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ,
              REML = FALSE
            )
          ))
          obj.func.opt <- obj.func[[1]]
        } else {
          obj.func <- suppressWarnings(logLik(
            MLPE.lmm2(
              resistance = cd,
              response = jl.inputs$response,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ,
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
    
    if(!is.null(GA.inputs$opt.digits)) {
      obj.func.opt <- round(obj.func.opt, GA.inputs$opt.digits)
      return(obj.func.opt)
    } else {
      return(obj.func.opt)
    }
  }