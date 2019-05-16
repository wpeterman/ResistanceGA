#' Optimize resistance surfaces individually
#'
#' Optimize all resistance surfaces that are located in the same directory individually. This optimization function is designed to be called from GA
#'
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. A vector with parameters specified in the order of resistance surfaces.These values are selected during optimization if called within GA function.
#' @param Resistance Resistance surface to be optimized. This should be an R raster object. If not specified, the function will attempt to find the a resistance surface from \code{GA.inputs}
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param jl.inputs Object created from running \code{\link[ResistanceGA]{jl.prep}} function. Defined if optimizing using CIRCUITSCAPE run in Julia
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max'). Default in 'max'
#' @param iter A counter for the number of surfaces that will be optimized
#' @param quiet Logical, if TRUE AICc and iteration duration will not be printed to the screen at the completion of each iteration.
#' @return Objective function value (either AIC, R2, or LL) from mixed effect model
#' @noRd
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Resistance.Opt_single <-
  function(PARM,
           Resistance,
           CS.inputs = NULL,
           gdist.inputs = NULL,
           jl.inputs = NULL,
           GA.inputs,
           Min.Max = 'max',
           iter = NULL,
           quiet = FALSE) {
    t1 <- Sys.time()
    
    
    # Modify resistance surface -----------------------------------------------
    
    
    if(is.null(iter)) {
      iter <- 1
    }
    
    method <- GA.inputs$method
    EXPORT.dir <- GA.inputs$Write.dir
    select.trans <- GA.inputs$select.trans
    r <- Resistance
    keep <- 1 # Indicator to keep or drop transformation
    
    
    # * Categorical -----------------------------------------------------------
    
    
    if (GA.inputs$surface.type[iter] == "cat") {
      PARM <- PARM / min(PARM)
      parm <- PARM
      df <-
        data.frame(id = unique(r), PARM) # Data frame with original raster values and replacement values
      r <- subs(r, df)
      
    } else {
      # * Continuous -----------------------------------------------------------
      
      ## Else continuous surface
      r <- SCALE(r, 0, 10)
      
      # Set equation for continuous surface
      equation <- floor(PARM[1]) # Parameter can range from 1-9.99
      
      # Read in resistance surface to be optimized
      SHAPE <- (PARM[2])
      Max.SCALE <- (PARM[3])
      
      if (equation %in% select.trans[[iter]]) {
        # Apply specified transformation
        rick.eq <- (equation == 2 ||
                      equation == 4 || equation == 6 || equation == 8)
        if (rick.eq == TRUE & SHAPE > 5) {
          equation <- 9
          keep <- 0 # Drop transformation
        }
        
        if (equation == 1) {
          r <- Inv.Rev.Monomolecular(r, parm = PARM)
          EQ <- "Inverse-Reverse Monomolecular"
          
        } else if (equation == 5) {
          r <- Rev.Monomolecular(r, parm = PARM)
          EQ <- "Reverse Monomolecular"
          
        } else if (equation == 3) {
          r <- Monomolecular(r, parm = PARM)
          EQ <- "Monomolecular"
          
        } else if (equation == 7) {
          r <- Inv.Monomolecular(r, parm = PARM)
          EQ <- "Inverse Monomolecular"
          
        } else if (equation == 8) {
          r <- Inv.Ricker(r, parm = PARM)
          EQ <- "Inverse Ricker"
          
        } else if (equation == 4) {
          r <- Ricker(r, parm = PARM)
          EQ <- "Ricker"
          
        } else if (equation == 6) {
          r <- Rev.Ricker(r, parm = PARM)
          EQ <- "Reverse Ricker"
          
        } else if (equation == 2) {
          r <- Inv.Rev.Ricker(r, parm = PARM)
          EQ <- "Inverse-Reverse Ricker"
          
        } else {
          r <- (r * 0) + 1 #  Distance
          EQ <- "Distance"
        } # End if-else
      } # Close select transformation
    } # Close surface type if-else
    
    ## If a surface was reclassified or transformed, apply the following
    if ((GA.inputs$surface.type[iter] == "cat") ||
        (equation %in% select.trans[[iter]])) {
      
      if(keep == 0) {
        ## Use -99999 as 'ga' maximizes
        obj.func.opt <- -99999
      } else {
        
        rclmat <- matrix(c(1e-100, 1e-6, 1e-06, 1e6, Inf, 1e6, -1, 0, 1e6),
                         ncol = 3,
                         byrow = TRUE)
        
        r <- reclassify(r, rclmat)
        
        
        # Calculate resistance ----------------------------------------------------
        
        # * CS --------------------------------------------------------------------
        
        if (!is.null(CS.inputs)) {
          
          if(!exists('obj.func.opt')) {
            CS.resist <- try(Run_CS(CS.inputs, r), TRUE)
          }
          
          if(exists('CS.resist') && isTRUE(class(CS.resist) == 'try-error')) {
            obj.func.opt <- -99999
          } 
          
          if(exists('CS.resist') && isTRUE(class(CS.resist) != 'try-error')) {  # Continue with iteration          
            
            # writeRaster(
            #   x = r,
            #   filename = paste0(EXPORT.dir, File.name, ".asc"),
            #   overwrite = TRUE
            # )
            # CS.resist <-
            #   Run_CS2(
            #     CS.inputs,
            #     GA.inputs,
            #     r = r,
            #     # EXPORT.dir = GA.inputs$Write.dir,
            #     File.name = File.name
            #   )
            
            dat <- data.frame(gd = CS.inputs$response,
                              cd = scale(CS.resist),
                              pop = CS.inputs$ID$pop1)
            
            fit.mod <- mlpe_rga(formula = gd ~ cd + (1 | pop),
                                data = dat,
                                ZZ = CS.inputs$ZZ,
                                REML = FALSE)
            
            # Run mixed effect model on each Circuitscape effective resistance
            if(lme4::fixef(fit.mod)['cd'] < 0) {
              obj.func.opt <- -99999
            } else {
              if (method == "AIC") {
                obj.func <- suppressWarnings(AIC(
                  fit.mod
                ))
                obj.func.opt <- obj.func * -1
              } else if (method == "R2") {
                obj.func <-
                  suppressWarnings(r.squaredGLMM(
                    fit.mod
                  ))
                obj.func.opt <- obj.func[[1]]
                
              } else {
                obj.func <- suppressWarnings(logLik(
                  fit.mod
                ))
                obj.func.opt <- obj.func[[1]]
              }
            } # Positive parameter value
          } # End Obj func ifelse
          
          #   if (method == "AIC") {
          #     obj.func <- suppressWarnings(AIC(
          #       MLPE.lmm2(
          #         resistance = CS.resist,
          #         response = CS.inputs$response,
          #         ID = CS.inputs$ID,
          #         ZZ = CS.inputs$ZZ,
          #         REML = FALSE
          #       )
          #     ))
          #     obj.func.opt <- obj.func * -1
          #   } else if (method == "R2") {
          #     obj.func <-
          #       suppressWarnings(r.squaredGLMM(
          #         MLPE.lmm2(
          #           resistance = CS.resist,
          #           response =
          #             CS.inputs$response,
          #           ID = CS.inputs$ID,
          #           ZZ = CS.inputs$ZZ,
          #           REML = FALSE
          #         )
          #       ))
          #     obj.func.opt <- obj.func[[1]]
          #     
          #   } else {
          #     obj.func <- suppressWarnings(logLik(
          #       MLPE.lmm2(
          #         resistance = CS.resist,
          #         response = CS.inputs$response,
          #         ID = CS.inputs$ID,
          #         ZZ = CS.inputs$ZZ,
          #         REML = FALSE
          #       )
          #     ))
          #     obj.func.opt <- obj.func[[1]]
          #   }
          # } # End Obj func ifelse
        } # End CS Loop 
        
        # * gdistance ------------------------------------------------------------
        
        
        if (!is.null(gdist.inputs)) {
          
          # if(cellStats(r, "mean") == 0) { # Skip iteration
          if(mean(r@data@values, na.rm = TRUE) == 0) { # Skip iteration            
            obj.func.opt <- -99999
            
          } 
          
          if(!exists('obj.func.opt')) {
            cd <- try(Run_gdistance(gdist.inputs, r), TRUE)
          }
          
          if(exists('cd') && isTRUE(class(cd) == 'try-error')) {
            obj.func.opt <- -99999
            # rm(cd, r)
            gc()
          } 
          
          if(exists('cd') && isTRUE(class(cd) != 'try-error')) {  # Continue with iteration
            
            # l.cd <- as.vector(cd)
            # 
            # if(length(l.cd) != length(gdist.inputs$response)) {
            #   cd <- Run_gdistance(gdist.inputs, r)
            # }
            
            dat <- data.frame(gd = gdist.inputs$response,
                              cd = scale(c(cd)),
                              pop = gdist.inputs$ID$pop1) 
            
            fit.mod <-  mlpe_rga(formula = gd ~ cd + (1 | pop),
                                 data = dat,
                                 ZZ = gdist.inputs$ZZ,
                                 REML = FALSE)
            
            if(lme4::fixef(fit.mod)['cd'] < 0) {
              obj.func.opt <- -99999
            } else {
              if (method == "AIC") {
                obj.func <- suppressWarnings(AIC(
                  fit.mod
                ))
                obj.func.opt <- obj.func * -1
              } else if (method == "R2") {
                obj.func <-
                  suppressWarnings(r.squaredGLMM(
                    fit.mod
                  ))
                obj.func.opt <- obj.func[[1]]
                
              } else {
                obj.func <- suppressWarnings(logLik(
                  fit.mod
                ))
                obj.func.opt <- obj.func[[1]]
              }
            } # Positive parameter value
            
          } # End keep loop
          
          # if (method == "AIC") {
          #   obj.func <- suppressWarnings(AIC(
          #     MLPE.lmm2(
          #       resistance = cd,
          #       response = gdist.inputs$response,
          #       ID = gdist.inputs$ID,
          #       ZZ = gdist.inputs$ZZ,
          #       REML = FALSE
          #     )
          #   ))
          #   obj.func.opt <- obj.func * -1
          # } else if (method == "R2") {
          #   obj.func <- suppressWarnings(r.squaredGLMM(
          #     MLPE.lmm2(
          #       resistance = cd,
          #       response =
          #         gdist.inputs$response,
          #       ID = gdist.inputs$ID,
          #       ZZ = gdist.inputs$ZZ,
          #       REML = FALSE
          #     )
          #   ))
          #   obj.func.opt <- obj.func[[1]]
          # } else {
          #   obj.func <- suppressWarnings(logLik(
          #     MLPE.lmm2(
          #       resistance = cd,
          #       response = gdist.inputs$response,
          #       ID = gdist.inputs$ID,
          #       ZZ = gdist.inputs$ZZ,
          #       REML = FALSE
          #     )
          #   ))
          #   obj.func.opt <- obj.func[[1]]
          # }
          rm(cd, r)
          gc()
          # } # Keep loop
        } # End gdistance Loop
        
        # * CS jl ------------------------------------------------------------
        
        if (!is.null(jl.inputs)) {
          
          ## Turned off 3/15/2019
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
            
            # dat <- data.frame(gd = jl.inputs$response,
            #                   cd = scale(cd),
            #                   pop = jl.inputs$ID$pop1)
            
            dat <- jl.inputs$df
            dat$cd <- scale(cd)
            
            fit.mod <-  mlpe_rga(formula = gd ~ cd + (1 | pop),
                                 data = dat,
                                 ZZ = jl.inputs$ZZ,
                                 REML = FALSE)
            
            if(lme4::fixef(fit.mod)['cd'] < 0) {
              obj.func.opt <- -99999
            } else {
              if (method == "AIC") {
                obj.func <- suppressWarnings(AIC(
                  fit.mod
                ))
                obj.func.opt <- obj.func * -1
              } else if (method == "R2") {
                obj.func <-
                  suppressWarnings(r.squaredGLMM(
                    fit.mod
                  ))
                obj.func.opt <- obj.func[[1]]
                
              } else {
                obj.func <- suppressWarnings(logLik(
                  fit.mod
                ))
                obj.func.opt <- obj.func[[1]]
              }
            } # Positive parameter value
            
            # if (method == "AIC") {
            #   obj.func <- suppressWarnings(AIC(
            #     MLPE.lmm2(
            #       resistance = cd,
            #       response = jl.inputs$response,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ,
            #       REML = FALSE
            #     )
            #   ))
            #   obj.func.opt <- obj.func * -1
            # } else if (method == "R2") {
            #   obj.func <- suppressWarnings(r.squaredGLMM(
            #     MLPE.lmm2(
            #       resistance = cd,
            #       response =
            #         jl.inputs$response,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ,
            #       REML = FALSE
            #     )
            #   ))
            #   obj.func.opt <- obj.func[[1]]
            # } else {
            #   obj.func <- suppressWarnings(logLik(
            #     MLPE.lmm2(
            #       resistance = cd,
            #       response = jl.inputs$response,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ,
            #       REML = FALSE
            #     )
            #   ))
            #   obj.func.opt <- obj.func[[1]]
            # }
          } # End Keep loop
        } # End Julia Loop
        
      } # End drop Loop
      
      rt <- Sys.time() - t1
      if (quiet == FALSE) {
        cat(paste0("\t", "Iteration took ", round(rt, digits = 2), " seconds"), "\n")
        #     cat(paste0("\t", EQ,"; ",round(SHAPE,digits=2),"; ", round(Max.SCALE,digits=2)),"\n")
        cat(paste0("\t", method, " = ", round(obj.func.opt, 4)), "\n", "\n")
        if (!is.null(iter)) {
          if (GA.inputs$surface.type[iter] != "cat") {
            cat(paste0("\t", EQ, " | Shape = ", PARM[2], " | Max = ", PARM[3]),
                "\n",
                "\n")
          }
        }
      } 
    } else { # End transformation loop
      
      ## Use -99999 as 'ga' maximizes
      obj.func.opt <- -99999
      
      rt <- Sys.time() - t1
      if (quiet == FALSE) {
        cat(paste0("\t", "Iteration took ", round(rt, digits = 2), " seconds"), "\n")
        cat(paste0("\t", method, " = ", obj.func.opt, "\n"))
        if (!is.null(iter)) {
          if (GA.inputs$surface.type[iter] != "cat") {
            cat(paste0("EXCLUDED TRANSFORMATION", "\n", "\n"))
          }
        }
      }
    }
    # rm(r)
    gc()
    if(!is.null(GA.inputs$opt.digits)) {
      obj.func.opt <- round(obj.func.opt, GA.inputs$opt.digits)
      return(obj.func.opt)
    } else {
      return(obj.func.opt)
    }
  }