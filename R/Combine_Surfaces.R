#' Combine multiple resistance surfaces together
#'
#' Combine multiple resistance surfaces into new composite surface based on specified parameters
#'
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Requires a vector with parameters specified in the order of resistance surfaces
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param jl.inputs Object created from running \code{\link[ResistanceGA]{jl.prep}} function. Defined if optimizing using CIRCUITSCAPE run in Julia
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function.
#' @param out Directory to write combined .asc file. Default = NULL and no files are exported
#' @param File.name Name of output .asc file. Default is the combination of all surfaces combined, separated by "."
#' @param rescale Locical. If TRUE (default), the values of the combined raster surface will be divided by the minimum value to create a resistance surface with a minimum value = 1.
#' @param p.contribution Logical. If TRUE, the function will return a list object containing (1) the combined raster surface and (2) the average contribution of each surface to the resistance values of the combined surface.
#' @details \code{PARM} is designed to accept the output of \code{MS_optim}. For continuous surfaces, there are three terms: 1) Transformation, 2) shape, and 3) maximum value. Transformation must be provided as a numeric value:\cr
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
#'
#' The Distance transformation sets all values equal to one. Because of the flexibility of the Ricker function to take a monomolecular shape (try \code{Plot.trans(PARM=c(10,100), Resistance=c(1,10), transformation="Ricker")} to see this), whenever a shape parameter >6 is selected in combination with a Ricker family transformation, the transformation reverts to a Distance transformation. In general, it seems that using a combination of intermediate Ricker and Monomolecular transformations provides the best, most flexible coverage of parameter space.
#' @return R raster object that is the sum all transformed and/or reclassified resistance surfaces provided
#' @usage Combine_Surfaces(PARM, 
#'                         CS.inputs = NULL, 
#'                         gdist.inputs = NULL, 
#'                         jl.inputs = NULL,
#'                         GA.inputs, 
#'                         out = NULL, 
#'                         File.name, 
#'                         rescale = TRUE, 
#'                         p.contribution = FALSE)
#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#' 
#' @examples  
#' ## Not run:
#' ## *** TO BE COMPLETED *** ##
#' 
#' ## End (Not run)

Combine_Surfaces <-
  function(PARM,
           CS.inputs = NULL,
           gdist.inputs = NULL,
           jl.inputs = NULL,
           GA.inputs,
           out = NULL,
           File.name = paste(GA.inputs$parm.type$name, collapse = "."),
           rescale = TRUE,
           p.contribution = FALSE) {
    
    if (!is.null(CS.inputs)) {
      ID <- CS.inputs$ID
      ZZ <- CS.inputs$ZZ
      response <- CS.inputs$response
      CS_Point.File <- CS.inputs$CS_Point.File
      CS.program <- CS.inputs$CS.program
      EXPORT.dir <- out
    }
    
    if (!is.null(gdist.inputs)) {
      ID <- gdist.inputs$ID
      ZZ <- gdist.inputs$ZZ
      response <- gdist.inputs$response
      samples <- gdist.inputs$samples
      EXPORT.dir <- out
    }
    
    if (!is.null(jl.inputs)) {
      ID <- jl.inputs$ID
      ZZ <- jl.inputs$ZZ
      response <- jl.inputs$response
      samples <- jl.inputs$samples
      EXPORT.dir <- out
    }
    
    
    ######
    select.trans <- GA.inputs$select.trans
    r <- GA.inputs$Resistance.stack
    keep <- 1
    
    
    for (i in 1:GA.inputs$n.layers) {
      if (GA.inputs$surface.type[i] == "cat") {
        parm <-
          PARM[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        
        ## Prevent NaN in parm
        if(is.na(sum(parm))) {
          parm <- replace(parm, values = rep(1,length(parm)))
          keep <- 0
        }
        
        parm <- parm / min(parm)
        if(max(parm) > GA.inputs$max.cat){
          parm <- SCALE.vector(parm, 1, GA.inputs$max.cat)
        }
        
        # parm <- parm / min(parm)
        # parm <- SCALE.vector(parm, 1, GA.inputs$max.cat)
        
        df <-
          data.frame(id = unique(r[[i]]), parm) # Data frame with original raster values and replacement values
        r[[i]] <- subs(r[[i]], df)
        
        # r[[i]] <- r[[i]]#-1 # Set minimum to 0
        
        
      } else {
        rast <- SCALE(data = r[[i]],
                      MIN = 0,
                      MAX = 10)
        parm <-
          PARM[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        
        # Prevent NAs and errors
        if(is.na(parm[1])) {
          parm[1] <- 9
          keep <- 0
        }
        
        if(is.na(parm[2])) {
          parm[2] <- 1
          parm[1] <- 9
          keep <- 0
        }
        
        if(is.na(parm[3])) {
          parm[3] <- 2
          parm[1] <- 9
          keep <- 0
        }
        
        # Set equation for continuous surface
        equation <- floor(parm[1]) # Parameter can range from 1-9.99
        
        SHAPE <-  (parm[2])
        Max.SCALE <- (parm[3])
        
        # Apply specified transformation
        rick.eq <- (equation == 2 ||
                      equation == 4 || equation == 6 || equation == 8)
        
        
        if (rick.eq == TRUE & SHAPE > 5) {
          equation <- 9
          keep <- 0
        }
        
        if (equation %in% select.trans[[i]] & keep == 1) {
          equation <- equation
          keep <- 1
        } else {
          equation <- 9
          keep <- 0
        }
        
        # Apply specified transformation
        if (equation == 1) {
          r[[i]] <- Inv.Rev.Monomolecular(rast, parm)
          EQ <- "Inverse-Reverse Monomolecular"
          
        } else if (equation == 5) {
          r[[i]] <- Rev.Monomolecular(rast, parm)
          EQ <- "Reverse Monomolecular"
          
        } else if (equation == 3) {
          r[[i]] <- Monomolecular(rast, parm)
          EQ <- "Monomolecular"
          
        } else if (equation == 7) {
          r[[i]] <- Inv.Monomolecular(rast, parm)
          EQ <- "Inverse Monomolecular"
          
        } else if (equation == 8) {
          r[[i]] <- Inv.Ricker(rast, parm)
          EQ <- "Inverse Ricker"
          
        } else if (equation == 4) {
          r[[i]] <- Ricker(rast, parm)
          EQ <- "Ricker"
          
        } else if (equation == 6) {
          r[[i]] <- Rev.Ricker(rast, parm)
          EQ <- "Reverse Ricker"
          
        } else if (equation == 2) {
          r[[i]] <- Inv.Rev.Ricker(rast, parm)
          EQ <- "Inverse-Reverse Ricker"
          
        } else {
          r[[i]] <- (rast * 0) + 1 #  Cancel layer...set to zero
        } # End if-else
      } # Close parameter type if-else
    } # Close layer loop
    
    
    # File.name <- File.name
    
    # ms.r <- sum(r) # Add all surfaces together
    multi_surface <- sum(r) # Add all surfaces together
    
    # If unused transformation applied, toss iteration
    if(keep == 0) {
      # ms.r <- (r[[1]] * 0)
      multi_surface <- (r[[1]] * 0)
    }
    
    
    
    if (is.null(out)) {
      if(p.contribution == TRUE) {
        cont.list <- vector(mode = 'list', length = GA.inputs$n.layers)
        
        for (i in 1:GA.inputs$n.layers) {
          p.cont <- r[[i]] / multi_surface
          # mean.cont <- cellStats(p.cont, mean, na.rm = TRUE)
          mean.cont <- mean(p.cont@data@values, na.rm = TRUE)
          cont.list[[i]] <- data.frame(surface = GA.inputs$layer.names[i], 
                                       mean = mean.cont)
        }
        
        if (keep != 0 && rescale == TRUE)
          multi_surface <-
            multi_surface / min(multi_surface@data@values, na.rm = TRUE) # Rescale to min of 1
        
        if (keep != 0 && max(multi_surface@data@values, na.rm = TRUE) > 1e6)
          multi_surface <-
            SCALE(multi_surface, 1, 1e6) # Rescale surface in case resistance are too high
        
        cont.df <- plyr::ldply(cont.list)
        
        # Work around for NA raster surfaces
        ## Memory issues?
        # if(is.na(cellStats(multi_surface, mean))) {
        #   multi_surface <- (GA.inputs$Resistance.stack[[1]] * 0)
        # }
        
        if(sum(!is.na(multi_surface@data@values)) == 0) {
          keep <- 0
          multi_surface <- (GA.inputs$Resistance.stack[[1]] * 0)
        }
        
        list.out <- list(percent.contribution = cont.df,
                         combined.surface = multi_surface)
        return(list.out)
        
      } else { # p.contribution = FALSE
        
        if (keep != 0 && rescale == TRUE)
          multi_surface <-
            multi_surface / min(multi_surface@data@values, na.rm = TRUE) # Rescale to min of 1
        
        if (keep != 0 && max(multi_surface@data@values, na.rm = TRUE) > 1e6)
          multi_surface <-
            SCALE(multi_surface, 1, 1e6) # Rescale surface in case resistance are too high
        
        rm(r)
        gc()
        return(multi_surface)
        
      }
      
    } else {   # Output directory specified
      
      if (keep != 0 && rescale == TRUE)
        multi_surface <-
          multi_surface / min(multi_surface@data@values, na.rm = TRUE) # Rescale to min of 1
      
      if (keep != 0 && max(multi_surface@data@values, na.rm = TRUE) > 1e6)
        multi_surface <-
          SCALE(multi_surface, 1, 1e6) # Rescale surface in case resistance are too high
      
      
      # Work around for NA raster surfaces
      ## Memory issues?
      # if(is.na(cellStats(multi_surface, mean))) {
      #   multi_surface <- (GA.inputs$Resistance.stack[[1]] * 0)
      # }
      
      if(sum(!is.na(multi_surface@data@values)) == 0) {
        keep <- 0
        multi_surface <- (GA.inputs$Resistance.stack[[1]] * 0)
      }
      
      writeRaster(
        x = multi_surface,
        filename = paste0(EXPORT.dir, File.name, ".asc"),
        overwrite = TRUE
      )
      
      rm(r)
      gc()
      return(multi_surface)
    }
  }
