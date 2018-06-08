## FOR ASSESSING AICc of fitted models, not exported
Resistance.Opt_AICc <-
  function(PARM,
           Resistance,
           CS.inputs = NULL,
           gdist.inputs = NULL,
           GA.inputs,
           Min.Max = 'max',
           iter = NULL,
           quiet = FALSE) {
    t1 <- proc.time()[3]
    
    
    EXPORT.dir <- GA.inputs$Write.dir
    
    r <- Resistance
    if (!is.null(iter)) {
      if (GA.inputs$surface.type[iter] == "cat") {
        PARM <- PARM / min(PARM)
        parm <- PARM
        df <-
          data.frame(id = unique(r), PARM) # Data frame with original raster values and replacement values
        r <- subs(r, df)
        
      } else {
        r <- SCALE(r, 0, 10)
        
        # Set equation for continuous surface
        equation <- floor(PARM[1]) # Parameter can range from 1-9.99
        
        # Read in resistance surface to be optimized
        SHAPE <- (PARM[2])
        Max.SCALE <- (PARM[3])
        
        # Apply specified transformation
        rick.eq <- (equation == 2 ||
                      equation == 4 ||
                      equation == 6 || equation == 8)
        if (rick.eq == TRUE & SHAPE > 5) {
          equation <- 9
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
      } # Close parameter type if-else
    } else {
      r <- SCALE(r, 0, 10)
      
      # Set equation for continuous surface
      equation <- floor(PARM[1]) # Parameter can range from 1-9.99
      
      # Read in resistance surface to be optimized
      SHAPE <- (PARM[2])
      Max.SCALE <- (PARM[3])
      
      # Apply specified transformation
      rick.eq <- (equation == 2 ||
                    equation == 4 || equation == 6 || equation == 8)
      if (rick.eq == TRUE & SHAPE > 5) {
        equation <- 9
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
    }
    File.name <- "resist_surface"
    if (cellStats(r, "max") > 1e6)
      r <-
      SCALE(r, 1, 1e6) # Rescale surface in case resistance are too high
    r <- reclassify(r, c(-Inf, 1e-06, 1e-06, 1e6, Inf, 1e6))
    
    
    if (!is.null(CS.inputs)) {
      writeRaster(
        x = r,
        filename = paste0(EXPORT.dir, File.name, ".asc"),
        overwrite = TRUE
      )
      CS.resist <-
        Run_CS2(
          CS.inputs,
          GA.inputs,
          r = r,
          EXPORT.dir = GA.inputs$Write.dir,
          File.name = File.name
        )
      
      # Replace NA with 0...a workaround for errors when two points fall within the same cell.
      # CS.resist[is.na(CS.resist)] <- 0
      
      # Run mixed effect model on each Circuitscape effective resistance
      AIC.stat <- suppressWarnings(AIC(
        MLPE.lmm2(
          resistance = CS.resist,
          response = CS.inputs$response,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ,
          REML = FALSE
        )
      ))
      ROW <- nrow(CS.inputs$ID)
      
    }
    
    if (!is.null(gdist.inputs)) {
      cd <- Run_gdistance(gdist.inputs, r)
      
      AIC.stat <- suppressWarnings(AIC(
        MLPE.lmm2(
          resistance = cd,
          response = gdist.inputs$response,
          ID = gdist.inputs$ID,
          ZZ = gdist.inputs$ZZ,
          REML = FALSE
        )
      ))
      ROW <- nrow(gdist.inputs$ID)
    }
    
    k <- length(PARM) + 1
    AICc <- (AIC.stat) + (((2 * k) * (k + 1)) / (ROW - k - 1))
    
    rt <- proc.time()[3] - t1
    if (quiet == FALSE) {
      cat(paste0("\t", "Iteration took ", round(rt, digits = 2), " seconds"), "\n")
      #     cat(paste0("\t", EQ,"; ",round(SHAPE,digits=2),"; ", round(Max.SCALE,digits=2)),"\n")
      cat(paste0("\t", "AICc = ", round(AICc, 4)), "\n")
      if (!is.null(iter)) {
        if (GA.inputs$surface.type[iter] != "cat") {
          cat(paste0("\t", EQ, " | Shape = ", PARM[2], " | Max = ", PARM[3]),
              "\n",
              "\n")
        }
      }
    }
    OPTIM.DIRECTION(Min.Max) * (AICc) # Function to be minimized/maximized
  }

# FUNCTIONS
OPTIM.DIRECTION <- function(x) {
  OPTIM <- ifelse(x == 'max', -1, 1)
  return(OPTIM)
}


Cont.Param <- function(PARM) {
  df <- data.frame(PARM[1], PARM[2])
  colnames(df) <- c("shape_opt", "max")
  row.names(df) <- NULL
  return(df)
}


read.matrix <-
  function(cs.matrix) {
    lower(read.table(cs.matrix)[-1, -1])
  }


read.matrix2 <- function(cs.matrix) {
  m <- read.table(cs.matrix)[-1, -1]
}


# Create ZZ matrix for mixed effects model
ZZ.mat <- function(ID) {
  Zl <-
    lapply(c("pop1", "pop2"), function(nm)
      Matrix::fac2sparse(ID[[nm]], "d", drop = FALSE))
  ZZ <- Reduce("+", Zl[-1], Zl[[1]])
  return(ZZ)
}

# Rescale function
SCALE.vector <- function(data, MIN, MAX, threshold = 1e-5) {
  if (abs(MIN - MAX) < threshold) {
    data[is.finite(data)] <- 0
    data
  } else {
    Mn = min(data)
    Mx = max(data)
    (MAX - MIN) / (Mx - Mn) * (data - Mx) + MAX
  }
}

# Define scaling function
# This will rescale from 1 to specified MAX
SCALE <- function(data, MIN, MAX, threshold = 1e-5) {
  if (abs(MIN - MAX) < threshold) {
    data[is.finite(raster::values(data))] <- 0
    data
  } else {
    Mn = cellStats(data, stat = 'min')
    Mx = cellStats(data, stat = 'max')
    (MAX - MIN) / (Mx - Mn) * (data - Mx) + MAX
  }
}

# Sample values for suggests
sv.cat <- function(levels, pop.size, min, max) {
  cat.starts <- matrix(nrow = pop.size, ncol = levels)
  for (r in 1:pop.size) {
    L <- list()
    for (i in 1:levels) {
      if (runif(1) < .5) {
        z <- runif(1)
      } else {
        z <- runif(1, min, max)
      }
      L[[i]] <- z
    }
    #   uz<-unlist(L)
    cat.starts[r, ] <- (unlist(L))
  }
  cat.starts[, 1] <- 1
  return(cat.starts)
}


# No Gaussian distribution
sv.cont.nG <- function(direction,
                       pop.size,
                       max,
                       min.scale,
                       max.scale,
                       scale = NULL, 
                       eqs) {
  inc <- c(1, 3)
  dec <- c(7, 5)
  peak <- c(2, 4, 6, 8)
  L <- list()
  
  if (!is.null(scale)) {
    cont.starts <- matrix(nrow = pop.size, ncol = 4)
    for (r in 1:pop.size) {
      scale.parm <- runif(1, min.scale, max.scale)
      if (runif(1) < .5 && direction == "Increase") {
        #       z1<-c(sample(inc,1)
        z <- Increase.starts.nG(sample(inc, 1))
        z <- c(z, scale.parm)
      } else if (runif(1) < .5 && direction == "Decrease") {
        z <- c(sample(dec, 1),
               runif(1, .01, 10),
               runif(1, 1, max),
               scale.parm)
      } else if (runif(1) < .5 && direction == "Peaked") {
        z <- c(sample(peak, 1),
               runif(1, .01, 10),
               runif(1, 1, max),
               scale.parm)
      } else {
        z <- c(sample(eqs, 1),
               # runif(1, 1, 9.99),
               runif(1, .01, 10),
               runif(1, 1, max),
               scale.parm)
      }
      cont.starts[r,] <- z
    }
  } else {
    cont.starts <- matrix(nrow = pop.size, ncol = 3)
    for (r in 1:pop.size) {
      if (runif(1) < .5 && direction == "Increase") {
        #       z1<-c(sample(inc,1)
        z <- Increase.starts.nG(sample(inc, 1))
      } else if (runif(1) < .5 && direction == "Decrease") {
        z <- c(sample(dec, 1), runif(1, .01, 10), runif(1, 1, max))
      } else if (runif(1) < .5 && direction == "Peaked") {
        z <- c(sample(peak, 1), runif(1, .01, 10), runif(1, 1, max))
      } else {
        z <- c(sample(eqs, 1),
               # runif(1, 1, 9.99), 
               runif(1, .01, 10), 
               runif(1, 1, max))
      }
      cont.starts[r,] <- z
    }
  }
  if(ncol(cont.starts) == 4) {
    rs <- sample(pop.size, floor(0.25 * pop.size), replace = F)
    cont.starts[rs, 4] <- 0.25
  }
  
  cont.starts
}

Increase.starts <- function(x) {
  if (x == 1) {
    z <- c(x, runif(1, .01, 10), runif(1, 1, 10), 1)
  } else {
    z <- c(x, runif(1, .01, 10), runif(1, 1, 100), 1)
  }
}

Increase.starts.nG <- function(x) {
  if (x == 1) {
    z <- c(x, runif(1, .01, 10), runif(1, 1, 10))
  } else {
    z <- c(x, runif(1, .01, 10), runif(1, 1, 100))
  }
}



unique <- raster::unique

eq.set <- function(include.list) {
  out <- vector(mode = "list", length = length(include.list))
  for (i in seq_along(include.list)) {
    if (include.list[[i]] == "A" | is.na(include.list[[i]])) {
      out <- 1:9
      return(out)
    } else if (include.list[[i]] == "M") {
      out <- c(1, 3, 5, 7, 9)
      return(out)
    } else if (include.list[[i]] == "R") {
      out <- c(2, 4, 6, 8, 9)
      return(out)
    } else if (!is.na(match(include.list[[i]], 1:9))) {
      out <- include.list
      return(out)
    } else {
      cat(
        "The specified transformations to assess are not valid. \n
        Please see Details of the GA.prep."
      )
    }
  }
}

get.EQ <- function(equation) {
  # Apply specified transformation
  if (is.numeric(equation)) {
    equation = floor(equation)
    if (equation == 1) {
      EQ <- "Inverse-Reverse Monomolecular"
      
    } else if (equation == 5) {
      EQ <- "Reverse Monomolecular"
      
    } else if (equation == 3) {
      EQ <- "Monomolecular"
      
    } else if (equation == 7) {
      EQ <- "Inverse Monomolecular"
      
    } else if (equation == 8) {
      EQ <- "Inverse Ricker"
      
    } else if (equation == 4) {
      EQ <- "Ricker"
      
    } else if (equation == 6) {
      EQ <- "Reverse Ricker"
      
    } else if (equation == 2) {
      EQ <- "Inverse-Reverse Ricker"
      
    } else {
      EQ <- "Distance"
    }
    
    (EQ)
  } else {
    if (equation == "Inverse-Reverse Monomolecular") {
      EQ <- 1
      
    } else if (equation == "Reverse Monomolecular") {
      EQ <- 5
      
    } else if (equation == "Monomolecular") {
      EQ <- 3
      
    } else if (equation == "Inverse Monomolecular") {
      EQ <- 7
      
    } else if (equation == "Inverse Ricker") {
      EQ <- 8
      
    } else if (equation == "Ricker") {
      EQ <- 4
      
    } else if (equation == "Reverse Ricker") {
      EQ <- 6
      
    } else if (equation == "Inverse-Reverse Ricker") {
      EQ <- 2
      
    } else {
      EQ <- 9
    }
    
    (EQ)
  }
}

Result.txt <-
  function(GA.results,
           GA.inputs,
           method,
           k,
           Run.Time,
           fit.stats,
           MLPE.coef = NULL,
           optim,
           aic,
           AICc,
           LL) {
    summary.file <-
      paste0(GA.inputs$Results.dir, "Multisurface_Optim_Summary.txt")
    # AICc<-GA.results@fitnessValue
    # AICc<-round(AICc,digits=4)
    ELITE <- floor(GA.inputs$percent.elite * GA.inputs$pop.size)
    #   mlpe.results<-MLPE.lmm_coef(GA.inputs$Results.dir,genetic.dist=CS.inputs$response)
    
    sink(summary.file)
    cat(paste0(
      "Summary from multisurface optimization run conducted on ",
      Sys.Date()
    ),
    "\n")
    cat("\n")
    cat(paste0("Optimized using: ", method), "\n")
    cat("\n")
    cat(paste0("Objective function: ", optim), "\n")
    cat("\n")
    cat(paste0("Surfaces included in optimization:"), "\n")
    cat(GA.inputs$parm.type$name, "\n")
    cat("\n")
    cat("Genetic Algorithm optimization settings:")
    cat("\n")
    cat(paste0("Type of genetic algorithm used: ", GA.results@type),
        "\n")
    cat(paste0("popSize at each iteration: ", GA.results@popSize),
        "\n")
    cat(paste0("Maximum number of iterations: ", GA.results@maxiter),
        "\n")
    cat(paste0(
      "Number individuals retained each generation (elistism): ",
      ELITE
    ),
    "\n")
    cat(paste0("Crossover probability: ", GA.results@pcrossover),
        "\n")
    cat(paste0("Mutation probability: ", GA.results@pmutation), "\n")
    cat("\n")
    cat(
      paste0(
        "The Genetic Algorithm completed after ",
        GA.results@iter,
        " iterations"
      ),
      "\n"
    )
    cat("\n")
    cat(paste0("k =  ", k), "\n")
    cat("\n")
    cat(paste0("Minimum AIC: ", aic), "\n")
    cat("\n")
    cat(paste0("AICc: ", AICc), "\n")
    cat("\n")
    cat(paste0("Pseudo marginal R-square (R2m): ", fit.stats[[1]]), "\n")
    cat(paste0("Pseudo conditional R-square (R2c): ", fit.stats[[2]]),
        "\n")
    cat("\n")
    cat(paste0("Log Likelihood: ", LL), "\n")
    cat("\n")
    cat(paste0("Optimized values for each surface:"), "\n")
    cat(GA.results@solution, "\n")
    cat("\n")
    cat(paste0("Optimization took ", Run.Time, " seconds to complete"),
        "\n")
    sink()
  }

##########################################################################################

get.name <- function(x) {
  nm <- deparse(substitute(x))
  return(nm)
}


# Transformation Eqs ------------------------------------------------------

Monomolecular <- function(r, parm) {
  parm[3] * (1 - exp(-1 * r / parm[2])) + 1 # Monomolecular
}

Inv.Monomolecular <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    R <- parm[3] * (exp(-1 * r / parm[2]))
    (R <- (R - cellStats(R, stat = "min")) + 1)
  } else {
    R <- parm[3] * (exp(-1 * r / parm[2]))
    (R <- (R - min(R)) + 1)
  }
}

Inv.Rev.Monomolecular <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Inv.Monomolecular(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Inv.Monomolecular(rev.rast, parm)
  }
}

Rev.Monomolecular <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Monomolecular(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Monomolecular(rev.rast, parm)
  }
}


Ricker <- function(r, parm) {
  parm[3] * r * exp(-1 * r / parm[2]) + 1 # Ricker
}

Inv.Ricker <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    R <- (-1 * parm[3]) * r * exp(-1 * r / parm[2]) - 1 # Ricker
    R <-
      SCALE(R,
            MIN = abs(cellStats(R, stat = 'max')),
            MAX = abs(cellStats(R, stat = 'min'))) # Rescale
  } else {
    R <- (-1 * parm[3]) * r * exp(-1 * r / parm[2]) - 1 # Ricker
    R <- SCALE.vector(R, MIN = abs(max(R)), MAX = abs(min(R))) # Rescale
  }
}

Inv.Rev.Ricker <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Inv.Ricker(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Inv.Ricker(rev.rast, parm)
  }
}

Rev.Ricker <- function(r, parm) {
  if (class(r) == "RasterLayer") {
    rev.rast <- SCALE((-1 * r), 0, 10)
    Ricker(rev.rast, parm)
  } else {
    rev.rast <- SCALE.vector((-1 * r), 0, 10)
    Ricker(rev.rast, parm)
  }
}

yn.question <- function(question, add_lines_before = TRUE) {
  choices <- c("Yes", "No", "New Subdirectory")
  if(add_lines_before) cat("------------------------\n")   
  the_answer <- menu(choices, title = question)
  
  if(the_answer == 1L) {
    return(TRUE)
  } else if(the_answer == 2L) {
    return(FALSE)
  } else if(the_answer == 3L){
    return(NA)
  }
  
  # ifelse(the_answer == 1L, TRUE, FALSE)   # returns TRUE or FALSE
}