#' Prepare and bundle input CIRCUITSCAPE model parameters to run in Julia
#'
#' This function will prepare objects needed for running optimization functions
#'
#' @param n.Pops The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix). Not necessary if only executing Julia run. If optimizing with select points only (see \code{pairs_to_include} below), still provide the full vector for the lower half of the square matrix. This will automatically be trimmed down based on the selected pairs to include.
#' @param CS_Point.File Provide a \code{\link[sp]{SpatialPoints}} object containing sample locations. Alternatively, specify the path to the Circuitscape formatted point file. See Circuitscape documentation for help.
#' @param covariates Data frame of additional covariates that you want included in the MLPE model during opitmization.
#' @param formula If covariates are included in the model, specify the R formula for the fixed effects portion of the MLPE model.
#' @param JULIA_HOME Path to the folder containing the Julia binary (See Details)
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
#' @param pairs_to_include Default is NULL. If you wish to use the advanced CIRCUITSCAPE setting mode to include or exclude certain pairs of sample locations, provide a vector consisting of 1 (keep) or 0 (drop) for each pairwise observation (see example). This is an option if you do not want to assess all pairiwse observations in the MLPE model. 
#' @param pop2ind Default is NULL. A vector with length equivalent to the number of populations sapled. Each value of the vector indicates the number of samples from each respective population.
#' @param nb Default is NULL. Optionally, specify a distance between points that should be considered as part of the same group or neighborhood. See Details.
#' @param parallel (Logical; Default = FALSE) Do you want to run CIRCUITSCAPE in parallel?
#' @param cores If `parallel = TRUE`, how many cores should be used for parallel processing?
#' @param cholmod (Logical; Default = TRUE). Should the cholmod solver be used? See details.
#' @param precision (Logical; Default = FALSE). Should experimental single precision method be used? See details.
#' @param run_test (Logical; Default = TRUE). Should a test of Julia Circuitscape be conducted? (This can take several seconds to complete)
#' @param write.files (Default = NULL). If a directory is specified, then the .ini and .asc files used in the CS.jl run will be exported.
#' @param write.criteria Criteria for writing .ini and .asc files. If a time in seconds is not specified, then all files will be written if a \code{write.files} directory is specified.
#' @param silent (Default = TRUE) No updates or logging of CIRCUITSCAPE will occur. May be useful to set to FALSE to debug. 
#' @param Julia_link Specify whether R should connect to Julia using the 'JuliaCall' package (Default), or the 'XRJulia' package
#' @param scratch Scratch directory for use if write access is limited. 
#' @param rm.files Should all temporary files be removed after Julia run (Default = TRUE)

#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Peterman.73@@osu.edu>
#' @usage jl.prep(n.Pops, 
#'                response = NULL, 
#'                CS_Point.File, 
#'                covariates = NULL,
#'                formula = NULL,
#'                JULIA_HOME = NULL,
#'                Neighbor.Connect = 8, 
#'                pairs_to_include = NULL, 
#'                pop2ind = NULL,
#'                nb = NULL,
#'                parallel = FALSE, 
#'                cores = NULL,
#'                cholmod = TRUE,
#'                precision = FALSE, 
#'                run_test = TRUE,
#'                write.files = NULL,
#'                write.criteria = NULL,
#'                silent = TRUE,
#'                Julia_link = 'JuliaCall',
#'                scratch = NULL,
#'                rm.files = TRUE)
#' @details 
#' This function requires that Julia is properly installed on your system. Upon first running of this function, the Circuitscape.jl library will be downloaded and tested. (see https://github.com/Circuitscape/Circuitscape.jl for more details). This may take some time.
#' 
#' Using \code{cholmod} (see https://github.com/Circuitscape/Circuitscape.jl)
#' "The cholesky decomposition is a direct solver method, unlike the algebraic multigrid method used by default in both the old and the new version. The advantage with this new direct method is that it can be much faster than the iterative solution, within a particular problem size. 
#' Word of caution: The cholesky decomposition is not practical to use beyond a certain problem size because of phenomenon called fill-in, which results in loss of sparsity and large memory consumption."
#' The cholmod solver can only be used when \code{precision} `= FALSE` (double precision).
#' 
#' If \code{precision} is TRUE, then the EXPERIMENTAL single precision method will be used. Single precision usually uses less memory, but is likely to reduce accuracy. NOTE: Preliminary testing of single precision mode in a Windows pc resulted in extremely slow runs.
#' 
#' \code{JULIA_HOME} is where the Julia binary files are stored. Usually in a `bin` directory within the Julia install directory.
#' 
#' When specifying a formula, provide it as: \code{response ~ covariate}.
#' the formula \code{response} will use the vector of values specified for the \code{response} parameter. Make sure that covariate names match variable names provided in \code{covariates}
#' 
#' If \code{nb} is specified, this indicates the maximum distance (in units of your \code{CS_Point.File} SpatialPoints object) that locations are considered to be part of the same neighborhood.
#' 
#' 
#' @examples  
#' ## Not run:
#' ## *** TO BE COMPLETED *** ##
#' 
#' ## End (Not run)
jl.prep <- function(n.Pops,
                    response = NULL,
                    CS_Point.File,
                    covariates = NULL,
                    formula = NULL,
                    JULIA_HOME = NULL,
                    Neighbor.Connect = 8,
                    pairs_to_include = NULL,
                    pop2ind = NULL,
                    nb = NULL,
                    parallel = FALSE,
                    cores = NULL,
                    cholmod = TRUE,
                    precision = FALSE,
                    run_test = TRUE,
                    write.files = NULL,
                    write.criteria = NULL,
                    silent = TRUE,
                    Julia_link = 'JuliaCall',
                    scratch = NULL,
                    rm.files = TRUE) {
  
  
  # Checks ------------------------------------------------------------------
  
  # Ensure covariates have same length as response
  if(!is.null(covariates) && !is.data.frame(covariates)) {
    stop("Please provide a data frame when specifying additional covariates")
  } 
  
  if(!is.null(covariates) && nrow(covariates) != length(response)) {
    stop("Response and covariates must have the same number of observations")
  } 
  
  if(Julia_link == 'XRJulia') {
    JULIA_HOME <- findJulia() 
  }
  
  if(!is.null(write.files)) {
    write.files <- paste0(normalizePath(write.files, winslash = '/'), "/")
    if(!dir.exists(write.files))
      stop("`write.files` directory does not exist")
  }
  
  # if(!is.null(pop2ind)) {
  #   if(!is.vector(pop2ind)) {
  #     stop("`pop2ind` must be a vector of population sizes")
  #   }
  #   keep.df <- expand.mat_vec(pop2ind, 
  #                             gd = response)
  # }
  
  # Setup Julia -------------------------------------------------------------
  
  if(Julia_link != 'XRJulia'){
    if(!is.null(JULIA_HOME)) {
      if(isTRUE(dir.exists(JULIA_HOME))) {
        julia_setup(JULIA_HOME = JULIA_HOME)
      } else {
        stop("Specified JULIA_HOME directory does not exist")
      } 
      
    } else {
      jl.setup <- try(julia_setup(), TRUE)
      JULIA_HOME <- XRJulia::findJulia()
      JULIA_HOME <- paste0(dirname(JULIA_HOME), "/")
      
      if(class(jl.setup) == "try-error")
        stop("Specified JULIA_HOME directory does not exist")
    }        
  }
  JULIA_HOME <- normalizePath(JULIA_HOME)
  
  
  # Determine if CIRCUITSCAPE package is installed
  # jl.cs <- try(julia_library("Circuitscape"), TRUE)
  
  # if(class(jl.cs) == "try-error"){
  #   stop(cat(paste("You must install the Julia CIRCUITSCAPE package!!!",
  #                  "https://github.com/Circuitscape/Circuitscape.jl", sep = "\n")))
  # }
  
  jl.cs <- julia_installed_package("Circuitscape")
  
  if(jl.cs == "nothing"){
    stop(cat(paste("You must install the Julia CIRCUITSCAPE package!!!",
                   "https://github.com/Circuitscape/Circuitscape.jl", sep = "\n")))
  }
  
  try(JuliaCall::julia_library("Circuitscape"), silent = T)
  JuliaCall::julia_library("Circuitscape")
  
  wd <- getwd()
  
  if(run_test == TRUE) {
    print("Test: Run Circuitscape from Julia")
    
    if(Sys.info()[['sysname']] == "Windows") {
      td <- paste0(normalizePath(tempdir(), winslash = "/"), "/")
      # setwd(JULIA_HOME)
      
    } else {
      td <- paste0(tempdir(),"/")
    }
    
    if(!is.null(scratch)) {
      if(Sys.info()[['sysname']] == "Linux") {
        scratch <- td <- paste0(normalizePath(scratch),"/")
      } else {
        scratch <- td <- normalizePath(scratch, winslash = "/")
      }
    }
    
    # if(!is.null(scratch)) {
    #   if(Sys.info()[['sysname']] == "Windows") {
    #     scratch <- td <- normalizePath(scratch)
    #   } else {
    #     scratch <- td <- scratch
    #   }
    # }
    
    write.table(ResistanceGA::samples[1:5,], 
                paste0(td,'samples.txt'), 
                quote = FALSE,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE)
    
    
    temp.ini <- tempfile(pattern = "", 
                         tmpdir = tempdir(),
                         fileext = ".ini")
    
    if(!is.null(scratch)) {
      temp.ini <- paste0(scratch, basename(temp.ini))
    }
    
    
    tmp.name <- basename(temp.ini) %>% strsplit(., '.ini') %>% unlist()
    
    # if(!is.null(scratch)) {
    #   writeRaster(resistance_surfaces$continuous,
    #               paste0(scratch,
    #                      tmp.name, '.asc'),
    #               overwirte = TRUE)
    #   
    #   write.CS_4.0(BATCH = paste0(scratch, tmp.name, ".ini"),
    #                OUT = paste0("output_file = ", scratch,"//", tmp.name, ".out"),
    #                HABITAT = paste0("habitat_file = ", scratch, tmp.name, '.asc'),
    #                LOCATION.FILE = paste0("point_file = ", scratch, 'samples.txt'),
    #                PARALLELIZE = FALSE,
    #                CORES = NULL,
    #                solver = 'cholmod',
    #                precision = FALSE,
    #                silent = silent
    #   )
    # } else {
    writeRaster(resistance_surfaces$continuous,
                paste0(td,
                       tmp.name, '.asc'),
                overwirte = TRUE)
    
    write.CS_4.0(BATCH = paste0(td, tmp.name, ".ini"),
                 OUT = paste0("output_file = ", td, tmp.name, ".out"),
                 HABITAT = paste0("habitat_file = ", td, tmp.name, '.asc'),
                 LOCATION.FILE = paste0("point_file = ", td, 'samples.txt'),
                 PARALLELIZE = FALSE,
                 CORES = NULL,
                 solver = 'cholmod',
                 precision = FALSE,
                 silent = silent
    )
    # }
    
    if(Julia_link == 'JuliaCall'){
      out <- julia_call('compute', temp.ini)[-1,-1]
      
    } else {
      cs.jl <- RJulia()
      cs.jl$Using("Circuitscape")
      
      cs.out <- cs.jl$Call("compute", temp.ini) 
      out <- as.matrix(read.table(paste0(scratch, "//", tmp.name, "_resistances.out"),
                                  quote="\"", comment.char=""))[-1,-1]
      # out <- read.delim(paste0(scratch, "/", tmp.name, "_resistances.out"), header = FALSE)[-1,-1]
      # out <- juliaGet(cs.out)[-1,-1] ## SLOW!!!
    }
    
    if(wd != getwd()) {
      setwd(wd)
    }
    
    
    if(dim(out)[1] == 5) {
      cat("\n"); cat("\n")
      
      cat("Test Passed")
      
      cat("\n"); cat("\n")
      
    } else {
      stop("Test Failed")
    }
    if(isTRUE(rm.files)) {
      unlink.list <- list.files(td,
                                pattern = tmp.name,
                                all.files = TRUE,
                                full.names = TRUE)
      del.files <- sapply(unlink.list, unlink, force = TRUE)
    }
    
    if(!is.null(scratch)){
      unlink.list2 <- list.files(scratch,
                                 pattern = tmp.name,
                                 all.files = TRUE,
                                 full.names = TRUE)
      del.files <- sapply(unlink.list2, unlink, force = TRUE)
    }
  } # End test
  
  
  # Format inputs -----------------------------------------------------------
  if(isTRUE(precision)) {
    precision <- 'single'
  } else {
    precision <- 'None'
  }
  
  if(isTRUE(cholmod) && (precision == 'single')) {
    stop(cat(paste0('\n', 
                    "jl.prep ERROR:", '\n',
                    "CHOLMOD solver only works when using double precision. Set either `cholmod = FALSE` OR `precision = FALSE` to proceed", 
                    '\n', '\n' ))
    )
  }
  
  if(isTRUE(cholmod)) {
    solver <- 'cholmod'
  } else {
    solver <- NULL
  }
  
  if(class(CS_Point.File) == 'SpatialPoints') {
    
    td <- paste0(tempdir(),"//")
    
    if(!is.null(scratch)) {
      if(Sys.info()[['sysname']] == "Windows") {
        td <- normalizePath(scratch)
      } else {
        td <- scratch
      }
    }
    
    site <- c(1:length(CS_Point.File))
    
    cs.txt <- data.frame(site, CS_Point.File)
    cs_coords <- cs.txt[,2:3]
    write.table(
      cs.txt,
      # file = sp_file,
      file = paste0(td, "sample_pts.txt"),
      col.names = F,
      row.names = F
    )
    
    # CS_Point.File <- sp_file
    CS_Point.File <- paste0(td, "sample_pts.txt")
  }
  
  if (grepl(".asc", x = CS_Point.File)) {
    CS_grid <- raster <- raster(CS_Point.File)
    CS_Point.txt <- rasterToPoints(x = CS_grid)
    site <- CS_Point.txt[, 3]
    CS_Point.txt <- data.frame(site, CS_Point.txt[, c(1, 2)])
    CS_Point.txt <- CS_Point.txt[order(site), ]
    cs_coords <- CS_Point.txt[, c(1, 2)]
    CS_Point.File <- sub(".asc", ".txt", x = CS_Point.File)
    write.table(
      CS_Point.txt,
      file = CS_Point.File,
      col.names = F,
      row.names = F
    )
  }
  
  if (!is.null(response)) {
    TEST.response <- ((is.vector(response)) || (ncol(response) == 1))
    if (TEST.response == FALSE) {
      stop("The object 'response' is not in the form of a single column vector")
    }
  }
  
  
  # Make data frame ---------------------------------------------------------
  # Make to-from population list
  if (!exists(x = "ID")) {
    if(!is.null(nb)){
      ID <- To.From.ID(sampled_pops = n.Pops,
                       pop_n = pop2ind,
                       spLoc = CS_Point.File,
                       nb = nb)
    } else {
      ID <- To.From.ID(sampled_pops = n.Pops,
                       pop_n = pop2ind,
                       spLoc = NULL,
                       nb = nb)
    }
  }
  
  ## * Update formula?
  suppressWarnings(ZZ <- ZZ.mat(ID))
  
  
  fmla <- formula
  df <- NULL
  
  if(!is.null(response)){
    if(!is.null(covariates)){
      if(!is.null(pop2ind)){
        keep.vec <- expand.keep(pop2ind)
        gd.mat <- matrix(0, n.Pops, n.Pops)
        
        
        cov.list <- vector('list', length(names(covariates)))
        for(i in 1:length(cov.list)){
          mat <- gd.mat
          mat[lower.tri(mat)] <- covariates[,i]
          
          cov.list[[i]] <- expand.mat(mat,
                                      pop2ind)
        }
        
        names(cov.list) <- names(covariates)
        cov.df <- as.data.frame(cov.list)
        covariates <- cov.df
        
        response <- response[keep.vec == 1]
        
        if(!is.null(nb)){
          
          df <- data.frame(gd = response,
                           cov.df,
                           pop = ID$pop,
                           grp = ID$pop1.pop,
                           cor.grp = ID$cor.grp)
        } else {
          df <- data.frame(gd = response,
                           cov.df,
                           pop = ID$pop,
                           grp = ID$pop1.pop)
        }
      } else {
        if(!is.null(nb)){
          df <- data.frame(gd = response,
                           covariates,
                           pop = ID$pop1,
                           cor.grp = ID$cor.grp)
        } else {
          df <- data.frame(gd = response,
                           covariates,
                           pop = ID$pop1)
        }
      }
      
    } else {
      if(!is.null(pop2ind)){
        keep.vec <- expand.keep(pop2ind)
        gd.mat <- matrix(0, n.Pops, n.Pops)
        response <- response[keep.vec == 1]
        
        if(!is.null(nb)){
          
          df <- data.frame(gd = response,
                           pop = ID$pop,
                           grp = ID$pop1.pop,
                           cor.grp = ID$cor.grp)
        } else {
          df <- data.frame(gd = response,
                           pop = ID$pop,
                           grp = ID$pop1.pop)
        } 
      } else {
        if(!is.null(nb)){
          df <- data.frame(gd = response,
                           pop = ID$pop1,
                           cor.grp = ID$cor.grp)
        } else {
          df <- data.frame(gd = response,
                           pop = ID$pop1)
        }
      }
    }
  }
  
  if(!is.null(fmla)) {
    if(!is.null(pop2ind)){
      if(!is.null(nb)){
        fmla <- update(fmla, gd ~ . + cd + (1 | pop) + (1 | grp) + (1 | cor.grp))
        
      } else {
        fmla <- update(fmla, gd ~ . + cd + (1 | pop) + (1 | grp))
        
      }
      
    } else {
      if(!is.null(nb)){
        fmla <- update(fmla, gd ~ . + cd + (1 | pop) + (1 | cor.grp))
      } else {
        fmla <- update(fmla, gd ~ . + cd + (1 | pop))
      }
    }
  } else {
    if(!is.null(pop2ind)){
      if(!is.null(nb)){
        fmla <- gd ~ cd + (1 | pop) + (1 | grp) + (1 | cor.grp)
        
      } else {
        fmla <- gd ~ cd + (1 | pop) + (1 | grp)
        
      }
      
    } else {
      if(!is.null(nb)){
        fmla <- gd ~ cd + (1 | pop) + (1 | cor.grp)
        
      } else {
        fmla <- gd ~ cd + (1 | pop)
        
      }
    }
  }
  
  
  
  # Pairs to include ---------------------------------------------------------
  pairs_to_include.file <- NULL
  keep <- pairs_to_include
  if(is.null(keep)){
    keep <-  rep(1, nrow(ID))
  }
  ID.keep <- ID
  covariates.keep <- covariates
  ZZ.keep <- ZZ
  response.keep <- response
  
  # if(!is.null(pop2ind)) {
  #   keep.vec <- expand.keep(pop2ind)
  #   response.keep <- response[keep == 1]
  #   covariates.keep <- covariates[keep == 1,]
  # }
  
  if (!is.null(pairs_to_include)) {
    
    keep <-  pairs_to_include
    
    df <- df[keep == 1,]
    ZZ.keep <- ZZ[,keep == 1]
    pop <- ID$pop1[keep == 1] ## Will need to be updated for pop2ind
    miss.pops <- as.character(pop)
    
    ## Reduce ZZ
    ZZ.keep <- ZZ.keep[rownames(ZZ.keep) %in% unique(miss.pops),]
    
    ## Reduce response & covariates
    response.keep <- response[keep == 1]
    covariates.keep <- covariates[keep == 1,]
    
    keep.id <- ID[keep == 1,]
    names(keep.id) <- c('mode', 'include')
    write.table(keep.id,
                file = paste0(td, "include_pairs.txt"),
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE)
    pairs_to_include.file <- paste0(td, "include_pairs.txt")
    ID.keep <- ID[keep == 1,]
    
  }
  
  ## Check slope
  if(!is.null(response)){
    cd <- c(dist(cs_coords))[keep==1]
    m <- lm(gd ~ cd, data = df)
    if(coef(m)[2] < 0){
      warning('Genetic distance decreases with distance. This is likely to result in a failed optimization.\nCheck your measure carefully and consider subtracting your values from 1 to reverse the relationship.')
    }
  }
  
  
  ##Return list  
  list(
    ID = ID.keep,
    ZZ = ZZ.keep,
    df = df,
    response = response.keep,
    covariates = covariates.keep,
    formula = fmla,
    CS_Point.File = CS_Point.File,
    Neighbor.Connect = Neighbor.Connect,
    n.Pops = n.Pops,
    pairs_to_include = pairs_to_include.file,
    pop2ind = pop2ind,
    nb = nb,
    parallel = parallel,
    cores = cores,
    solver = solver,
    precision = precision,
    JULIA_HOME = JULIA_HOME,
    write.files = write.files,
    write.criteria = write.criteria,
    silent = silent,
    Julia_link = Julia_link,
    scratch = scratch,
    keep = keep,
    response.all = response,
    ID.all = ID,
    ZZ.all = ZZ,
    covariates.all = covariates,
    rm.files = rm.files,
    opt.input = 'jl'
  )
}
