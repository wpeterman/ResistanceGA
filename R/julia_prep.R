#' Prepare and bundle input CIRCUITSCAPE model parameters to run in Julia
#'
#' This function will prepare objects needed for running optimization functions
#'
#' @param n.Pops The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix). Not necessary if only executing Julia run.
#' @param CS_Point.File Provide a \code{\link[sp]{SpatialPoints}} object containing sample locations. Alternatively, specify the path to the Circuitscape formatted point file. See Circuitscape documentation for help.
#' @param covariates Data frame of additional covariates that you want included in the MLPE model during opitmization.
#' @param JULIA_HOME Path to the folder containing the Julia binary (See Details)
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
#' @param pairs_to_include Default is NULL. If you wish to use the advanced CIRCUITSCAPE setting mode to include or exclude certain pairs of sample locations, provide the path to the properly formatted "pairs_to_include.txt" file here. Currently only "include" method is supported.
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
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage jl.prep(n.Pops, 
#' response, 
#' CS_Point.File, 
#' covariates = NULL,
#' JULIA_HOME,
#' Neighbor.Connect, 
#' pairs_to_include, 
#' platform, 
#' parallel, 
#' cores,
#' cholmod,
#' precision, 
#' run_test,
#' write.files = NULL,
#' write.criteria = NULL,
#' silent = TRUE,
#' Julia_link = 'JuliaCall',
#' scratch = NULL,
#' rm.files = TRUE)
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
jl.prep <- function(n.Pops,
                    response = NULL,
                    CS_Point.File,
                    covariates = NULL,
                    JULIA_HOME = NULL,
                    Neighbor.Connect = 8,
                    pairs_to_include = NULL,
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
  
  # Determine if CIRCUITSCAPE package is installed
  # if(Julia_link == 'JuliaCall'){
  if(julia_installed_package("Circuitscape") == 'nothing') {
    julia_install_package('Circuitscape')
    julia_call('Pkg.test', "Circuitscape")
    
    
    # Run test ----------------------------------------------------------------
    
  } else {
    if(Julia_link == 'JuliaCall'){
      julia_library("Circuitscape")
    }
    
    if(run_test == TRUE) {
      print("Running test #1: Run Circuitscape from Julia")
      
      if(Sys.info()[['sysname']] == "Windows") {
        td <- paste0(tempdir(),"\\")
      } else {
        td <- paste0(tempdir(),"/")
      }
      write.table(samples, 
                  paste0(td,'samples.txt'), 
                  quote = FALSE,
                  sep = "\t",
                  row.names = FALSE,
                  col.names = FALSE)
      
      temp.ini <- tempfile(pattern = "", 
                           tmpdir = tempdir(),
                           fileext = ".ini")
      
      tmp.name <- basename(temp.ini) %>% strsplit(., '.ini') %>% unlist()
      
      writeRaster(resistance_surfaces$continuous,
                  paste0(td,
                         tmp.name, '.asc'),
                  overwirte = TRUE)
      
      if(!is.null(scratch)) {
        write.CS_4.0(BATCH = paste0(td, tmp.name, ".ini"),
                     OUT = paste0("output_file = ", scratch,"/", tmp.name, ".out"),
                     HABITAT = paste0("habitat_file = ", td, tmp.name, '.asc'),
                     LOCATION.FILE = paste0("point_file = ", td, 'samples.txt'),
                     PARALLELIZE = FALSE,
                     CORES = NULL,
                     solver = 'cholmod',
                     precision = FALSE,
                     silent = silent
        )
      } else {
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
      }
      
      if(Julia_link == 'JuliaCall'){
        out <- julia_call('compute', temp.ini)[-1,-1]
        
      } else {
        cs.jl <- RJulia()
        cs.jl$Using("Circuitscape")
        
        cs.out <- cs.jl$Call("compute", temp.ini) 
        out <- as.matrix(read.table(paste0(scratch, "/", tmp.name, "_resistances.out"),
                                    quote="\"", comment.char=""))[-1,-1]
        # out <- read.delim(paste0(scratch, "/", tmp.name, "_resistances.out"), header = FALSE)[-1,-1]
        # out <- juliaGet(cs.out)[-1,-1] ## SLOW!!!
      }
      
      
      
      if(dim(out)[1] == 25) {
        cat("\n"); cat("\n")
        
        cat("Test #1: Passed")
        
        cat("\n"); cat("\n")
        
      } else {
        stop("Test #1: Failed")
      }
      
      #   
      #   cat("Running test #2: Run Circuitscape from Julia in parallel")
      #   
      #   write.CS_4.0(BATCH = paste0(td, tmp.name, ".ini"),
      #                OUT = paste0("output_file = ", td, tmp.name, ".out"),
      #                HABITAT = paste0("habitat_file = ", td, tmp.name, '.asc'),
      #                LOCATION.FILE = paste0("point_file = ", td, 'samples.txt'),
      #                PARALLELIZE = TRUE,
      #                CORES = 2,
      #                solver = 'cholmod',
      #                precision = NULL
      #   )
      #   
      #   out <- julia_call('compute', temp.ini)[-1,-1]
      #   
      #   
      #   if(dim(out)[1] == 25) {
      #     cat("\n"); cat("\n")
      #     
      #     print("Test #2: Passed")
      #   } else {
      #     stop("Test #2: Failed")
      #   } 
      #   
      unlink.list <- list.files(td,
                                pattern = tmp.name,
                                all.files = TRUE,
                                full.names = TRUE)
      del.files <- sapply(unlink.list, unlink, force = TRUE)
      
      if(!is.null(scratch)){
        unlink.list2 <- list.files(scratch,
                                   pattern = tmp.name,
                                   all.files = TRUE,
                                   full.names = TRUE)
        del.files <- sapply(unlink.list2, unlink, force = TRUE)
      }
      
    } # End Julia setup
  } # End test
  
  
  # Format inputs -----------------------------------------------------------
  if(isTRUE(precision)) {
    precision <- 'single'
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
    td <- tempdir()
    # sp_file <- tempfile(pattern = "sample_pts_", 
    #                     tmpdir = tempdir(),
    #                     fileext = ".txt")
    
    # sp_name <- basename(sp_file) %>% strsplit(., '.txt') %>% unlist()
    
    site <- c(1:length(CS_Point.File))
    
    cs.txt <- data.frame(site, CS_Point.File)
    
    write.table(
      cs.txt,
      # file = sp_file,
      file = paste0(td, "\\samples_pts.txt"),
      col.names = F,
      row.names = F
    )
    
    # CS_Point.File <- sp_file
    CS_Point.File <- paste0(td, "\\samples_pts.txt")
  }
  
  if (grepl(".asc", x = CS_Point.File)) {
    CS_grid <- raster <- raster(CS_Point.File)
    CS_Point.txt <- rasterToPoints(x = CS_grid)
    site <- CS_Point.txt[, 3]
    CS_Point.txt <- data.frame(site, CS_Point.txt[, c(1, 2)])
    CS_Point.txt <- CS_Point.txt[order(site), ]
    CS_Point.File <- sub(".asc", ".txt", x = CS_Point.File)
    write.table(
      CS_Point.txt,
      file = CS_Point.File,
      col.names = F,
      row.names = F
    )
  }
  
  if (!is.null(response)) {
    TEST.response <- (is.vector(response) || ncol(response) == 1)
    if (TEST.response == FALSE) {
      stop("The object 'response' is not in the form of a single column vector")
    }
  }
  
  if (!is.null(pairs_to_include)) {
    if (!file.exists(pairs_to_include)) {
      stop("The specified pairs_to_include file does not exist")
    }
    toMatch <- c("min", "max")
    if (grep(
      read.table(
        file = pairs_to_include,
        header = F,
        sep = "\t"
      )[1, 1],
      pattern = paste(toMatch, collapse = "|")
    )) {
      # Function to make list of observations to include
      Min.MAX <-
        read.table(file = pairs_to_include,
                   header = F,
                   sep = "\t")[c(1, 2), c(1, 2)]
      MIN <- Min.MAX[which(Min.MAX[, 1] == "min"), 2]
      MAX <- Min.MAX[which(Min.MAX[, 1] == "max"), 2]
      PTI <-
        read.table(file = pairs_to_include,
                   header = F,
                   sep = "\t")[-c(1:3), ]
      site <- PTI[, 1]
      PTI <- as.matrix(PTI[, -1])
      
      p_t_i <- list()
      count <- 0
      for (i in 1:(ncol(PTI) - 1)) {
        for (j in (i + 1):ncol(PTI)) {
          if (PTI[j, i] >= MIN && PTI[j, i] <= MAX) {
            count <- count + 1
            p_t_i[[count]] <- data.frame(i, j)
          } # close if statement
        } # close j loop
      } # close i loop
      
      ID <- plyr::ldply(p_t_i, .fun = identity)
      colnames(ID) <- c("pop1", "pop2")
      #           ID<-arrange(tmp2,as.numeric(pop1),as.numeric(pop2))
      n1 <- table(ID$pop1)[[1]]
      p1 <- ID[n1, 1]
      p2 <- ID[n1, 2]
      ID[n1, 1] <- p2
      ID[n1, 2] <- p1
      ID$pop1 <- factor(ID$pop1)
      ID$pop2 <- factor(ID$pop2)
    } # close function
  } # close pairs to include statement
  
  # Make to-from population list
  if (!exists(x = "ID")) {
    ID <- To.From.ID(n.Pops)
  }
  suppressWarnings(ZZ <- ZZ.mat(ID))
  
  # if(!is.null(covariates)) {
  #   formula = NULL
  # } else {
  #   formula = lme4::lmer(response ~ . + (1|pop1), data = covariates) 
  # } mlpe_rga(response ~ . -pop + (1|pop), data=covariates)
  
  list(
    ID = ID,
    ZZ = ZZ,
    response = response,
    covariates = covariates,
    CS_Point.File = CS_Point.File,
    Neighbor.Connect = Neighbor.Connect,
    n.Pops = n.Pops,
    pairs_to_include = pairs_to_include,
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
    rm.files = rm.files
  )
}