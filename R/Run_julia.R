###############################################################################
############ RUN JULIA ############
###############################################################################
#' Run Julia version of CIRCUITSCAPE from R
#'
#' Execute julia CS from R
#'
#' @param jl.inputs Object created from running \code{\link[ResistanceGA]{jl.prep}} function
#' @param r Accepts two types of inputs. Provide either the path to the resistance surface file (.asc) or specify an R RasterLayer object
#' @param CurrentMap Logical. If TRUE, the cumulative current resistance map will be generated during the CS run (Default = FALSE)
#' @param full.mat Logical (Default = FALSE). If TRUE, the full distance matrix will be generated as an R object, rather than just the lower half of the distance matrix.
#' @param EXPORT.dir Directory where CS results will be written. Only specify if Circuitscape current map outputs are requested. It is critical that there are NO SPACES in the directory, as this may cause the function to fail.
#' @param output Specifiy either "matrix" or "raster". "matrix" will return the lower half of the pairwise resistance matrix (default), while "raster" will return a \code{raster} object of the cumulative current map. The raster map can only be returned if \code{CurrentMap=TRUE}
#' @param CS_Point.File Provide a \code{\link[sp]{SpatialPoints}} object containing sample locations. Alternatively, specify the path to the Circuitscape formatted point file. See Circuitscape documentation for help. Only necessary to specify if \code{jl.inputs} are not specified.
#' @param pairs_to_include Default is NULL. If you wish to use the advanced CIRCUITSCAPE setting mode to include or exclude certain pairs of sample locations, provide a vector consisting of 1 (keep) or 0 (drop) for each pairwise observation (see example). This is an option if you do not want to assess all pairiwse observations. Only necessary to specify if \code{jl.inputs} are not specified.
#' @param parallel (Logical; Default = FALSE) Do you want to run CIRCUITSCAPE in parallel? Only necessary to specify if \code{jl.inputs} are not specified.
#' @param cores If `parallel = TRUE`, how many cores should be used for parallel processing? Only necessary to specify if \code{jl.inputs} are not specified.
#' @param cholmod (Logical; Default = TRUE). Should the cholmod solver be used? See details. Only necessary to specify if \code{jl.inputs} are not specified.
#' @param precision (Logical; Default = FALSE). Should experimental single precision method be used? See details. Only necessary to specify if \code{jl.inputs} are not specified.
#' @param JULIA_HOME Path to the folder containing the Julia binary (See Details). Only necessary to specify if \code{jl.inputs} are not specified.
#' #' @param Julia_link Specify whether R should connect to Julia using the 'JuliaCall' package or the 'XRJulia' package. Will Default to using 'JuliaCall' if not specified here or in \code{jl.inputs}
#' @param rm.files Should all temporary files be removed after Julia run (Default = TRUE). Only necessary to specify if \code{jl.inputs} are not specified.
#' @param scratch Scratch directory for use if write access is limited. Must be specified if raster results are desired for outputs.
#' @return Vector of CIRCUITSCAPE resistance distances (lower half of resistance matrix) OR a full square distance matrix if `full.mat` = TRUE. Alternatively, a raster object of the cumulative current map can be returned when \code{CurrentMap = TRUE} and \code{output = "raster"}. The `full.mat` cannot be requested if only select pairs are being analyzed.
#' 
#' @examples 
#' ## To do
#' 
#' @usage Run_CS.jl(jl.inputs = NULL,
#' r,
#' CurrentMap = FALSE,
#' full.mat = FALSE,
#' EXPORT.dir = NULL,
#' output = "matrix",
#' CS_Point.File = NULL,
#' pairs_to_include = NULL,
#' parallel = FALSE, 
#' cores = NULL,
#' cholmod = TRUE,
#' precision = FALSE, 
#' JULIA_HOME = NULL,
#' Julia_link = NULL,
#' rm.files = TRUE,
#' scratch = NULL)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Run_CS.jl <-
  function(jl.inputs = NULL,
           r,
           CurrentMap = FALSE,
           full.mat = FALSE,
           EXPORT.dir = NULL,
           output = "matrix",
           CS_Point.File = NULL,
           pairs_to_include = NULL,
           parallel = FALSE,
           cores = NULL,
           cholmod = TRUE,
           precision = FALSE,
           JULIA_HOME = NULL,
           Julia_link = NULL,
           rm.files = TRUE,
           scratch = NULL) {
    
    if(!is.null(jl.inputs$rm.files)) {
      rm.files <- jl.inputs$rm.files
    }
    
    if(is.null(Julia_link) & !is.null(jl.inputs)) {
      Julia_link <- jl.inputs$Julia_link
    } else if(is.null(Julia_link) & is.null(jl.inputs)) {
      Julia_link <- 'JuliaCall'
    }
    
    if(Julia_link == 'XRJulia') {
      JULIA_HOME <- findJulia()
    } 
    
    if(is.null(jl.inputs) & is.null(JULIA_HOME)) {
      stop("Specify either `jl.inputs` or `JULIA_HOME` and `CS_Point.File`!")
    }
    
    if((is.null(CS_Point.File) | is.null(JULIA_HOME)) & is.null(jl.inputs)) {
      stop("Both `JULIA_HOME` and `CS_Point.File` must be specified!")
    }
    
    if(is.null(jl.inputs)) {
      jl.inputs <- jl.prep(n.Pops = length(CS_Point.File),
                           CS_Point.File = CS_Point.File,
                           pairs_to_include = pairs_to_include,
                           JULIA_HOME = JULIA_HOME,
                           scratch = scratch,
                           Julia_link = Julia_link,
                           parallel = parallel,
                           cores = cores,
                           cholmod = cholmod,
                           precision = precision,
                           run_test = FALSE)
    }
    
    scratch <- jl.inputs$scratch
    
    if(isTRUE(full.mat) & !is.null(jl.inputs$pairs_to_include)) {
      warning("The full matrix generated will have -1 values. These indicate pairs that were excluded from the analysis.")
    }
    
    
    if(!is.null(jl.inputs)) {
      if(is.null(jl.inputs$write.files)) {
        write.files <- NULL
      } else {
        write.files <- jl.inputs$write.files
      }
      
      if(is.null(jl.inputs$write.criteria)) {
        write.criteria <- NULL
      } else {
        write.criteria <- jl.inputs$write.criteria
      }
      
      JULIA_HOME <- jl.inputs$JULIA_HOME
    }
    
    if(Julia_link == 'JuliaCall') {
      julia_setup(JULIA_HOME = JULIA_HOME)
      julia_library("Circuitscape")
    } else {
      juliaEval("using Circuitscape")
    }
    
    
    if(!exists('r')) {
      stop("Missing variable r: Please specify a 'RasterLayer' object or provide the path to a raster file!")
    }
    
    if (class(r)[1] != 'RasterLayer') {
      R <- raster(r)
      File.name <- basename(r)
      File.name <- sub(".asc", "", File.name)
      names(R) <- File.name
    } else {
      R = r
      
    }
    
    if(is.null(EXPORT.dir)) {   
      if(Sys.info()[['sysname']] == "Windows") {
        EXPORT.dir <- paste0(tempdir(),"\\")
      } else {
        EXPORT.dir <- paste0(tempdir(),"/")
      }
    }
    
    if (CurrentMap == FALSE) {
      File.name <- names(R)
      MAP = "write_cum_cur_map_only = False"
      CURRENT.MAP = "write_cur_maps = 0"
      
    } else {
      File.name <- names(R)
      MAP = "write_cum_cur_map_only = True"
      CURRENT.MAP = "write_cur_maps = True"
    }
    
    ######
    
    if (max(R@data@values, na.rm = TRUE) > 1e6)
      R <-
      SCALE(R, 1, 1e6) # Rescale surface in case resistances are too high
    R <- reclassify(R, c(-Inf, 0, 1))
    
    rm.rast <- temp_rast <- tempfile(pattern = "raster_", 
                                     tmpdir = tempdir(),
                                     fileext = ".asc") 
    
    tmp.name <- basename(temp_rast) %>% strsplit(., '.asc') %>% unlist()
    
    if(CurrentMap == FALSE) {
      writeRaster(
        x = R,
        filename = temp_rast,
        overwrite = TRUE
      )
    } else {
      writeRaster(
        x = R,
        # filename = paste0(EXPORT.dir, File.name, '.asc'),
        filename = temp_rast,
        overwrite = TRUE
      )
      
      # temp_rast <- paste0(EXPORT.dir, File.name, '.asc')
      # tmp.name <- File.name
    }    
    
    # Modify and write Circuitscape.ini file
    if(jl.inputs$Neighbor.Connect == 4) {
      connect <- "True"
    } else {
      connect <- "False"
    }
    
    if (is.null(jl.inputs$pairs_to_include)) {
      PAIRS_TO_INCLUDE <-
        paste0("included_pairs_file = (Browse for a file with pairs to include or exclude)")
      PAIRS <- paste0("use_included_pairs = False")
    } else {
      PAIRS_TO_INCLUDE <-
        paste0("included_pairs_file = ", jl.inputs$pairs_to_include)
      PAIRS <- paste0("use_included_pairs = True")
    }
    if(!is.null(scratch)) {
      write.CS_4.0(
        BATCH = paste0(EXPORT.dir, tmp.name, ".ini"),
        OUT = paste0("output_file = ", scratch, "/", tmp.name, ".out"),
        HABITAT = paste0("habitat_file = ", temp_rast),
        LOCATION.FILE = paste0("point_file = ", jl.inputs$CS_Point.File),
        CONNECTION = paste0("connect_four_neighbors_only =", connect),
        MAP = MAP,
        CURRENT.MAP = CURRENT.MAP,
        PAIRS_TO_INCLUDE = PAIRS_TO_INCLUDE,
        PAIRS = PAIRS,
        PARALLELIZE = jl.inputs$parallel,
        CORES = jl.inputs$cores,
        solver = jl.inputs$solver,
        precision = jl.inputs$precision,
        silent = jl.inputs$silent
      )
    } else {
      write.CS_4.0(
        BATCH = paste0(EXPORT.dir, tmp.name, ".ini"),
        OUT = paste0("output_file = ", EXPORT.dir, tmp.name, ".out"),
        HABITAT = paste0("habitat_file = ", temp_rast),
        LOCATION.FILE = paste0("point_file = ", jl.inputs$CS_Point.File),
        CONNECTION = paste0("connect_four_neighbors_only =", connect),
        MAP = MAP,
        CURRENT.MAP = CURRENT.MAP,
        PAIRS_TO_INCLUDE = PAIRS_TO_INCLUDE,
        PAIRS = PAIRS,
        PARALLELIZE = jl.inputs$parallel,
        CORES = jl.inputs$cores,
        solver = jl.inputs$solver,
        precision = jl.inputs$precision,
        silent = jl.inputs$silent
      )
    }
    
    
    # Run CIRCUITSCAPE.jl -----------------------------------------------------
    rt <- NULL
    
    if(Julia_link == 'JuliaCall') {
      if(!is.null(write.criteria)) {
        t1 <- proc.time()[3]
        out <- julia_call('compute', paste0(EXPORT.dir, tmp.name, ".ini"))[-1,-1]
        rt <- proc.time()[3] - t1
      } else {
        out <- julia_call('compute', paste0(EXPORT.dir, tmp.name, ".ini"))[-1,-1]
      }
    } else { # use XRJulia
      cs.jl <- RJulia()
      cs.jl$Using("Circuitscape")
      
      ini.file <- paste0(EXPORT.dir, tmp.name, ".ini")
      cs.out <- cs.jl$Call("compute", ini.file) 
      Sys.sleep(0.5)
      out <- as.matrix(read.table(paste0(scratch, "/", tmp.name, "_resistances.out"),
                                  quote="\"", comment.char=""))[-1,-1]
      # out <- read.delim(paste0(scratch, "/", tmp.name, "_resistances.out"), header = FALSE)[-1,-1]
      # out <- juliaGet(cs.out)[-1,-1] ## Slow!
    }
    
    
    if (output == "raster" & CurrentMap == TRUE) {
      if(EXPORT.dir == paste0(tempdir(), "\\") & is.null(scratch)) {
        stop("Specify an `EXPORT.dir` or `scratch` directory to write CIRCUITSCAPE results to")
      } 
      
      if(EXPORT.dir != paste0(tempdir(), "\\")) {
        rast <- raster(paste0(EXPORT.dir, "/", tmp.name, "_cum_curmap.asc"))
      } else {
        rast <- raster(paste0(scratch, "/", tmp.name, "_cum_curmap.asc"))
      }
      
      # NAME <- basename(rast@file@name)
      # NAME <- sub("^([^.]*).*", "\\1", NAME)
      names(rast) <- File.name
      
      ## Remove files
      if(isTRUE(rm.files)) {
        unlink.list <- list.files(EXPORT.dir, 
                                  pattern = tmp.name,
                                  all.files = TRUE,
                                  full.names = TRUE)
        unlink.list2 <- list.files(EXPORT.dir, 
                                   pattern = basename(temp_rast),
                                   all.files = TRUE,
                                   full.names = TRUE)
        
        del.files <- sapply(unlink.list, unlink, force = TRUE)
        del.files <- sapply(unlink.list2, unlink, force = TRUE)
        unlink(rm.rast, force = TRUE)
        
        if(!is.null(scratch)){
          unlink.list2 <- list.files(scratch,
                                     pattern = tmp.name,
                                     all.files = TRUE,
                                     full.names = TRUE)
          del.files <- sapply(unlink.list2, unlink, force = TRUE)
        }
      }
      
      return(rast)
      
    } else {
      if(full.mat == FALSE) {
        cs.matrix <- lower(out)
        cs.matrix <- cs.matrix[cs.matrix != -1]
      } else {
        cs.matrix <- out
      }
      # Replace NA with 0...a workaround for errors when two points fall within the same cell.
      if (any(is.na(cs.matrix))) {
        cs.matrix[is.na(cs.matrix)] <- 0
        message(
          " NA values generated by CIRCUITSCAPE \n Check point file to see if multiple points share the same raster cell!"
        )
      }
      
      if(!is.null(write.files)) {
        if(!is.null(rt)) {
          if(rt > write.criteria){
            copy.files <- list(paste0(EXPORT.dir, tmp.name, ".ini"),
                               temp_rast)
            file.copy(copy.files, write.files)
          }
        } else {
          copy.files <- list(paste0(EXPORT.dir, tmp.name, ".ini"),
                             temp_rast)
          file.copy(copy.files, write.files)
        }
      }
      
      ## Remove files
      if(isTRUE(rm.files)) {
        unlink.list <- list.files(EXPORT.dir, 
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
          unlink(rm.rast, force = TRUE)
          
        }
      }
      
      asc.files <- list.files(EXPORT.dir, 
                              pattern = '.asc',
                              all.files = TRUE,
                              full.names = TRUE)
      
      make.times <- file.mtime(asc.files)
      
      time.diff <- Sys.time() - make.times
      
      unlink(asc.files[time.diff > 30], force = TRUE)
      
      return(cs.matrix)
    }
  } # End Julia function