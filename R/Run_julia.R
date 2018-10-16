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
#' @param JULIA_HOME Path to the folder containing the Julia binary (See Details). Only necessary to specify if \code{jl.inputs} are not specified.
#' #' @param Julia_link Specify whether R should connect to Julia using the 'JuliaCall' package or the 'XRJulia' package. Will Default to using 'JuliaCall' if not specified here or in \code{jl.inputs}
#' @return Vector of CIRCUITSCAPE resistance distances (lower half of resistance matrix) OR a full square distance matrix if `full.mat` = TRUE. Alternatively, a raster object of the cumulative current map can be returned when \code{CurrentMap = TRUE} and \code{output = "raster"}.
#' @usage Run_CS.jl(jl.inputs = NULL,
#' r,
#' CurrentMap = FALSE,
#' full.mat = FALSE,
#' EXPORT.dir = NULL,
#' output = "matrix",
#' CS_Point.File = NULL,
#' JULIA_HOME = NULL,
#' Julia_link = NULL)

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
           JULIA_HOME = NULL,
           Julia_link = NULL) {
    
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
                           JULIA_HOME = JULIA_HOME,
                           run_test = FALSE)
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
      EXPORT.dir <- paste0(tempdir(), "\\")
    } else {
      # if(CurrentMap == FALSE) {
      #   EXPORT.dir <- paste0(tempdir(), "\\")
      # } else {
      EXPORT.dir
      # }
    }
    
    if (CurrentMap == FALSE) {
      File.name <- names(R)
      MAP = "write_cum_cur_map_only = False"
      CURRENT.MAP = "write_cur_maps = False"
      
    } else {
      File.name <- names(R)
      MAP = "write_cum_cur_map_only = True"
      CURRENT.MAP = "write_cur_maps = 1"
    }
    
    ######
    
    if (max(R@data@values, na.rm = TRUE) > 1e6)
      R <-
        SCALE(R, 1, 1e6) # Rescale surface in case resistances are too high
    R <- reclassify(R, c(-Inf, 0, 1))
    
    temp_rast <- tempfile(pattern = "raster_", 
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
        filename = paste0(EXPORT.dir, File.name, '.asc'),
        overwrite = TRUE
      )
      
      temp_rast <- paste0(EXPORT.dir, File.name, '.asc')
      tmp.name <- File.name
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
      out <- read.delim(paste0(EXPORT.dir, tmp.name, "_resistances.out"), header = FALSE)[-1,-1]
      # out <- juliaGet(cs.out)[-1,-1] ## Slow!
    }
    
    
    if (output == "raster" & CurrentMap == TRUE) {
      rast <- raster(paste0(EXPORT.dir, tmp.name, "_cum_curmap.asc"))
      # NAME <- basename(rast@file@name)
      # NAME <- sub("^([^.]*).*", "\\1", NAME)
      names(rast) <- File.name
      (rast)
    } else {
      if(full.mat == FALSE) {
        cs.matrix <- lower(out)
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
      
      
      unlink.list <- list.files(EXPORT.dir, 
                                pattern = tmp.name,
                                all.files = TRUE,
                                full.names = TRUE)
      
      del.files <- sapply(unlink.list, unlink)
      
      return(cs.matrix)
    }
  } # End Julia function