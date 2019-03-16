###############################################################################
############ RUN CIRCUITSCAPE ############
###############################################################################
#' Run CIRCUITSCAPE in R
#'
#' Execute CS from R
#'
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param r Accepts two types of inputs. Provide either the path to the resistance surface file (.asc) or specify an R RasterLayer object
#' @param CurrentMap Logical. If TRUE, the cumulative current resistance map will be generated during the CS run (Default = FALSE)
#' @param full.mat Logical (Default = FALSE). If TRUE, the full distance matrix will be generated as an R object, rather than just the lower half of the distance matrix.
#' @param EXPORT.dir Directory where CS results should be written.Only specify if Circuitscape current map outputs are requested. It is critical that there are NO SPACES in the directory, as this may cause the function to fail.
#' @param output Specifiy either "matrix" or "raster". "matrix" will return the lower half of the pairwise resistance matrix (default), while "raster" will return a \code{raster} object of the cumulative current map. The raster map can only be returned if \code{CurrentMap=TRUE}
#' @param hidden Logical. If TRUE (Default), then no output from CIRCUITSCAPE will be printed to the console. Only set to FALSE when trying to troubleshoot/debug code.
#' @return Vector of CIRCUITSCAPE resistance distances (lower half of "XXX_resistances.out") OR a full square distance matrix if `full.mat` = TRUE. Alternatively, a raster object of the cumulative current map can be returned when \code{CurrentMap = TRUE} and \code{output = "raster"}.
#' @usage Run_CS(CS.inputs, 
#' r, 
#' CurrentMap = FALSE, 
#' full.mat = FALSE,
#' EXPORT.dir = NULL, 
#' output = "matrix", 
#' hidden = TRUE)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Run_CS <-
  function(CS.inputs,
           r,
           CurrentMap = FALSE,
           full.mat = FALSE,
           EXPORT.dir = NULL,
           output = "matrix",
           hidden = TRUE) {
    
    save.dir <- EXPORT.dir
    
    if(!exists('r')) {
      stop("Missing variable r: Please specify a 'RasterLayer' object or provide the path a raster file!")
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
      save.dir <- EXPORT.dir
    } else {
      save.dir <- EXPORT.dir
      EXPORT.dir <- paste0(tempdir(), "\\")
    }
    
    if (CurrentMap == FALSE) {
      File.name <- names(R)
      MAP = "write_cum_cur_map_only = False"
      CURRENT.MAP = "write_cur_maps = 0"
      
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
        filename = paste0(save.dir, File.name, '.asc'),
        overwrite = TRUE
      )
      
      temp_rast <- paste0(save.dir, File.name, '.asc')
      tmp.name <- File.name
    }
    
    # Modify and write Circuitscape.ini file
    if(CS.inputs$Neighbor.Connect == 4){
      connect <- "True"
    } else {
      connect <- "False"
    }
    
    if (is.null(CS.inputs$pairs_to_include)) {
      PAIRS_TO_INCLUDE <-
        paste0("included_pairs_file = (Browse for a file with pairs to include or exclude)")
      PAIRS <- paste0("use_included_pairs = False")
    } else {
      PAIRS_TO_INCLUDE <-
        paste0("included_pairs_file = ", CS.inputs$pairs_to_include)
      PAIRS <- paste0("use_included_pairs = True")
    }
    
    write.CS_4.0(
      BATCH = paste0(EXPORT.dir, tmp.name, ".ini"),
      OUT = paste0(paste0("output_file = ", EXPORT.dir), tmp.name, ".out"),
      HABITAT = paste0("habitat_file = ", temp_rast),
      LOCATION.FILE = paste0("point_file = ", CS.inputs$CS_Point.File),
      CONNECTION = paste0("connect_four_neighbors_only=", connect),
      MAP = MAP,
      CURRENT.MAP = CURRENT.MAP,
      PAIRS_TO_INCLUDE = PAIRS_TO_INCLUDE,
      PAIRS = PAIRS
    )
    
    # Execute CS ------------------------------------------------------------------
    
    # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
    # hidden = TRUE
    # Run Circuitscape
    if (grepl("cs_run.exe", CS.inputs$CS.program)) {
      CS.ini <- paste0(EXPORT.dir, tmp.name, ".ini")
      CS.Run.output <-
        system(paste(CS.inputs$CS.program, CS.ini), hidden)
    } else {
      CS.ini <- paste0(EXPORT.dir, tmp.name, ".ini")
      CS.Run.output <-
        system(paste(CS.inputs$CS.program, CS.ini), hidden)
    }
    
    
    
    
    #########################################
    # Run mixed effect model on each Circuitscape effective resistance
    
    CS.results <- paste0(EXPORT.dir, tmp.name, "_resistances.out")
    
    if (output == "raster" & CurrentMap == TRUE) {
      rast <- raster(paste0(EXPORT.dir, tmp.name, "_cum_curmap.asc"))
      NAME <- basename(rast@file@name)
      NAME <- sub("^([^.]*).*", "\\1", NAME)
      names(rast) <- NAME
      (rast)
    } else {
      if(full.mat == FALSE) {
        cs.matrix <- read.matrix(CS.results)
      } else {
        cs.matrix <- read.matrix2(CS.results)
      }
      # Replace NA with 0...a workaround for errors when two points fall within the same cell.
      if (any(is.na(cs.matrix))) {
        cs.matrix[is.na(cs.matrix)] <- 0
        message(
          " NA values generated by CIRCUITSCAPE \n Check point file to see if multiple points share the same raster cell!"
        )
      }
      
      unlink.list <- list.files(EXPORT.dir, 
                                pattern = tmp.name,
                                all.files = TRUE,
                                full.names = TRUE)
      
      del.files <- sapply(unlink.list, unlink)
      
      return(cs.matrix)
    }
  }
################################

Run_CS2 <-
  function(CS.inputs,
           r,
           EXPORT.dir = NULL,
           full.mat = FALSE) {

    # Modify and write Circuitscape.ini file
    #############################################################################################
    if(is.null(EXPORT.dir)) {
      EXPORT.dir <- paste0(tempdir(), "\\")
    } else {
      if(CurrentMap == FALSE) {
        EXPORT.dir <- paste0(tempdir(), "\\")
      } else {
        EXPORT.dir
      }
    }
    
    temp_rast <- tempfile(pattern = "raster_", 
                          tmpdir = tempdir(),
                          fileext = ".asc") 
    
    File.name <- basename(temp_rast) %>% strsplit(., '.asc') %>% unlist()
    
      writeRaster(
        x = r,
        filename = temp_rast,
        overwrite = TRUE
      )
    
    
    ifelse(CS.inputs$Neighbor.Connect == 4,
           connect <- "True",
           connect <- "False")
    
    if (is.null(CS.inputs$pairs_to_include)) {
      PAIRS_TO_INCLUDE <-
        paste0("included_pairs_file = (Browse for a file with pairs to include or exclude)")
      PAIRS <- paste0("use_included_pairs = False")
    } else {
      PAIRS_TO_INCLUDE <-
        paste0("included_pairs_file = ", CS.inputs$pairs_to_include)
      PAIRS <- paste0("use_included_pairs = True")
    }
    
    write.CS_4.0(
      BATCH = paste0(EXPORT.dir, File.name, ".ini"),
      OUT = paste0(paste0("output_file = ", EXPORT.dir), File.name, ".out"),
      HABITAT = paste0("habitat_file = ", paste0(EXPORT.dir, File.name, ".asc")),
      LOCATION.FILE = paste0("point_file = ", CS.inputs$CS_Point.File),
      CONNECTION = paste0("connect_four_neighbors_only=", connect),
      MAP = "write_cum_cur_map_only = False",
      CURRENT.MAP = "write_cur_maps = 0",
      PAIRS_TO_INCLUDE = PAIRS_TO_INCLUDE,
      PAIRS = PAIRS
    )
    ##########################################################################################
    # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
    hidden = TRUE
    # Run Circuitscape
    if (grepl("cs_run.exe", CS.inputs$CS.program)) {
      CS.ini <- paste0(EXPORT.dir, File.name, ".ini")
      CS.Run.output <-
        system(paste(CS.inputs$CS.program, CS.ini), hidden)
    } else {
      CS.ini <- paste0(EXPORT.dir, File.name, ".ini")
      CS.Run.output <-
        system(paste(CS.inputs$CS.program, CS.ini), hidden)
    }
    
    CS.results <- paste0(EXPORT.dir, File.name, "_resistances.out")
    
    if(full.mat == FALSE) {
      cs.matrix <- read.matrix(CS.results)
    } else {
      cs.matrix <- read.matrix2(CS.results)
    }
    
    # Replace NA with 0...a workaround for errors when two points fall within the same cell.
    if (any(is.na(cs.matrix))) {
      cs.matrix[is.na(cs.matrix)] <- 0
      message(
        " NA values generated by CIRCUITSCAPE \n Check point file to see if multiple points share the same raster cell!"
      )
    }
    
    unlink.list <- list.files(EXPORT.dir, 
                              pattern = tmp.name,
                              all.files = TRUE,
                              full.names = TRUE)
    
    del.files <- sapply(unlink.list, unlink)
    return(cs.matrix)
  }