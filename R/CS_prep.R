#' Prepare and bundle input CIRCUITSCAPE model parameters
#'
#' This function will prepare objects needed for running optimization functions
#'
#' @param n.Pops The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param CS_Point.File The path to the Circuitscape formatted point file. See Circuitscape documentation for help.
#' @param CS.program The path to the CIRCUITSCAPE executable file (cs_run.exe) on a Windows PC. If using a Linux or Mac system, provide the full path to the "csrun.py" file. See details below.
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
#' @param pairs_to_include Default is NULL. If you wish to use the advanced CIRCUITSCAPE setting mode to include or exclude certain pairs of sample locations, provide the path to the properly formatted "pairs_to_include.txt" file here. Currently only "include" method is supported.
#' @param platform What computing platform are you using ("pc", "other"). This code has only been tested on Windows PC!!!
#' @param parallel (Logical) If using Linux / Ubuntu, do you want to run CIRCUITSCAPE in parallel?
#' @param cores If using Linux / Ubuntu and `parallel = TRUE`, how many cores should be used for parallel processing?

#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage CS.prep(n.Pops, 
#' response, 
#' CS_Point.File, 
#' CS.program, 
#' Neighbor.Connect, 
#' pairs_to_include, 
#' platform, 
#' parallel, 
#' cores)
#' @details \code{CS.program} Example of path to CIRCUITSCAPE executible on Windows:
#'
#' '"C:/Program Files/Circuitscape/cs_run.exe"'
#'
#' ***NOTE: Double quotation used***
#' This is the current default for \code{CS.program}, but the directory may need to be changed depending upon your installation of CIRCUITSCAPE
#'
#' Linux
#' To call CIRCUITSCAPE from R on with Linux, first change file permissions from the command line terminal (shortcut: control + alt+ t) ` sudo chomod 755 /usr/local/bin/csrun.py `
#' Then specify \code{CS.program} as `csrun.py`
#'
#' Only with Linux, \code{parallel} can be set to \code{TRUE}, and the number of cores to run in parallel can be specified with \code{cores}.
#'
#' The Linux and Mac versions are in development. Please let me know if you encounter errors.

CS.prep <- function(n.Pops,
                    response = NULL,
                    CS_Point.File,
                    CS.program = '"C:/Program Files/Circuitscape/cs_run.exe"',
                    Neighbor.Connect = 8,
                    pairs_to_include = NULL,
                    platform = 'pc',
                    parallel = FALSE,
                    cores = NULL) {
  CS.exe_Test <- gsub("\"", "", CS.program)
  # Error messages
  if (!file.exists(CS_Point.File)) {
    stop("The specified CS_Point.File does not exist")
  }
  if (platform == 'pc') {
    if (!file.exists(gsub("\"", "", CS.program))) {
      stop("The specified path to 'cs_run.exe' is incorrect")
    }
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
  
  if (platform == 'pc') {
    parallel = FALSE
    cores = NULL
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
  list(
    ID = ID,
    ZZ = ZZ,
    response = response,
    CS_Point.File = CS_Point.File,
    CS.program = CS.program,
    Neighbor.Connect = Neighbor.Connect,
    n.Pops = n.Pops,
    platform = platform,
    pairs_to_include = pairs_to_include,
    parallel = parallel,
    cores = cores
  )
}