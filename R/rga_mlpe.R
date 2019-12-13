# Run Mixed effects models using formula interface
#' Run maximum likelihood population effects mixed effects model (MLPE)
#'
#' Runs MLPE as detailed by Clarke et al. (2002). This is a general function for flexibly fitting MLPE models using the standard \code{lme4} formula interface
#'
#' @param formula \code{lme4} style mixed effects model equation.
#' @param data data frame containing vectors of values from the lower half of genetic/resistance distance matrices
#' @param REML Logical. If TRUE, mixed effects model will be fit using restricted maximum likelihood. Default = FALSE
#' @param ZZ The sparse matrix object for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @param keep A vector consisting of 1 (keep) or 0 (drop) for each pairwise observation in \code{data}. An option if you do not want to assess all pairiwse observations in the MLPE model.
#' @return A lmer object from the fitted model
#' @details An AIC value will only be returned if \code{REML = FALSE}. The random effect must be the population vector \code{pop1} generated from the function \code{\link[ResistanceGA]{To.From.ID}}.
#' 
#' A generalized MLPE model can be fit if an appropriate \code{\link[stats]{family}}. It is not possible to use REML when fitting a generalized model.
#' 
#' @examples  
#' # Create square 'distance' matrices
#' y <- matrix(rnorm(25), 5)
#' y.pois <- matrix(rpois(25, 5), 5)
#' x <- matrix(rnorm(25), 5)
#' 
#' # Create to-from object (5 populations sampled)
#' id <- To.From.ID(5)
#' 
#' # Create data frame
#' df <- data.frame(y = lower(y),
#'                  y.pois = lower(y.pois),
#'                  x = lower(x),
#'                  pop = id$pop1)
#' 
#' # Fit MLPE model
#' out <- mlpe_rga(formula = y ~ x + (1 | pop),
#'                 REML = TRUE,
#'                 data = df)
#'                 
#' # Fit generalized MLPE model
#' out.pois <- mlpe_rga(formula = y.pois ~ x + (1 | pop),
#'                      family = poisson,
#'                      data = df)
#'                                 
#' # Fit model with only select pairs
#'  keep <- c(1,0,1,1,1,1,0,0,1,0)
#'  
#'  out.select <- mlpe_rga(formula = y ~ x + (1 | pop),
#'                         data = df,
#'                         keep = keep)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage mlpe_rga(formula,
#'                 data, 
#'                 REML = FALSE, 
#'                 ZZ = NULL,
#'                 keep = NULL)
#' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

mlpe_rga <-
  function(formula,
           data,
           REML = FALSE,
           ZZ = NULL,
           keep = NULL,
           ...) {
    
    
    if(class(formula) != 'formula') {
      formula <- as.formula(formula)
    }
    
    if(is.null(ZZ)) {
      obs <- 0.5 * (sqrt((8 * nrow(data)) + 1) + 1)
      ID <- To.From.ID(obs)
      ZZ <- ZZ.mat(ID)
    }
    
    if(!is.null(keep)) {
      data <- data[keep == 1,]
      ZZ <- ZZ[,keep == 1]
      pop <- ID$pop1[keep == 1]
      miss.pops <- as.character(pop)
      
      ## Reduce ZZ
      ZZ <- ZZ[rownames(ZZ) %in% unique(miss.pops),]
    }
    
    
    # Fit merMod --------------------------------------------------------------
    
    
    # > glmer -----------------------------------------------------------------
    args <- list(...)
    
    if(any("family" %in% names(args))) {
      if(isTRUE(REML)) {
        cat("REML will be ignored when fitting a generalized MLPE model.")
      }
      mod <-
        glFormula(formula,
                  data = data,
                  ...)
      mod$reTrms$Zt <- ZZ
      dfun <- do.call(mkGlmerDevfun, mod)
      opt <- optimizeGlmer(dfun)
      
    } else {
      
      # > lmer ---------------------------------------------------------------
      
      mod <-
        lFormula(formula,
                 data = data,
                 REML = REML)
      mod$reTrms$Zt <- ZZ
      dfun <- do.call(mkLmerDevfun, mod)
      opt <- optimizeLmer(dfun)
      
    }
    
    MOD <-
      (mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    
    return(MOD)
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


lower <- function(matrix) {
  if (is.vector(matrix) == TRUE ||
      dim(matrix)[1] != dim(matrix)[2]) {
    warning("Must provide square distance matrix with no column or row names")
  }
  lm <- matrix[lower.tri(matrix)]
  return(lm)
}