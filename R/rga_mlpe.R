# Run Mixed effects models using formula interface
#' Run maximum likelihood population effects mixed effects model (MLPE)
#'
#' Runs MLPE as detailed by Clarke et al. (2002). This is a general function for flexibly fitting MLPE models using the standard \code{lme4} formula interface
#'
#' @param formula \code{lme4} style mixed effects model equation.
#' @param data data frame containing vectors of values from the lower half of genetic/resistance distance matrices
#' @param REML Logical. If TRUE, mixed effects model will be fit using restricted maximum likelihood. Default = FALSE
#' @param ZZ The sparse matrix object for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @return A lmer object from the fitted model
#' @details An AIC value will only be returned if \code{REML = FALSE}. The random effect must be the population vector \code{pop1} generated from the function \code{\link[ResistanceGA]{To.From.ID}}.
#' @examples  
#' # Create square 'distance' matrices
#' y <- matrix(rnorm(16), 4)
#' x <- matrix(rnorm(16), 4)
#' 
#' # Create to-from object (4 populations sampled)
#' id <- To.From.ID(4)
#' 
#' # Create data frame
#' df <- data.frame(y = lower(y),
#'                  x = lower(x),
#'                  pop = id$pop1)
#' 
#' # Fit MLPE model
#' out <- mlpe_rga(formula = y ~ x + (1 | pop),
#'                 data = df)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage mlpe_rga(formula, data, REML = FALSE, ZZ = NULL)
#' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

mlpe_rga <-
  function(formula,
           data,
           REML = FALSE,
           ZZ = NULL) {
    
    
    if(class(formula) != 'formula') {
      formula <- as.formula(formula)
    }
    
    if(is.null(ZZ)) {
      obs <- 0.5 * (sqrt((8 * nrow(data)) + 1) + 1)
      ID <- To.From.ID(obs)
      ZZ <- ZZ.mat(ID)
    }
    
    # Fit model
    mod <-
      lFormula(formula,
               data = data,
               REML = REML)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun, mod)
    opt <- optimizeLmer(dfun)
    MOD <-
      (mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    return(MOD)
  }


# Make to-from population list
To.From.ID <- function(POPS) {
  tmp <- matrix(nrow = POPS, ncol = POPS)
  dimnames(tmp) <- list(1:POPS, 1:POPS)
  tmp2 <-
    as.data.frame(which(row(tmp) < col(tmp), arr.ind = TRUE))
  tmp2[[2]] <- dimnames(tmp)[[2]][tmp2$col]
  tmp2[[1]] <- dimnames(tmp)[[2]][tmp2$row]
  colnames(tmp2) <- c("pop1", "pop2")
  as.numeric(tmp2$pop1)
  as.numeric(tmp2$pop2)
  ID <- plyr::arrange(tmp2, as.numeric(pop1), as.numeric(pop2))
  p1 <- ID[POPS - 1, 1]
  p2 <- ID[POPS - 1, 2]
  ID[POPS - 1, 1] <- p2
  ID[POPS - 1, 2] <- p1
  ID$pop1 <- factor(ID$pop1)
  ID$pop2 <- factor(ID$pop2)
  return(ID)
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