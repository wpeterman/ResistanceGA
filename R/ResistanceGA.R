# About file
#' @docType package
#' @name ResistanceGA-package
#' @title  About this R package
#' @description  This package contains functions to optimize resistance surfaces using Genetic Algorithms. Continuous and categorical surfaces can be optimized, and multiple surfaces can be simultaneously optimized to create novel resistance surfaces.
#' @details  \tabular{ll}{
#'    Package: \tab ResistanceOptimization\cr
#'    Type: \tab Package\cr
#'    License: \tab >=GPL-2\cr
#'  }
#'  This package provides functions to prepare data and execute a number of functions to optimize continuous and categorical resistance surfaces using CIRCUITSCAPE and Genetic Algorithms within R. You must have CIRCUITSCAPE (4.0-Beta or higher) installed to run these functions. Output from functions in this package include: a summary table with AIC, AICc, conditional and marginal R2 values, and log-likelihood values for each optimized surface, parameters that optimized each of the top models, coefficients from the fitted mixed effects models, plots of the fitted response curves, diagnostic plots of model fit, and Circuitscape outputs for each of the optimized resistance surfaces. Resistance surfaces can also be optimized using least cost paths, which are implemented using the `gdistance` package.
#'   
#'  *** Use of this package to run CIRCUITSCAPE is limited to Windows machines due its use of the Circuitscape .exe file. Anyone interested in adapting the code to accommodate command-line execution on other platforms is welcome to do so. 
#'  
#'   In order to use this package with CIRCUITSCAPE, you must have CIRCUITSCAPE v4.0 or greater installed.
#'   
#'   Official release: \url{http://www.circuitscape.org/downloads}
#' 
#' @import raster GA lme4 ggplot2 gdistance
#' @importFrom utils combn menu
#' @importFrom ggExtra removeGrid ggMarginal
#' @importFrom Matrix fac2sparse
#' @importFrom plyr arrange rbind.fill ldply create_progress_bar progress_text 
#' @importFrom dplyr mutate group_by summarise filter tally left_join
#' @importFrom akima interp
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom magrittr "%>%"
#' @importFrom plyr arrange rbind.fill ldply
#' @importFrom spatstat as.im blur as.matrix.im
#' @importFrom grDevices dev.off tiff topo.colors
#' @importFrom graphics abline filled.contour par
#' @importFrom stats AIC lm logLik qqline qqnorm resid residuals runif as.formula sigma
#' @importFrom utils file_test read.csv read.delim read.table write.table
#' 
#' @references Please cite: 
#' Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology 23:2402--2413. \href{http://goo.gl/RJb6Go}{Peterman et al.}
#' 
#' Peterman, W. E. 2018. ResistanceGA: An R package for the optimization of resistance surfaces using genetic algorithms. Methods in Ecology and Evolution doi:10.1111/2041-210X.12984. \href{https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12984}{"MEE Publication"}
#' 
#' @author Bill Peterman \email{bill.peterman@@gmail.com}
#' 
NULL

#' Simulated resistance surfaces
#' 
#' A raster stack containing three raster surfaces
#' 
#' \itemize{#'    
#'    \item categorical. 3-class categorical resistance surface
#'    \item continuous. Continuous resistance surface
#'    \item feature. 2-class categorical (feature) resistance surface
#'    }
#' 
#' @docType data
#' @name resistance_surfaces
#' @format RasterStack object of length 3
#' @usage data(resistance_surfaces)
#' @keywords datasets
#' 
NULL

#' Example sample location file to run CIRCUITSCAPE
#' 
#' A data frame that can be saved as a .txt file for running examples in the vignette 
#' 
#' @docType data
#' @name samples
#' @format A 25 x 3 data frame
#' @usage data(samples)
#' @description  Sample file to be used with examples in the vignette
#' @keywords datasets
NULL