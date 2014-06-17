# About file
#' @docType package
#' @name ResistanceGA-package
#' @title  About this R package
#' @description  This package contains functions to optimize resistance surfaces using Genetic Algorithms. Continuous and categorical surfaces can be optimized, and multiple surfaces can be simultaneously optimized to create novel resistance surfaces.
#' @details  \tabular{ll}{
#'    Package: \tab ResistanceOptimization\cr
#'    Type: \tab Package\cr
#'    License: \tab GPL-3\cr
#'  }
#'  This package provides functions to prepare data and execute a number of functions to optimize continuous and categorical resistance surfaces using CIRCUITSCAPE and Genetic Algorithms within R. You must have CIRCUITSCAPE (4.0-Beta or higher) installed to run these functions. Output from functions in this package include: an AICc model rank table comparing all optimized surfaces, parameters that optimized each of the top models, coefficients from the fitted mixed effects models, plots of the fitted response curves, diagnostic plots of model fit, and Circuitscape outputs for each of the optimized resistance surfaces.
#'   
#'  ***This code can currently only be used with Windows machines due its use of the Circuitscape .exe file. Anyone interested in adapting the code to accommodate command-line execution on other platforms is welcome to do so.
#'  
#'   In order to use this package, you must have CIRCUITSCAPE v4.0 or greater installed.
#'   
#'   Official release: \url{http://www.circuitscape.org/downloads}
#' 
#' @import ggplot2 raster GA lme4 gdistance
#' @importFrom Matrix fac2sparse
#' @importFrom akima interp
#' @importFrom plyr arrange rbind.fill
#' @references Please cite: Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology 23:2402--2413. \href{http://goo.gl/RJb6Go}{Peterman et al.}
#' 
#' @author Bill Peterman \email{bill.peterman@@gmail.com}
#' 
NULL

#' Simulated resistance surfaces
#' 
#' A raster stack containing four raster surfaces
#' 
#' \itemize{
#'    \item continuous1. Continuous resistance surface
#'    \item continuous2. Continuous resistance surface
#'    \item categorical. 3-class categorical resistance surface
#'    \item feature. 2-class categorical (feature) resistance surface
#'    }
#' 
#' @docType data
#' @name resistance_surfaces
#' @format Raster stack object of length 4
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