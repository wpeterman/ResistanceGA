# About file
#' @docType package
#' @name ResistanceGA-package
#' @title  About this R package
#' @description  This package contains functions to optimize resistance surfaces using Genetic Algorithms. Continuous and categorical surfaces can be optimized, and multiple surfaces can be simultaneously optimized to create novel resistance surfaces.
#' @details  \tabular{ll}{
#'    Package: \tab ResistanceGA\cr
#'    Type: \tab Package\cr
#'    License: \tab >=GPL-2\cr
#'  }
#'  This package provides functions to prepare data and execute a number of functions to optimize continuous and categorical resistance surfaces using CIRCUITSCAPE, \href{https://github.com/Circuitscape/Circuitscape.jl}{CIRCUITSCAPE written in Julia} and Genetic Algorithms within R. You must have CIRCUITSCAPE (4.0-Beta or higher) or Julia installed to run these functions. The continued development and support of functions is primarily occurring with `gdistsance` and Julia implementations. Output from functions in this package include: a summary table with AIC, AICc, conditional and marginal R2 values, and log-likelihood values for each optimized surface, parameters that optimized each of the top models, coefficients from the fitted mixed effects models, plots of the fitted response curves, diagnostic plots of model fit, and Circuitscape outputs for each of the optimized resistance surfaces. Resistance surfaces can also be optimized using least cost paths, which are implemented using the `gdistance` package.
#'   
#'  *** Use of this package to run CIRCUITSCAPE is limited to Windows machines due its use of the Circuitscape .exe file. However, Julia can be installed on any operating system, making this a more versatile option. 
#'  
#' 
#' @import raster GA lme4 ggplot2 gdistance 
#' @importFrom utils combn menu
#' @importFrom ggExtra removeGrid ggMarginal
#' @importFrom Matrix fac2sparse drop0
#' @importFrom plyr arrange rbind.fill ldply create_progress_bar progress_text 
#' @importFrom dplyr mutate group_by summarise filter tally left_join dense_rank
#' @importFrom akima interp
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom plyr arrange rbind.fill ldply
#' @importFrom spatstat as.im blur as.matrix.im
#' @importFrom spdep dnearneigh nb2mat
#' @importFrom grDevices dev.off tiff topo.colors
#' @importFrom graphics abline filled.contour par
#' @importFrom stats AIC lm logLik qqline qqnorm resid residuals runif as.formula sigma
#' @importFrom utils file_test read.csv read.delim read.table write.table
#' @importFrom JuliaCall julia_setup julia_library julia_call julia_installed_package julia_install_package
#' @importFrom XRJulia findJulia juliaEval RJulia juliaGet
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %do% %dopar% foreach
#' 
#' @references Please cite: 
#' Peterman, W.E., G.M. Connette, R.D. Semlitsch, and L.S. Eggert. 2014. Ecological resistance surfaces predict fine-scale genetic differentiation in a terrestrial woodland salamander. Molecular Ecology 23:2402--2413. \href{http://goo.gl/RJb6Go}{Peterman et al.}
#' 
#' Peterman, W. E. 2018. ResistanceGA: An R package for the optimization of resistance surfaces using genetic algorithms. Methods in Ecology and Evolution 9:1638--1647. \href{https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12984}{"MEE Publication"}
#' 
#' @author Bill Peterman \email{bill.peterman@@gmail.com}
#' 
NULL


#' Data file. Simulated true resistance surfaces
#' 
#' A raster stack containing three raster surfaces generated from the \code{RandomFields} package. These surfaces were used with individual-based genetic simulations conducted using \code{PopGenReport}.
#' 
#' \itemize{    
#'    \item cat_true. 4-class categorical resistance surface. Derived from reclass: 1 = 1, 2 = 75, 3 = 200, 4 = 25, 5 = 1.
#'    \item cont_true. Continuous resistance surface. Derived from inverse monomolecular transformation with a shape parameter of 2.75, and maximum resistance parameter of 100.
#'    \item multi_true. A composite surface created by summing transformed categorical and continuous surfaces, rescaled to a minimum value of 1.
#'    }
#' 
#' @docType data
#' @name raster_true
#' @format RasterStack object of length 3
#' @usage data(raster_true)
#' @keywords datasets
#' 
NULL


#' Data file. Simulated original resistance surfaces
#' 
#' A raster stack containing two raster surfaces generated from the \code{RandomFields} package. These surfaces are the raw, original surfaces (prior to transformation) used with individual-based genetic simulations conducted using \code{PopGenReport}. These are the inputs for optimization examples.
#' 
#' \itemize{    
#'    \item cat_orig. 5-class categorical resistance surface with equal representation of each categorcal level.
#'    \item cont_orig. Continuous resistance surface.
#'    }
#' 
#' @docType data
#' @name raster_orig
#' @format RasterStack object of length 2
#' @usage data(raster_orig)
#' @keywords datasets
#' 
NULL


#' Data file. A list of three \code{sp} objects 
#' 
#' 50 Point locations in each list that correspond to sample points (populations) used in genetic simulations
#' 
#'  \itemize{    
#'    \item sample_cat. Sample locations used in categorical surface genetic simulation example
#'    \item sample_cont. Sample locations used in continuous surface genetic simulation example
#'    \item sample_multi. Sample locations used in multi surface genetic simulation example
#'    }
#' 
#' @docType data
#' @name sample_pops
#' @format A list of length 3
#' @usage data(sample_pops)
#' @description  Sample file to be used with examples in the vignette
#' @keywords datasets
NULL


#' Data file. A list of three matrices 
#' 
#' Each matrix depicts the pairwise genetic distance (measured as chord distance) between sample locations
#' 
#'  \itemize{    
#'    \item Dc_cat. Pairwise chord distance matrix simulated across the categorical resistance surface
#'    \item Dc_cont. Pairwise chord distance matrix simulated across the continuous resistance surface
#'    \item Dc_multi. Pairwise chord distance matrix simulated across the multivariate resistance surface
#'    }
#' 
#' @docType data
#' @name Dc_list
#' @format A list of length 3
#' @usage data(Dc_list)
#' @description  Sample file to be used with examples in the vignette
#' @keywords datasets
NULL


#' Data file. A list of three matrices 
#' 
#' Each matrix depicts the pairwise effective resistance distance (measured using `gdistance`) between sample locations across `resist_true` surfaces
#' 
#'  \itemize{      
#'    \item resist_cat. Pairwise effective distance matrix calculated across the categorical resistance surface
#'    \item resist_cont. Pairwise effective distance matrix calculated across the continuous resistance surface
#'    \item resist_cont. Pairwise effective distance matrix calculated across the multivariate resistance surface
#'    }
#' 
#' @docType data
#' @name resist_list
#' @format A list of length 3
#' @usage data(resist_list)
#' @description  Sample file to be used with examples in the vignette
#' @keywords datasets
NULL



#' Data file. Example sample location file
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

#' Data file. Simulated resistance surfaces
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