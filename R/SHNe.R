#' Calculate distance-based Spatial-Heterogeneity-in-Ne (SHNe) metrics
#'
#' Use formulas from Prunier et al. to calculate di and dhm
#'
#' @param n.samples Specify the number of samples (populations) included in your analysis
#' @param pop.size Provide a vector of values that indicates Ne for each sampled population. These values should be ordered the same as your sample populations
#' @param output Default is 'df', will return a data frame with distance measures. If 'matrix' is specified, then a list containing a matrix for each distance measure is produced.
#' @return Either a data frame or a list of matrices containing the SHNe metrics
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @references Prunier, J. G., V. Dubut, L. Chikhi, and S. Blanchet. 2017. Contribution of spatial heterogeneity in effective population sizes to the variance in pairwise measures of genetic differentiation. Methods in Ecology and Evolution 8:1866-1877.
#' 
#' #' @examples  
#' ## Create vector of Ne
#' Ne <- rep(c(5, 25, 100), each = 3)
#' SHNe.out <- SHNe(n.samples = length(Ne),
#'                  pop.size = Ne
#'                  )

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage SHNe(n.samples,
#'             pop.size,
#'             output = 'df')
#' 
SHNe <- function(n.samples, 
                 pop.size, 
                 output = 'df') {
  di <- dhm <- matrix(0, n.samples, n.samples)
  
  if(!is.vector(pop.size)) {
    pop.size <- pop.size[[1]]
  }
  
  for(i in 1:n.samples){
    for(j in 1:n.samples) {
      if (j <= i) {
        next
      } #else {
      dhm[j,i] <- -((2 * pop.size[i] * pop.size[j]) / (pop.size[i] + pop.size[j]))
      di[j,i] <- ((pop.size[i] + pop.size[j]) / (pop.size[i] * pop.size[j]))
      # }
    } # Close j
  } # Close i
  
  if(output == 'df') {
    out  <- data.frame(dhm = lower(dhm),
                       di = lower(di))
  } else {
    out <- list(dhm = dhm,
                di = di)
  }
  return(out)
} # End function