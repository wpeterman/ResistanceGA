#' Make a vector of the lower half of a square distance matrix
#'
#' This function will prepare and compile objects and commands for optimization with the GA package
#'
#' @param matrix Square distance matrix with no row names.
#' @return A vector of the lower half of the matrix
#'
#' @details This is a convenience function to obtain the lower half of a matrix, which is required as input for several other functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>

lower <- function(matrix) {
  if (is.vector(matrix) == TRUE ||
      dim(matrix)[1] != dim(matrix)[2]) {
    warning("Must provide square distance matrix with no column or row names")
  }
  lm <- matrix[lower.tri(matrix)]
  return(lm)
}
