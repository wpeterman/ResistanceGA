# Make to-from population list
#' Convenience function to make to-from object needed for specifying MLPE random effects
#'
#' @param POPS Provide an integer value representing the number of individuals or populations that were sampled
#' @return A data frame with two columns 
#' @details This object indicates the population pairs that distance values were calculated between. Note: Distance values must be taken from the lower half of a distance matrix. You can use \code{\link[ResistanceGA]{lower}} to obtain these values. When specifying a MLPE model using \code{\link[ResistanceGA]{mlpe_rga}}, \code{pop1} from \code{To.From.ID} must be added to the data frame containing model data, and must be specified as the random effect.
#' 
#' @examples To.From.ID(3)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage To.From.ID(POPS)

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