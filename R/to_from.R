# Make to-from population list
#' Convenience function to make to-from object needed for specifying MLPE random effects
#'
#' @param sampled_pops Provide an integer value representing the number of individuals or populations that were sampled
#' @param pop_n Provide a vector of length equal to `sampled_pops`, indicating the number of individuals sampled at each population. Only needed if going from a population- to an individual-based analysis
#' @param spLoc \code{\link[sp]{SpatialPoints}} object of sample locations
#' @param nb Distance between points that should be considered as part of the same group or neighborhood. If \code{spLoc} is specified, \code{nb} must also be specified. See Details.
#' @param ... Additional arguments. Currently only `dist` argument is accepted.
#' @return A data frame with two columns 
#' @details This function creates an object indicating the population pairs that distance values were calculated between. Note: Distance values must be taken from the lower half of a distance matrix. You can use \code{\link[ResistanceGA]{lower}} to obtain these values. When specifying a MLPE model using \code{\link[ResistanceGA]{mlpe_rga}}, \code{pop1} from \code{To.From.ID} must be added to the data frame containing model data, and must be specified as the random effect.
#' 
#' By providing \code{spLoc}, you can add an additional random effect to account for spatial proximity of sample locations. You must also specify \code{nb}, which indicates the maximum distance (in units of your \code{spLoc} object) that locations are considered to be part of the same neighborhood.
#' 
#' @examples 
#' ## Standard use
#' To.From.ID(sampled_pops = 4)
#' 
#' ## If going from population to individual-level 
#' To.From.ID(sampled_pops = 4,
#'            pop_n = c(1,2,3,1))
#'          
#' ## A vector of pairwise distances can also be expanded            
#' To.From.ID(sampled_pops = 4,
#'            pop_n = c(1,2,3,1),
#'            dist = runif(6)) 
#'            
#' ## Spatial location incorporated
#' spLoc <- SpatialPoints(cbind(c(1,2,9,10), 
#'                              c(1,2,9,10)))
#'           
#' To.From.ID(sampled_pops = 4,
#'            spLoc = spLoc,
#'            nb = 3)              

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage To.From.ID(sampled_pops,
#'                   pop_n = NULL,
#'                   spLoc = NULL,
#'                   nb = NULL,
#'                   ...)
#' 

To.From.ID <- function(sampled_pops, 
                       pop_n = NULL,
                       spLoc = NULL,
                       nb = NULL,
                       ...) { 
  
  args <- as.list(substitute(list(...)))[-1L]
  
  if(missing("sampled_pops")) { 
    if(any("POPS" %in% names(args))) {
      sampled_pops <- eval(args[['POPS']])
    }
    warning("Argument `POPS` deprecated. Use `sampled_pops`")
  }
  
  
  if(!is.null(spLoc) & is.null(nb)) {
    stop("You must specify `nb` when providing 'spLoc'")
  }
  
  if(!is.null(pop_n)) {
    tmp <- matrix(0, nrow = sampled_pops, ncol = sampled_pops)
    colmat <- col(tmp)
    rowmat <- row(tmp)
    
    colmat[sampled_pops, 1] <- sampled_pops
    rowmat[sampled_pops, 1] <- 1
    
    colmat <- colmat[rep(1:nrow(colmat), times = pop_n), 
                     rep(1:ncol(colmat), times = pop_n)]
    
    rowmat <- rowmat[rep(1:nrow(rowmat), times = pop_n), 
                     rep(1:ncol(rowmat), times = pop_n)]
    
    ## Add correlation group
    if(!is.null(spLoc)) {
      nb_ <- spdep::dnearneigh(spLoc, d1 = 0, d2 = nb)
      nb.mat <- spdep::nb2mat(nb_, style = "B", zero.policy = TRUE)
      nb.mat_ind <- nb.mat[rep(1:nrow(nb.mat), times = pop_n),
                           rep(1:ncol(nb.mat), times = pop_n)]
      

      pop.ID <- data.frame(pop1 = factor(lower(colmat)),
                           pop2 = factor(lower(rowmat)),
                           corr_ = lower(nb.mat_ind),
                           cor.grp = factor(lower(colmat) * lower(nb.mat_ind)))
    } else {
      pop.ID <- data.frame(pop1 = factor(lower(colmat)),
                           pop2 = factor(lower(rowmat)))
    }
    
    pop <- rep(1:sum(pop_n), 
               length.out = length(pop.ID$pop2))
    
    keep <- ifelse(apply(pop.ID, 1, function(x) x[1] == x[2]), 0, 1)
    
    grp <- paste0(pop.ID$pop1, "_", pop.ID$pop2)
    
    if(any("dist" %in% names(args))){
      dist.mat <- tmp
      dist.mat[lower.tri(dist.mat)] <- eval(args[['dist']])
      
      dist.mat <- dist.mat[rep(1:nrow(dist.mat), times = pop_n), 
                           rep(1:ncol(dist.mat), times = pop_n)]
      
      dist <- lower(dist.mat)
      
      ## Individual to-from  
      tmp <- matrix(0, nrow = sampled_pops, ncol = sampled_pops)
      colmat <- col(tmp)
      rowmat <- row(tmp)
      
      colmat <- colmat[rep(1:nrow(colmat), times = pop_n), 
                       rep(1:ncol(colmat), times = pop_n)]
      
      rowmat <- rowmat[rep(1:nrow(rowmat), times = pop_n), 
                       rep(1:ncol(rowmat), times = pop_n)]
      
      ind.colmat <- col(rowmat)
      ind.rowmat <- row(rowmat)
      
      ind.colmat[sum(pop_n), 1] <- sum(pop_n)
      ind.rowmat[sum(pop_n), 1] <- 1
      
      ind.mat <- col(rowmat)
      
      ind <- lower(ind.mat)
      
      if(!is.null(spLoc)) {
        
        ind.ID <- data.frame(pop1.ind = factor(lower(ind.colmat)),
                             pop2.ind = factor(lower(ind.rowmat)),
                             pop1.pop = pop.ID$pop1,
                             pop2.pop = pop.ID$pop2,
                             cor.grp = pop.ID$cor.grp,
                             corr_ = pop.ID$corr_,
                             dist = dist,
                             pop = pop)
        
      } else {
        ind.ID <- data.frame(pop1.ind = factor(lower(ind.colmat)),
                             pop2.ind = factor(lower(ind.rowmat)),
                             pop1.pop = pop.ID$pop1,
                             pop2.pop = pop.ID$pop2,
                             dist = dist,
                             pop = pop)
      }
      
    } else {
      ## Individual to-from         
      tmp <- matrix(0, nrow = sampled_pops, ncol = sampled_pops)
      colmat <- col(tmp)
      rowmat <- row(tmp)
      
      colmat <- colmat[rep(1:nrow(colmat), times = pop_n), 
                       rep(1:ncol(colmat), times = pop_n)]
      
      rowmat <- rowmat[rep(1:nrow(rowmat), times = pop_n), 
                       rep(1:ncol(rowmat), times = pop_n)]
      
      ind.colmat <- col(rowmat)
      ind.rowmat <- row(rowmat)
      
      ind.colmat[sum(pop_n), 1] <- sum(pop_n)
      ind.rowmat[sum(pop_n), 1] <- 1
      
      if(!is.null(spLoc)) {
        ind.ID <- data.frame(pop1.ind = factor(lower(ind.colmat)),
                             pop2.ind = factor(lower(ind.rowmat)),
                             pop1.pop = pop.ID$pop1,
                             pop2.pop = pop.ID$pop2,
                             cor.grp = pop.ID$cor.grp,
                             corr_ = pop.ID$corr_,
                             pop.df = pop,
                             cor.df = cor.df)
      } else {
        ind.ID <- data.frame(pop1.ind = factor(lower(ind.colmat)),
                             pop2.ind = factor(lower(ind.rowmat)),
                             pop1.pop = pop.ID$pop1,
                             pop2.pop = pop.ID$pop2,
                             pop = pop)
      }
      
      
    }
    
    ind.ID <- ind.ID[keep == 1,]   
    return(ind.ID)
    
  } else { # End expanded to-from
    tmp <- matrix(nrow = sampled_pops, ncol = sampled_pops)
    colmat <- col(tmp)
    rowmat <- row(tmp)
    
    colmat[sampled_pops, 1] <- sampled_pops
    
    rowmat[sampled_pops, 1] <- 1
    
    ## Add correlation group
    if(!is.null(spLoc)) {
      nb_ <- spdep::dnearneigh(spLoc, d1 = 0, d2 = nb)
      nb.mat <- spdep::nb2mat(nb_, style = "B", zero.policy = TRUE)
      l_nb.mat <- lower(nb.mat)
      
      ID <- data.frame(pop1 = factor(lower(colmat)),
                       pop2 = factor(lower(rowmat)),
                       corr_ = lower(nb.mat),
                       cor.grp = factor(lower(colmat)))
    } else {
      ID <- data.frame(pop1 = factor(lower(colmat)),
                       pop2 = factor(lower(rowmat)))
    }
    return(ID)
  } # End standard to-from
} # End Function


# Original function -------------------------------------------------------

To.From.ID_ <- function(POPS) {
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