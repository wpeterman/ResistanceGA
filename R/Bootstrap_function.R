#' Run bootstrap on optimized resistance surfaces
#'
#' This is a 'pseudo' bootstrap procedure that will subsample a specified proportion of the sample locations/individuals, and will refit the MLPE model using the previously optimized resistance surface.
#'
#' @param mod.names A vector of the model names to be assessed
#' @param dist.mat A list containing all distance matrices from optimized resistance surfaces
#' @param n.parameters A vector the length of mod.names, specifying the number of parameters in each model
#' @param sample.prop Proportion of observations to be sampled each iteration (Default = 0.75)
#' @param iters Number of bootstrap iterations to be conducted
#' @param obs Total number of observations (populations or individuals) in your original analysis
#' @param rank.method What metric should be used to rank models during bootstrap analysis? c('AIC', 'AICc', 'R2', 'LL'). Default = 'AICc'
#' @param genetic.mat Genetic distance matrix without row or column names.
#' @return A data frame reporting the average model weight, average rank, number of time a model was the top model in the set, and the frequency a model was best.
#' 
#' @details This is a 'pseudo'bootstrap procedure that subsamples distance and genetic matrices, refits the MLPE model for each surface. AICc is calculated based on the number of parameters specified. Ranking of models during the bootstrap analysis is based on the specified \code{rank.method}, which defaults to 'AICc'. The objective of this procedure is to identify the surfaces that is top ranked across all bootstrap iterations.
#' 
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage Resist.boot(mod.names, 
#'                    dist.mat, 
#'                    n.parameters, 
#'                    sample.prop, 
#'                    iters, obs, 
#'                    rank.method, 
#'                    genetic.mat)
#'

Resist.boot <-
  function (mod.names,
            dist.mat,
            n.parameters,
            sample.prop = 0.75,
            iters,
            obs,
            rank.method = 'AICc',
            genetic.mat) {

    options(warn = -1)
    progress_bar <- plyr::progress_text()
    progress_bar$init(iters * length(mod.names))
    
    sample.n <- floor(sample.prop * obs)
    ID <- To.From.ID(sample.n)
    ZZ <- ZZ.mat(ID)
    
    k.mod <- data.frame(surface = mod.names, k = n.parameters)
    
    sample.list <-
      replicate(iters, sample(obs, size = sample.n, replace = F), simplify = F)
    AIC.tab.list <- vector(mode = "list", length = iters)
   
    names(dist.mat) <- mod.names
    
    for (i in 1:iters) {
      # bootsrap iteration
      samp <- sort(sample.list[[i]])
      genetic.samp <- lower(genetic.mat[samp, samp])
      AICc.tab <- vector(mode = "list", length = length(dist.mat))
      
      for (j in seq_along(dist.mat)) {
        # composite model
        dat <- lower(dist.mat[[j]][samp, samp])
        
        AICc <- boot.AICc(
          response = genetic.samp,
          resistance = dat,
          ID = ID,
          ZZ = ZZ,
          k = n.parameters[j],
          obs = length(samp)
        )
        
        mod.aic <- data.frame(mod.names[j], n.parameters[j], AICc)
        names(mod.aic) <- c("surface", "k", "AIC","AICc", "R2m", "LL")
        
        AICc.tab[[j]] <- mod.aic
        progress_bar$step()
        
      } # Close composite loop
      # Calculate Delta AICc, weight, and rank
      AICc.tab <- plyr::ldply(AICc.tab, "identity")
      
      AICc.tab <- AICc.tab %>% dplyr::mutate(., AIC = AIC) %>%
        dplyr::mutate(., AICc = AICc) %>%
        plyr::mutate(., delta = AICc - min(AICc)) %>%
        dplyr::mutate(., weight = (exp(-0.5 * delta)) / sum(exp(-0.5 * delta))) %>%
        dplyr::mutate(., iteration = i) %>%
        dplyr::mutate(., LL = LL) %>%
        dplyr::mutate(., R2m = R2m)
      
      if(rank.method == 'LL') {
        AICc.tab <- dplyr::mutate(AICc.tab, rank = dense_rank(desc(LL)))
      } else if (rank.method == 'AIC') {
        AICc.tab <- dplyr::mutate(AICc.tab, rank = dense_rank(AIC))
      } else if(rank.method == 'R2'){
        AICc.tab <- dplyr::mutate(AICc.tab, rank = dense_rank(desc(R2m)))
      } else {
        AICc.tab <- dplyr::mutate(AICc.tab, rank = dense_rank(AICc))
      }
      
      AIC.tab.list[[i]] <- AICc.tab
    } # Close iteration loop
    
    # Get average weight, rank, & R2
    group.list <- AIC.tab.list %>% plyr::ldply(.) %>% dplyr::group_by(., surface)
    boot.avg <-
      group.list %>% 
      dplyr::summarise(.,
                       avg.AIC = mean(AIC),
                       avg.AICc = mean(AICc),
                       avg.weight = mean(weight),
                       avg.rank = mean(rank),
                       avg.R2m = mean(R2m),
                       avg.LL = mean(LL)) %>% 
      arrange(., avg.rank)
    
    Freq_Percent <-
      group.list %>%  
      dplyr::filter(., rank == 1) %>% 
      dplyr::tally(.) %>%  
      dplyr::mutate(Percent.top = (100 * n) / sum(n))
    
    boot.avg <- dplyr::left_join(boot.avg, Freq_Percent, "surface") %>% 
      dplyr::left_join(., k.mod, "surface")
    boot.avg[is.na(boot.avg)] <- 0
    # boot.avg <- as.data.frame(boot.avg)
    
    return(boot.avg)
    
    options(warn = 0)
    
  } # End function

#!#!#!#!#!#!#!#!#!#!#!#!
# Other necessary functions
#!#!#!#!#!#!#!#!#!#!#!#!

# Make to-from population list
To.From.ID <- function(POPS) {
  tmp <- matrix(nrow = POPS, ncol = POPS)
  dimnames(tmp) <- list(1:POPS, 1:POPS)
  tmp2 <- as.data.frame(which(row(tmp) < col(tmp), arr.ind = TRUE))
  tmp2[[2]] <- dimnames(tmp)[[2]][tmp2$col]
  tmp2[[1]] <- dimnames(tmp)[[2]][tmp2$row]
  colnames(tmp2) <- c("pop1", "pop2")
  as.numeric(tmp2$pop1)
  as.numeric(tmp2$pop2)
  ID <- plyr::arrange(tmp2, as.numeric(pop1), as.numeric(pop2))
  #   ID<-tmp2[with(tmp2, order(pop1, pop2)), ]
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

# Bootstrap MLPE
boot.AICc <- function(response, resistance, ID, ZZ, k, obs) {
  resistance <- scale(resistance, center = TRUE, scale = TRUE)
  dat <- data.frame(ID, resistance = resistance, response = response)
  colnames(dat) <- c("pop1", "pop2", "resistance", "response")
  
  mod <- lFormula(response ~ resistance + (1 | pop1),
                  data = dat,
                  REML = F)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun, mod)
  opt <- optimizeLmer(dfun)
  fit.mod <- mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr)
  R.sq <- MuMIn::r.squaredGLMM(fit.mod)[[1]]
  mod.AIC <- AIC(fit.mod)
  LL <- logLik(fit.mod)
  AICc <- mod.AIC + ((2 * k * (k + 1)) / (obs - k - 1))
  out <- data.frame(AIC = mod.AIC, 
                    AICc = AICc, 
                    R2m = R.sq,
                    LL = LL)
  return(out)
}