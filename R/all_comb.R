#' Analyze all combinations
#'
#' A wrapper function to run both \code{\link[ResistanceGA]{SS_optim}} and \code{\link[ResistanceGA]{MS_optim}} to optimize all combinations of resistance surfaces with the Genetic Algorithms package \pkg{GA}. Following optimization, \code{\link[ResistanceGA]{Resist.boot}} is run to conduct a bootstrap analysis.  This function can only be used when optimizing resistance surface with least cost path or commute distance (\pkg{gdistance}).
#'
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function. Be sure that the \code{Results.dir} has been been correctly specified as "all.comb"
#' @param results.dir Directory to write and save analysis results. This should be an empty directory. Any existing files located in this directory will be deleted!
#' @param max.combination The maximum number of surfaces to include in the all combinations analysis (Default = 4). Alternatively, specify a vector with the minimum and maximum number of surfaces to combine (e.g., c(2,4). If the minimum > 1, then the single surface optimization will be skipped.
#' @param iters Number of bootstrap iterations to be conducted (Default = 1000)
#' @param sample.prop Proportion of observations to be sampled each iteration (Default = 0.75)
#' @param replicate The number of times to replicate the GA optimization process for each surface (Default = 1)
#' @param nlm Logical, if TRUE, the final step of optimization will use nlm to fine-tune parameter estimates. This may lead to overfitting in some cases. (Default = FALSE)
#' @param dist_mod Logical, if TRUE, a Distance model will be calculated and added to the output table (default = TRUE)
#' @param null_mod Logical, if TRUE, an intercept-only model will be calculated and added to the output table (Default = TRUE)

#' @return This function optimizes resistance surfaces in isolation using \code{\link[ResistanceGA]{SS_optim}}, followed by multisurface optimization using \code{\link[ResistanceGA]{MS_optim}}, and then conducts a bootstrap analysis.
#' @usage all_comb(gdist.inputs, 
#'                 GA.inputs, 
#'                 results.dir,
#'                 max.combination = 4,
#'                 iters = 1000,
#'                 replicate = 1,
#'                 sample.prop = 0.75,
#'                 nlm = FALSE,
#'                 dist_mod = TRUE,
#'                 null_mod = TRUE)
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @export
all_comb <- function(gdist.inputs,
                     GA.inputs,
                     results.dir,
                     max.combination = 4,
                     iters = 1000,
                     replicate = 1,
                     sample.prop = 0.75,
                     nlm = FALSE,
                     dist_mod = TRUE,
                     null_mod = TRUE) {
  
  if(!exists('results.dir')) 
    return(cat("ERROR: An empty results directory must be specified"))
  
  if(!exists('gdist.inputs')) 
    return(cat("ERROR: Please specify gdist.inputs"))
  
  if(!exists('GA.inputs')) 
    return(cat("ERROR: Please specify GA.inputs"))
  
  if(!is.null(GA.inputs$Results.dir) & 
     !is.null(GA.inputs$Write.dir) &
     !is.null(GA.inputs$Plots.dir)) {
    return(cat("ERROR: Please correctly specify the `Results.dir` as 'all.comb' when running GA.prep"))
  }
  
  if(length(max.combination) > 2) {
    return(cat("ERROR: Please specify either a single value or a vector with the minimum and maximum value"))
  }
  
  dir.files <- list.files(results.dir)
  if(length(dir.files) != 0) {
    q <- yn.question(cat(
      paste0("This function is about to delete all files and folders in '", results.dir, "'"),
      '\n', '\n',
      paste0("Do you want to proceed? Select 1 (Yes), 2 (No), or 3 (create a new subdirectory),  then press [Enter]")))
    
    # if(q == FALSE) return(cat("Function stopped"))
    
    if(is.na(q)) { # Create subdir
      dir.NAME <- floor(as.numeric(as.POSIXct(Sys.time())))
      dir.create(path = paste0(results.dir, "all_comb_", dir.NAME))
      results.dir <- paste0(results.dir, "all_comb_", dir.NAME, "/")
    } else if(q == FALSE) { # Stop function
      return(cat("Function stopped"))
    } else { # Remove exisiting folder
      unlink(dir(results.dir, 
                 full.names = TRUE),
             recursive = TRUE,
             force = T
      )
    }
  }
  
 
  
  # if(length(max.combination) == 2) {
  #   if(max.combination[2] > GA.inputs$n.layers) {
  #     return(cat("ERROR: Please specify a maximum combination that is less than or equal to the number of raster layers in the analysis"))
  #   }
  # }

  
  # Create combination list -------------------------------------------------
  mc <- max.combination
  
  if(length(max.combination) == 2) {
    if(mc[1] == 1) {
      min.combination <- 2
      max.combination <- mc[2]
      ss <- TRUE
    } else {
      min.combination <- mc[1]
      max.combination <- mc[2]
      ss <- FALSE
    } 
  } else {
    min.combination <- 2
    ss <- TRUE
  }
  
  if(max.combination > GA.inputs$n.layers) {
    max.combination <- GA.inputs$n.layers
  }
  comb.list <- vector(mode = "list", length = (max.combination - 1))
  
  
  list.count <- 0
  surface.count <- 0
  for(i in min.combination:max.combination) {
    list.count <- list.count + 1
    comb.list[[list.count]] <- t(combn(1:GA.inputs$n.layers, i))
    if(is.null(nrow(comb.list[[list.count]]))) {
      n.comb <- 1
    } else {
      n.comb <- nrow(comb.list[[list.count]])
    }
    surface.count <- surface.count + n.comb
  }
  
  all.combs <- list()
  comb.names <- list()
  row.index <- 0
  for(i in 1:length(comb.list)){
    combs <- comb.list[[i]]
    
    if(is.null(nrow(comb.list[[i]]))) {
      t.combs <- 1
    } else {
      t.combs <- nrow(comb.list[[i]])
    }
    
    for(j in 1:t.combs) {
      row.index <- row.index + 1
      all.combs[[row.index]] <- combs[j,]
      c.names <- GA.inputs$layer.names[combs[j,]]
      comb.names[[row.index]] <- paste(c.names, collapse = ".")
    }
  }
  
  GA.input_orig <- GA.inputs
  
  Results <- vector(mode = 'list', length = replicate)
  # Begin Replicate Loop --------------------------------------------------
  for(i in 1:replicate){
    # Skip if min combination > 1
    if(ss == FALSE) {
      ss.results <- NULL
      AICc.tab <- NULL
      dir.create(paste0(results.dir,'rep_',i))
    } else {  # Do single surface optimization
      dir.create(paste0(results.dir,'rep_',i))
      dir.create(paste0(results.dir,'rep_',i, "/", "single.surface"))
      dir.create(paste0(results.dir,'rep_',i, "/", "single.surface/", "Plots"))
      
      # Single Surface optimization -----------------------------------------------------
      if(!is.null(GA.inputs$scale)) {
        
        # * Single Surface: scaled --------------------------------------------------------
        
        # Update GA.input directories
        GA.inputs$Plots.dir <- paste0(results.dir,
                                      'rep_',i, 
                                      "/",
                                      "single.surface/",
                                      "Plots/")
        
        GA.inputs$Results.dir <- paste0(results.dir,
                                        'rep_',i, 
                                        "/", 
                                        "single.surface/")
        
        ss.results <- SS_optim.scale(gdist.inputs = gdist.inputs,
                                     GA.inputs = GA.inputs,
                                     nlm = nlm,
                                     dist_mod = dist_mod,
                                     null_mod = null_mod)
        
        AICc.tab <- ss.results$AICc 
      } else {
        # * Single Surface --------------------------------------------------------
        
        # Update GA.input directories
        GA.inputs$Plots.dir <- paste0(results.dir,
                                      'rep_',i, 
                                      "/",
                                      "single.surface/",
                                      "Plots/")
        
        GA.inputs$Results.dir <- paste0(results.dir,
                                        'rep_',i, 
                                        "/", 
                                        "single.surface/")
        
        ss.results <- SS_optim(gdist.inputs = gdist.inputs,
                               GA.inputs = GA.inputs,
                               nlm = nlm,
                               dist_mod = dist_mod,
                               null_mod = null_mod)
        
        AICc.tab <- ss.results$AICc
      }
    }
    
    
    # Multisurface optimization -----------------------------------------------
    
    # * Multisurface: scaled ----------------------------------------------------------
    
    if(!is.null(GA.inputs$scale)) {
      
      ms.cd <- vector(mode = 'list',
                      length = length(all.combs))
      
      ms.k <- vector(mode = 'list',
                     length = length(all.combs))
      
      AICc.tab_list <- vector(mode = 'list',
                              length = length(all.combs))
      
      ms.results <- vector(mode = "list", length = length(all.combs))
      
      if(is.null(ss.results)) {
        n_ss.cd <- 0
        all.cd <- ms.cd
      } else {
        n_ss.cd <- length(ss.results$cd)
        all.cd <- c(ss.results$cd,
                    ms.cd)
      }
      
      
      for(j in 1:length(all.combs)) {
        dir.create(paste0(results.dir,'rep_',i, "/", comb.names[[j]]))
        # dir.create(paste0(results.dir,'rep_',i, "/", comb.names[[j]],"/Plots"))
        
        # Select raster surfaces
        r.vec <- 1:GA.input_orig$n.layers
        drop.vec <- r.vec[!(r.vec %in% all.combs[[j]])]
        asc.comb <- dropLayer(GA.input_orig$Resistance.stack, drop.vec)
        
        if(!is.null(GA.input_orig$inputs$select.trans)) {
          s.trans <- GA.input_orig$inputs$select.trans[all.combs[[j]]]
        } else {
          s.trans <- GA.input_orig$inputs$select.trans
        }
        
        if(is.null(GA.input_orig$inputs$scale)) {
          ms.scale <- FALSE
        } else {
          ms.scale <- TRUE
        }
        
        if(sum(GA.input_orig$inputs$scale.surfaces) == 0) {
          sc.surf <- NULL
        } else {
          sc.surf <- GA.input_orig$inputs$scale.surfaces[all.combs[[j]]]
        }
        
        # sc.surf <- GA.input_orig$inputs$scale.surfaces[all.combs[[j]]]
        
        # Update GA.input 
        GA.inputs <- GA.prep(ASCII.dir = asc.comb,
                             Results.dir = 'all.comb',
                             min.cat = GA.input_orig$inputs$min.cat,
                             max.cat = GA.input_orig$inputs$max.cat,
                             max.cont = GA.input_orig$inputs$max.cont,
                             min.scale = GA.input_orig$inputs$min.scale,
                             max.scale = GA.input_orig$inputs$max.scale,
                             cont.shape = NULL,
                             select.trans = s.trans,
                             method = GA.input_orig$inputs$method,
                             scale = ms.scale,
                             scale.surfaces = sc.surf,
                             k.value = GA.input_orig$inputs$k.value,
                             pop.mult = GA.input_orig$inputs$pop.mult,
                             percent.elite = GA.input_orig$inputs$percent.elite,
                             type = GA.input_orig$inputs$type,
                             pcrossover = GA.input_orig$inputs$pcrossover,
                             pmutation = GA.input_orig$inputs$pmutation,
                             maxiter = GA.input_orig$inputs$maxiter,
                             run = GA.input_orig$inputs$run,
                             keepBest = GA.input_orig$inputs$keepBest,
                             population = GA.input_orig$inputs$population,
                             selection = GA.input_orig$inputs$selection,
                             crossover = GA.input_orig$inputs$crossover,
                             mutation = GA.input_orig$inputs$mutation,
                             pop.size = GA.input_orig$inputs$pop.size,
                             parallel = GA.input_orig$inputs$parallel,
                             seed = GA.input_orig$inputs$seed,
                             quiet = GA.input_orig$inputs$quiet
        )
        
        # Update GA.input directories
        GA.inputs$Plots.dir <- paste0(results.dir,
                                      'rep_',i, 
                                      "/", 
                                      comb.names[[j]],
                                      "/")
        
        GA.inputs$Results.dir <- paste0(results.dir,
                                        'rep_',i, 
                                        "/", 
                                        comb.names[[j]],
                                        "/")
        
        ms.results[[j]] <- MS_optim.scale(gdist.inputs = gdist.inputs,
                                          GA.inputs = GA.inputs)
        
        all.cd[[(j + n_ss.cd)]] <- ms.results[[j]]$cd[[1]]
        
        ms.k[[j]] <- ms.results[[j]]$k
        
        AICc.tab_list[[j]] <- ms.results[[j]]$AICc.tab
      } # End combination loop
      
      # Convert combination lists to data frames
      all.k <- rbind(ss.results$k,
                     plyr::ldply(ms.k))
      
      names(all.cd) <- all.k$surface
      
      if(is.null(ss.results)) {
        all.AICc <- plyr::ldply(AICc.tab_list)
      } else {
        all.AICc <- rbind(ss.results$AICc,
                          plyr::ldply(AICc.tab_list))  
      }
      
      
      all.AICc <- all.AICc %>% 
        dplyr::arrange(., AICc) %>%
        dplyr::mutate(., delta.AICc = AICc - min(AICc)) %>%
        dplyr::mutate(., weight = (exp(-0.5 * delta.AICc)) / sum(exp(-0.5 * delta.AICc))) %>%
        as.data.frame()
      
      
      # * Bootstrap: scaled -----------------------------------------------------
      
      obs <- gdist.inputs$n.Pops
      genetic.mat <- matrix(0, obs, obs)
      genetic.mat[lower.tri(genetic.mat)] <- gdist.inputs$response
      
      boot.results <- Resist.boot(mod.names = all.k[,1],
                                  dist.mat = all.cd,
                                  n.parameters = all.k[,2],
                                  sample.prop = sample.prop,
                                  iters = iters,
                                  obs = obs,
                                  genetic.mat = genetic.mat)      
    } else {     # End scaled optimization
      
      # Optimization without scaling --------------------------------------------
      
      # * Multisurface ----------------------------------------------------------
      
      ms.cd <- vector(mode = 'list',
                      length = length(all.combs))
      
      ms.k <- vector(mode = 'list',
                     length = length(all.combs))
      
      AICc.tab_list <- vector(mode = 'list',
                              length = length(all.combs))
      
      ms.results <- vector(mode = "list", length = length(all.combs))
      
      n_ss.cd <- length(ss.results$cd)
      all.cd <- c(ss.results$cd,
                  ms.cd)
      
      for(j in 1:length(all.combs)) {
        dir.create(paste0(results.dir,'rep_',i, "/", comb.names[[j]]))
        # dir.create(paste0(results.dir,'rep_',i, "/", comb.names[[j]],"/Plots"))
        
        # Select raster surfaces
        r.vec <- 1:GA.input_orig$n.layers
        drop.vec <- r.vec[!(r.vec %in% all.combs[[j]])]
        asc.comb <- dropLayer(GA.input_orig$Resistance.stack, drop.vec)
        
        if(!is.null(GA.input_orig$inputs$select.trans)) {
          s.trans <- GA.input_orig$inputs$select.trans[all.combs[[j]]]
        } else {
          s.trans <- GA.input_orig$inputs$select.trans
        }
        
        if(is.null(GA.input_orig$inputs$scale)) {
          ms.scale <- FALSE
        } else {
          ms.scale <- TRUE
        }
        
        if(sum(GA.input_orig$inputs$scale.surfaces) == 0) {
          sc.surf <- NULL
        } else {
          sc.surf <- GA.input_orig$inputs$scale.surfaces[all.combs[[j]]]
        }
        
        # sc.surf <- GA.input_orig$inputs$scale.surfaces[all.combs[[j]]]
        
        # Update GA.input 
        GA.inputs <- GA.prep(ASCII.dir = asc.comb,
                             Results.dir = 'all.comb',
                             min.cat = GA.input_orig$inputs$min.cat,
                             max.cat = GA.input_orig$inputs$max.cat,
                             max.cont = GA.input_orig$inputs$max.cont,
                             min.scale = GA.input_orig$inputs$min.scale,
                             max.scale = GA.input_orig$inputs$max.scale,
                             cont.shape = NULL,
                             select.trans = s.trans,
                             method = GA.input_orig$inputs$method,
                             scale = ms.scale,
                             scale.surfaces = sc.surf,
                             k.value = GA.input_orig$inputs$k.value,
                             pop.mult = GA.input_orig$inputs$pop.mult,
                             percent.elite = GA.input_orig$inputs$percent.elite,
                             type = GA.input_orig$inputs$type,
                             pcrossover = GA.input_orig$inputs$pcrossover,
                             pmutation = GA.input_orig$inputs$pmutation,
                             maxiter = GA.input_orig$inputs$maxiter,
                             run = GA.input_orig$inputs$run,
                             keepBest = GA.input_orig$inputs$keepBest,
                             population = GA.input_orig$inputs$population,
                             selection = GA.input_orig$inputs$selection,
                             crossover = GA.input_orig$inputs$crossover,
                             mutation = GA.input_orig$inputs$mutation,
                             pop.size = GA.input_orig$inputs$pop.size,
                             parallel = GA.input_orig$inputs$parallel,
                             seed = GA.input_orig$inputs$seed,
                             quiet = GA.input_orig$inputs$quiet
        )
        
        # Update GA.input directories
        GA.inputs$Plots.dir <- paste0(results.dir,
                                      'rep_',i, 
                                      "/", 
                                      comb.names[[j]],
                                      "/")
        
        GA.inputs$Results.dir <- paste0(results.dir,
                                        'rep_',i, 
                                        "/", 
                                        comb.names[[j]],
                                        "/")
        
        ms.results[[j]] <- MS_optim(gdist.inputs = gdist.inputs,
                                    GA.inputs = GA.inputs)
        
        all.cd[[(j + n_ss.cd)]] <- ms.results[[j]]$cd[[1]]
        
        ms.k[[j]] <- ms.results[[j]]$k
        
        AICc.tab_list[[j]] <- ms.results[[j]]$AICc.tab
      } # End combination loop
      
      # Convert combination lists to data frames
      all.k <- rbind(ss.results$k,
                     plyr::ldply(ms.k))
      
      names(all.cd) <- all.k$surface
      
      if(is.null(ss.results)) {
        all.AICc <- plyr::ldply(AICc.tab_list)
      } else {
        all.AICc <- rbind(ss.results$AICc,
                          plyr::ldply(AICc.tab_list))  
      }
      
      all.AICc <- all.AICc %>% 
        dplyr::arrange(., AICc) %>%
        dplyr::mutate(., delta.AICc = AICc - min(AICc)) %>%
        dplyr::mutate(., weight = (exp(-0.5 * delta.AICc)) / sum(exp(-0.5 * delta.AICc))) %>%
        as.data.frame()
      
      
      # * Bootstrap -----------------------------------------------------
      
      obs <- gdist.inputs$n.Pops
      genetic.mat <- matrix(0, obs, obs)
      genetic.mat[lower.tri(genetic.mat)] <- gdist.inputs$response
      
      boot.results <- Resist.boot(mod.names = all.k[,1],
                                  dist.mat = all.cd,
                                  n.parameters = all.k[,2],
                                  sample.prop = sample.prop,
                                  iters = iters,
                                  obs = obs,
                                  genetic.mat = genetic.mat)
      
    } # End scaled if-else
    
    # Write AICc table and Boot Results to replicate results directory
    write.table(x = all.AICc,
                paste0(results.dir,'rep_',i,"/",
                       "All_Combinations_Summary.csv"),
                row.names = F,
                col.names = T,
                sep = ",")
    
    write.table(x = as.data.frame(boot.results),
                paste0(results.dir,'rep_',i,"/",
                       "Bootstrap_Results.csv"),
                row.names = F,
                col.names = T,
                sep = ",")
    
    if(replicate > 1) {
      Results[[i]] <- list(summary.table = all.AICc,
                           boot.results = boot.results,
                           all.k = all.k,
                           all.cd = all.cd,
                           genetic.dist_mat = genetic.mat,
                           ss.results = ss.results,
                           ms.results = ms.results
                           )
      names(Results)[i] <- paste0('rep_',i)
      
    } else {
      Results <- list(summary.table = all.AICc,
                      boot.results = boot.results,
                      all.k = all.k,
                      all.cd = all.cd,
                      genetic.dist_mat = genetic.mat,
                      ss.results = ss.results,
                      ms.results = ms.results
                      )
    }
    
  } # Close replicate loop
  
  return(Results)
  
} # End function
