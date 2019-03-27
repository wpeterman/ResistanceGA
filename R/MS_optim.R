#' Simultaneous optimization of multiple resistance surfaces
#'
#' Optimize multiple resistance surfsaces simultaneously using genetic algorithms
#'
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param jl.inputs Object created from running \code{\link[ResistanceGA]{jl.prep}} function. Defined if optimizing using CIRCUITSCAPE run in Julia
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @return This function optimizes multiple resistance surfaces, returning a Genetic Algorithm (GA) object with summary information. Diagnostic plots of model fit are output to the "Results/Plots" folder that is automatically generated within the folder containing the optimized ASCII files. A text summary of the optimization settings and results is printed to the results folder.
#' @usage MS_optim(CS.inputs, 
#'              gdist.inputs, 
#'              jl.inputs,
#'              GA.inputs)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
MS_optim <- function(CS.inputs = NULL,
                     gdist.inputs = NULL,
                     jl.inputs = NULL,
                     GA.inputs) {
  if (!is.null(GA.inputs$scale)) {
    stop(
      "This function should NOT be used if you intend to apply kernel smoothing to your resistance surfaces"
    )
  }
  k.value <- GA.inputs$k.value
  
  
  # Circuitscape ------------------------------------------------------------
  
  
  if (!is.null(CS.inputs)) {
    # if (GA.inputs$parallel != FALSE) {
    #   warning(
    #     "\n CIRCUITSCAPE cannot be optimized in parallel. \n Ignoring parallel arguement. \n If you want to optimize in parallel, use least cost paths or commute time with gdistance.",
    #     immediate. = TRUE
    #   )
    # }
    t1 <- proc.time()[3]
    
    multi.GA_nG <- ga(
      type = "real-valued",
      fitness = Resistance.Opt_multi,
      population = GA.inputs$population,
      selection = GA.inputs$selection,
      mutation = GA.inputs$mutation,
      pcrossover = GA.inputs$pcrossover,
      crossover = GA.inputs$crossover,
      pmutation = GA.inputs$pmutation,
      Min.Max = GA.inputs$Min.Max,
      GA.inputs = GA.inputs,
      CS.inputs = CS.inputs,
      lower = GA.inputs$ga.min,
      upper = GA.inputs$ga.max,
      popSize = GA.inputs$pop.size,
      maxiter = GA.inputs$maxiter,
      optim = GA.inputs$optim,
      optimArgs = GA.inputs$optimArgs,
      parallel = GA.inputs$parallel,
      run = GA.inputs$run,
      keepBest = GA.inputs$keepBest,
      seed = GA.inputs$seed,
      suggestions = GA.inputs$SUGGESTS,
      monitor = GA.inputs$monitor,
      quiet = GA.inputs$quiet
    )
    rt <- proc.time()[3] - t1
    
    if(dim(multi.GA_nG@solution)[1] > 1) {
      multi.GA_nG@solution <- t(as.matrix(multi.GA_nG@solution[1,]))
    }
    
    Opt.parm <- GA.opt <- multi.GA_nG@solution
    for (i in 1:GA.inputs$n.layers) {
      if (GA.inputs$surface.type[i] == "cat") {
        ga.p <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        parm <- ga.p / min(ga.p)
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
        
      } else {
        parm <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
      }
    }
    multi.GA_nG@solution <- Opt.parm
    
    
    RAST <-
      Combine_Surfaces(
        PARM = multi.GA_nG@solution,
        CS.inputs = CS.inputs,
        GA.inputs = GA.inputs,
        rescale = TRUE,
        p.contribution = TRUE
      )
    
    p.cont <- RAST$percent.contribution
    RAST <- RAST$combined.surface
    
    NAME <- paste(GA.inputs$parm.type$name, collapse = ".")
    names(RAST) <- NAME
    
    cd <- Run_CS(CS.inputs,
                 r = RAST,
                 full.mat = TRUE,
                 EXPORT.dir = GA.inputs$Results.dir)
    
    write.table(
      cd,
      file = paste0(GA.inputs$Results.dir, NAME, "_csResistMat.csv"),
      sep = ",",
      row.names = F,
      col.names = F
    )
    writeRaster(RAST,
                paste0(GA.inputs$Results.dir, NAME, ".asc"),
                overwrite = TRUE)
    
    ifelse(length(unique(RAST)) > 15,
           type <- "continuous",
           type <- "categorical")
    
    Diagnostic.Plots(
      resistance.mat = lower(cd),
      genetic.dist = CS.inputs$response,
      plot.dir = GA.inputs$Plots.dir,
      type = "continuous",
      name = NAME,
      ID = CS.inputs$ID,
      ZZ = CS.inputs$ZZ
    )
    
    dat <- CS.inputs$df
    dat$cd <- scale(lower(cd))
    fit.mod <- mlpe_rga(formula = CS.inputs$formula,
                        data = dat,
                        ZZ = CS.inputs$ZZ,
                        REML = FALSE)
    fit.mod_REML <- mlpe_rga(formula = CS.inputs$formula,
                             data = dat,
                             ZZ = CS.inputs$ZZ,
                             REML = TRUE)
    
    aic <- suppressWarnings(AIC(
      fit.mod
    ))
    
    fit.stats <- suppressWarnings(r.squaredGLMM(
      fit.mod
    ))
    
    LL <- logLik(
      fit.mod
    )[[1]]
    
    MLPE.model <- fit.mod
    # fit.stats <- r.squaredGLMM(
    #   MLPE.lmm(
    #     REML = F,
    #     resistance = lower(cd),
    #     pairwise.genetic = CS.inputs$response,
    #     ID = CS.inputs$ID,
    #     ZZ = CS.inputs$ZZ
    #   )
    # )
    # aic <- AIC(
    #   MLPE.lmm(
    #     REML = F,
    #     resistance = lower(cd),
    #     pairwise.genetic = CS.inputs$response,
    #     ID = CS.inputs$ID,
    #     ZZ = CS.inputs$ZZ
    #   )
    # )
    # 
    # LL <-
    #   logLik(
    #     MLPE.lmm(
    #       resistance = lower(cd),
    #       pairwise.genetic = CS.inputs$response,
    #       REML = F,
    #       ID = CS.inputs$ID,
    #       ZZ = CS.inputs$ZZ
    #     )
    #   )
    # 
    # MLPE.model <-
    #   MLPE.lmm(
    #     resistance = lower(cd),
    #     pairwise.genetic = CS.inputs$response,
    #     REML = F,
    #     ID = CS.inputs$ID,
    #     ZZ = CS.inputs$ZZ
    #   )
    
    if (k.value == 1) {
      k <- 2
    } else if (k.value == 2) {
      k <- GA.inputs$parm.type$n.parm[i] + 
        length(lme4::fixef(fit.mod)) - 1
      
    } else if (k.value == 3) {
      k <- GA.inputs$parm.type$n.parm[i] + 
        length(GA.inputs$layer.names) + 
        length(lme4::fixef(fit.mod)) - 1
      
    } else {
      k <- length(GA.inputs$layer.names[i]) + 
        length(lme4::fixef(fit.mod)) - 1
      
    }
    
    n <- CS.inputs$n.Pops
    AICc <- (-2 * LL) + (2 * k) + ((2 * k) * (k + 1)) / (n - k - 1)
    
    
    # Get parameter estimates
    MLPE.results <- MLPE.lmm_coef(
      resistance = GA.inputs$Results.dir,
      genetic.dist = CS.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "cs",
      ID = CS.inputs$ID,
      ZZ = CS.inputs$ZZ
    )
    
    Result.txt(
      GA.results = multi.GA_nG,
      GA.inputs = GA.inputs,
      method = "CIRCUITSCAPE",
      Run.Time = rt,
      fit.stats = fit.stats,
      optim = GA.inputs$method,
      k = k,
      aic = aic,
      AICc = AICc,
      LL = LL[[1]]
    )
    
    write.table(p.cont, file = paste0(GA.inputs$Results.dir, "Percent_Contribution.csv"), sep = ",",
                row.names = F,
                col.names = T)
    
    # save(multi.GA_nG, 
    #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
    
    # file.remove(list.files(GA.inputs$Write.dir, full.names = TRUE))
    
    # unlink(GA.inputs$Write.dir, recursive = T, force = T)
    
    k.df <- data.frame(surface = NAME, k = k)
    
    cd.list <- list(as.matrix(cd))
    names(cd.list) <- NAME
    
    AICc.tab <- data.frame(surface = NAME,
                           obj = multi.GA_nG@fitnessValue,
                           k = k,
                           AIC = aic,
                           AICc = AICc,
                           R2m = fit.stats[[1]],
                           R2c = fit.stats[[2]],
                           LL = LL)
    
    colnames(AICc.tab) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        "k",
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL"
      )
    
    out <- list(GA.summary = multi.GA_nG,
                MLPE.model = MLPE.model,
                AICc.tab = AICc.tab,
                cd = cd.list,
                percent.contribution = p.cont,
                k = k.df)
    
    return(out)
  }
  
  
  # gdistance ---------------------------------------------------------------
  
  if (!is.null(gdist.inputs)) {
    
    
    # MLPE with Covariates ----------------------------------------------------
    
    if(!is.null(gdist.inputs$covariates)) { 
      #  * Island GA -------------------------------------------------------------
      if(isTRUE(GA.inputs$gaisl)) {
        t1 <- proc.time()[3]
        multi.GA_nG <- gaisl(
          type = "real-valued",
          fitness = Resistance.Opt_multi.cov,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.inputs,
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          numIslands = GA.inputs$numIslands,
          migrationRate = GA.inputs$migrationRate,
          migrationInterval = GA.inputs$migrationInterval,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          # keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      } else {
        # * Standard GA -------------------------------------------------------------
        
        t1 <- proc.time()[3]
        multi.GA_nG <- ga(
          type = "real-valued",
          fitness = Resistance.Opt_multi.cov,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.inputs,
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      }
      
    } # End covariates
    
    # MLPE no Covariates ------------------------------------------------------
    if(is.null(gdist.inputs$covariates)) {
      #  * Island GA -------------------------------------------------------------
      if(isTRUE(GA.inputs$gaisl)) {
        t1 <- proc.time()[3]
        multi.GA_nG <- gaisl(
          type = "real-valued",
          fitness = Resistance.Opt_multi,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.inputs,
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          numIslands = GA.inputs$numIslands,
          migrationRate = GA.inputs$migrationRate,
          migrationInterval = GA.inputs$migrationInterval,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          # keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      } else {
        # * Standard GA -------------------------------------------------------------
        
        t1 <- proc.time()[3]
        multi.GA_nG <- ga(
          type = "real-valued",
          fitness = Resistance.Opt_multi,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          mutation = GA.inputs$mutation,
          pcrossover = GA.inputs$pcrossover,
          crossover = GA.inputs$crossover,
          pmutation = GA.inputs$pmutation,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.inputs,
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$ga.min,
          upper = GA.inputs$ga.max,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          suggestions = GA.inputs$SUGGESTS,
          quiet = GA.inputs$quiet
        )
        rt <- proc.time()[3] - t1
        
      }
    } # End no covariates
    
    
    
    multi.GA_nG.o <- multi.GA_nG
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, paste(GA.inputs$parm.type$name, collapse = "."), "_full.rds"))
    
    if(dim(multi.GA_nG@solution)[1] > 1) {
      multi.GA_nG@solution <- t(as.matrix(multi.GA_nG@solution[1,]))
    }
    
    Opt.parm <- GA.opt <- multi.GA_nG@solution
    
    for (i in 1:GA.inputs$n.layers) {
      if (GA.inputs$surface.type[i] == "cat") {
        ga.p <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        parm <- ga.p / min(ga.p)
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
        
      } else {
        parm <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
      }
    }
    multi.GA_nG@solution <- Opt.parm
    
    RAST <-
      Combine_Surfaces(
        PARM = multi.GA_nG@solution,
        gdist.inputs = gdist.inputs,
        GA.inputs = GA.inputs,
        rescale = TRUE,
        p.contribution = TRUE
      )
    
    p.cont <- RAST$percent.contribution
    RAST <- RAST$combined.surface
    
    NAME <- paste(GA.inputs$parm.type$name, collapse = ".")
    names(RAST) <- NAME
    cd <- Run_gdistance(gdist.inputs, RAST)
    dat <- gdist.inputs$df
    dat$cd <- scale(c(cd))
    
    write.table(
      as.matrix(cd),
      file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,"_distMat.csv"),
      sep = ",",
      row.names = F,
      col.names = F
    )
    writeRaster(RAST,
                paste0(GA.inputs$Results.dir, NAME, ".asc"),
                overwrite = TRUE)
    
    ifelse(length(unique(RAST)) > 15,
           type <- "continuous",
           type <- "categorical")
    
    Diagnostic.Plots(
      resistance.mat = cd,
      genetic.dist = gdist.inputs$response,
      plot.dir = GA.inputs$Plots.dir,
      type = type,
      name = NAME,
      ID = gdist.inputs$ID,
      ZZ = gdist.inputs$ZZ
    )
    
    if(!is.null(gdist.inputs$covariates)) { 
      MLPE.results <- NULL
    } else {
      # Get parameter estimates
      MLPE.results <- MLPE.lmm_coef(
        resistance = GA.inputs$Results.dir,
        genetic.dist = gdist.inputs$response,
        out.dir = GA.inputs$Results.dir,
        method = "gd",
        ID = gdist.inputs$ID,
        ZZ = gdist.inputs$ZZ
      )
    }
    
    fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                        data = dat,
                        ZZ = gdist.inputs$ZZ,
                        REML = FALSE)
    fit.mod_REML <- mlpe_rga(formula = gdist.inputs$formula,
                             data = dat,
                             ZZ = gdist.inputs$ZZ,
                             REML = TRUE)
    
    aic <- suppressWarnings(AIC(
      fit.mod
    ))
    
    fit.stats <- suppressWarnings(r.squaredGLMM(
      fit.mod
    ))
    
    LL <- logLik(
      fit.mod
    )[[1]]
    
    MLPE.model <- fit.mod
    
    # MLPE.model <- MLPE.lmm(
    #   resistance = cd,
    #   pairwise.genetic = gdist.inputs$response,
    #   REML = F,
    #   ID = gdist.inputs$ID,
    #   ZZ = gdist.inputs$ZZ
    # )
    
    if (k.value == 1) {
      k <- 2
    } else if (k.value == 2) {
      k <- GA.inputs$parm.type$n.parm[i] + 
        length(lme4::fixef(fit.mod)) - 1
      
    } else if (k.value == 3) {
      k <- GA.inputs$parm.type$n.parm[i] + 
        length(GA.inputs$layer.names) + 
        length(lme4::fixef(fit.mod)) - 1
      
    } else {
      k <- length(GA.inputs$layer.names[i]) + 
        length(lme4::fixef(fit.mod)) - 1
      
    }
    
    n <- gdist.inputs$n.Pops
    AICc <- (-2 * LL) + (2 * k) + ((2 * k) * (k + 1)) / (n - k - 1)
    
    Result.txt(
      GA.results = multi.GA_nG.o,
      GA.inputs = GA.inputs,
      method = gdist.inputs$method,
      Run.Time = rt,
      fit.stats = fit.stats,
      optim = GA.inputs$method,
      k = k,
      aic = aic,
      AICc = AICc,
      LL = LL[[1]],
      fit.mod_REML = fit.mod_REML
    )
    
    write.table(p.cont, file = paste0(GA.inputs$Results.dir, "Percent_Contribution.csv"), sep = ",",
                row.names = F,
                col.names = T)
    
    # save(multi.GA_nG, 
    #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
    
    # unlink(GA.inputs$Write.dir, recursive = T, force = T)
    
    k.df <- data.frame(surface = NAME, k = k)
    
    cd.list <- list(as.matrix(cd))
    names(cd.list) <- NAME
    
    AICc.tab <- data.frame(surface = NAME,
                           obj = multi.GA_nG@fitnessValue,
                           k = k,
                           AIC = aic,
                           AICc = AICc,
                           R2m = fit.stats[[1]],
                           R2c = fit.stats[[2]],
                           LL = LL)
    
    colnames(AICc.tab) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        "k",
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL"
      )
    
    out <- list(GA.summary = multi.GA_nG.o,
                MLPE.model = MLPE.model,
                MLPE.model_REML = fit.mod_REML,
                AICc.tab = AICc.tab,
                cd = cd.list,
                percent.contribution = p.cont,
                k = k.df)
    
    return(out)
  }
  
  
  # Julia ---------------------------------------------------------------
  # MLPE with Covariates ----------------------------------------------------
  
  if(!is.null(gdist.inputs$covariates)) { 
    #  * Island GA -------------------------------------------------------------
    if(isTRUE(GA.inputs$gaisl)) {
      # stop("Optimization with covariates is not currently supported with gaisl!")
      
      t1 <- proc.time()[3]
      multi.GA_nG <- gaisl(
        type = "real-valued",
        fitness = Resistance.Opt_multi.cov,
        population = GA.inputs$population,
        selection = GA.inputs$selection,
        mutation = GA.inputs$mutation,
        pcrossover = GA.inputs$pcrossover,
        crossover = GA.inputs$crossover,
        pmutation = GA.inputs$pmutation,
        Min.Max = GA.inputs$Min.Max,
        GA.inputs = GA.inputs,
        jl.inputs = jl.inputs,
        lower = GA.inputs$ga.min,
        upper = GA.inputs$ga.max,
        popSize = GA.inputs$pop.size,
        maxiter = GA.inputs$maxiter,
        numIslands = GA.inputs$numIslands,
        migrationRate = GA.inputs$migrationRate,
        migrationInterval = GA.inputs$migrationInterval,
        optim = GA.inputs$optim,
        optimArgs = GA.inputs$optimArgs,
        parallel = GA.inputs$parallel,
        run = GA.inputs$run,
        # keepBest = GA.inputs$keepBest,
        seed = GA.inputs$seed,
        monitor = GA.inputs$monitor,
        suggestions = GA.inputs$SUGGESTS,
        quiet = GA.inputs$quiet
      )
      rt <- proc.time()[3] - t1
      
    } else {
      # * Standard GA -------------------------------------------------------------
      
      t1 <- proc.time()[3]
      multi.GA_nG <- ga(
        type = "real-valued",
        fitness = Resistance.Opt_multi.cov,
        population = GA.inputs$population,
        selection = GA.inputs$selection,
        mutation = GA.inputs$mutation,
        pcrossover = GA.inputs$pcrossover,
        crossover = GA.inputs$crossover,
        pmutation = GA.inputs$pmutation,
        Min.Max = GA.inputs$Min.Max,
        GA.inputs = GA.inputs,
        jl.inputs = jl.inputs,
        lower = GA.inputs$ga.min,
        upper = GA.inputs$ga.max,
        popSize = GA.inputs$pop.size,
        maxiter = GA.inputs$maxiter,
        optim = GA.inputs$optim,
        optimArgs = GA.inputs$optimArgs,
        parallel = GA.inputs$parallel,
        run = GA.inputs$run,
        keepBest = GA.inputs$keepBest,
        seed = GA.inputs$seed,
        monitor = GA.inputs$monitor,
        suggestions = GA.inputs$SUGGESTS,
        quiet = GA.inputs$quiet
      )
      rt <- proc.time()[3] - t1
      
    }
    
  } # End covariates

  # MLPE no Covariates ------------------------------------------------------
  
  if (!is.null(jl.inputs)) {
    #  * Island GA -------------------------------------------------------------
    if(isTRUE(GA.inputs$gaisl)) {
      t1 <- proc.time()[3]
      multi.GA_nG <- gaisl(
        type = "real-valued",
        fitness = Resistance.Opt_multi,
        population = GA.inputs$population,
        selection = GA.inputs$selection,
        mutation = GA.inputs$mutation,
        pcrossover = GA.inputs$pcrossover,
        crossover = GA.inputs$crossover,
        pmutation = GA.inputs$pmutation,
        Min.Max = GA.inputs$Min.Max,
        GA.inputs = GA.inputs,
        jl.inputs = jl.inputs,
        lower = GA.inputs$ga.min,
        upper = GA.inputs$ga.max,
        popSize = GA.inputs$pop.size,
        maxiter = GA.inputs$maxiter,
        numIslands = GA.inputs$numIslands,
        migrationRate = GA.inputs$migrationRate,
        migrationInterval = GA.inputs$migrationInterval,
        optim = GA.inputs$optim,
        optimArgs = GA.inputs$optimArgs,
        parallel = GA.inputs$parallel,
        run = GA.inputs$run,
        # keepBest = GA.inputs$keepBest,
        seed = GA.inputs$seed,
        monitor = GA.inputs$monitor,
        suggestions = GA.inputs$SUGGESTS,
        quiet = GA.inputs$quiet
      )
      rt <- proc.time()[3] - t1
      
    } else {
      # * Standard GA -------------------------------------------------------------
      
      t1 <- proc.time()[3]
      multi.GA_nG <- ga(
        type = "real-valued",
        fitness = Resistance.Opt_multi,
        population = GA.inputs$population,
        selection = GA.inputs$selection,
        mutation = GA.inputs$mutation,
        pcrossover = GA.inputs$pcrossover,
        crossover = GA.inputs$crossover,
        pmutation = GA.inputs$pmutation,
        Min.Max = GA.inputs$Min.Max,
        GA.inputs = GA.inputs,
        jl.inputs = jl.inputs,
        lower = GA.inputs$ga.min,
        upper = GA.inputs$ga.max,
        popSize = GA.inputs$pop.size,
        maxiter = GA.inputs$maxiter,
        optim = GA.inputs$optim,
        optimArgs = GA.inputs$optimArgs,
        parallel = GA.inputs$parallel,
        run = GA.inputs$run,
        keepBest = GA.inputs$keepBest,
        seed = GA.inputs$seed,
        monitor = GA.inputs$monitor,
        suggestions = GA.inputs$SUGGESTS,
        quiet = GA.inputs$quiet
      )
      rt <- proc.time()[3] - t1
      
    }
    
    multi.GA_nG.o <- multi.GA_nG
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, paste(GA.inputs$parm.type$name, collapse = "."), "_full.rds"))
    
    if(dim(multi.GA_nG@solution)[1] > 1) {
      multi.GA_nG@solution <- t(as.matrix(multi.GA_nG@solution[1,]))
    }
    
    Opt.parm <- GA.opt <- multi.GA_nG@solution
    for (i in 1:GA.inputs$n.layers) {
      if (GA.inputs$surface.type[i] == "cat") {
        ga.p <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        parm <- ga.p / min(ga.p)
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
        
      } else {
        parm <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
      }
    }
    multi.GA_nG@solution <- Opt.parm
    
    RAST <-
      Combine_Surfaces(
        PARM = multi.GA_nG@solution,
        jl.inputs = jl.inputs,
        GA.inputs = GA.inputs,
        rescale = TRUE,
        p.contribution = TRUE
      )
    
    p.cont <- RAST$percent.contribution
    RAST <- RAST$combined.surface
    
    NAME <- paste(GA.inputs$parm.type$name, collapse = ".")
    names(RAST) <- NAME
    
    cd <- suppressWarnings(Run_CS.jl(jl.inputs, RAST, full.mat = TRUE))
    cd.l <- lower(cd)
    cd.l <- cd.l[cd.l != -1]
    cd.l <- scale(cd.l)
    dat <- jl.inputs$df
    dat$cd <- cd.l
    # dat$cd <- scale(lower(cd)[which(lower(cd) != -1)])
    
    write.table(
      cd,
      file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
      sep = ",",
      row.names = F,
      col.names = F
    )
    writeRaster(RAST,
                paste0(GA.inputs$Results.dir, NAME, ".asc"),
                overwrite = TRUE)
    
    ifelse(length(unique(RAST)) > 15,
           type <- "continuous",
           type <- "categorical")
    
    Diagnostic.Plots(
      resistance.mat = dat$cd,
      genetic.dist = jl.inputs$response,
      plot.dir = GA.inputs$Plots.dir,
      type = type,
      name = NAME,
      ID = jl.inputs$ID,
      ZZ = jl.inputs$ZZ
    )
    
    # Get parameter estimates
    MLPE.results <- MLPE.lmm_coef(
      formula = jl.inputs$formula,
      inputs = dat,
      resistance = GA.inputs$Results.dir,
      genetic.dist = jl.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "jl",
      ID = jl.inputs$ID,
      ZZ = jl.inputs$ZZ
    )
    
    # fit.stats <- r.squaredGLMM(
    #   MLPE.lmm(
    #     resistance = lower(cd),
    #     pairwise.genetic = jl.inputs$response,
    #     REML = F,
    #     ID = jl.inputs$ID,
    #     ZZ = jl.inputs$ZZ
    #   )
    # )
    # 
    # aic <- AIC(
    #   MLPE.lmm(
    #     resistance = lower(cd),
    #     pairwise.genetic = jl.inputs$response,
    #     REML = F,
    #     ID = jl.inputs$ID,
    #     ZZ = jl.inputs$ZZ
    #   )
    # )
    # 
    # LL <- logLik(
    #   MLPE.lmm(
    #     resistance = lower(cd),
    #     pairwise.genetic = jl.inputs$response,
    #     REML = F,
    #     ID = jl.inputs$ID,
    #     ZZ = jl.inputs$ZZ
    #   )
    # )
    # 
    # MLPE.model <- MLPE.lmm(
    #   resistance = lower(cd),
    #   pairwise.genetic = jl.inputs$response,
    #   REML = F,
    #   ID = jl.inputs$ID,
    #   ZZ = jl.inputs$ZZ
    # )
    
    fit.mod <- mlpe_rga(formula = jl.inputs$formula,
                        data = dat,
                        ZZ = jl.inputs$ZZ,
                        REML = FALSE)
    fit.mod_REML <- mlpe_rga(formula = jl.inputs$formula,
                             data = dat,
                             ZZ = jl.inputs$ZZ,
                             REML = TRUE)
    
    aic <- suppressWarnings(AIC(
      fit.mod
    ))
    
    fit.stats <- suppressWarnings(r.squaredGLMM(
      fit.mod
    ))
    
    LL <- logLik(
      fit.mod
    )[[1]]
    
    MLPE.model <- fit.mod
    
    if (k.value == 1) {
      k <- 2
    } else if (k.value == 2) {
      k <- GA.inputs$parm.type$n.parm[i] + 
        length(lme4::fixef(fit.mod)) - 1
      
    } else if (k.value == 3) {
      k <- GA.inputs$parm.type$n.parm[i] + 
        length(GA.inputs$layer.names) + 
        length(lme4::fixef(fit.mod)) - 1
      
    } else {
      k <- length(GA.inputs$layer.names[i]) + 
        length(lme4::fixef(fit.mod)) - 1
      
    }
    
    n <- jl.inputs$n.Pops
    AICc <- (-2 * LL) + (2 * k) + ((2 * k) * (k + 1)) / (n - k - 1)
    
    Result.txt(
      GA.results = multi.GA_nG.o,
      GA.inputs = GA.inputs,
      method = "CIRCUITSCAPE.jl",
      Run.Time = rt,
      fit.stats = fit.stats,
      optim = GA.inputs$method,
      k = k,
      aic = aic,
      AICc = AICc,
      LL = LL[[1]],
      fit.mod_REML = fit.mod_REML
    )
    
    write.table(p.cont, file = paste0(GA.inputs$Results.dir, "Percent_Contribution.csv"), sep = ",",
                row.names = F,
                col.names = T)
    
    # save(multi.GA_nG, 
    #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
    
    saveRDS(multi.GA_nG, 
            file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
    
    # unlink(GA.inputs$Write.dir, recursive = T, force = T)
    
    k.df <- data.frame(surface = NAME, k = k)
    
    cd.list <- list(as.matrix(cd))
    names(cd.list) <- NAME
    
    AICc.tab <- data.frame(surface = NAME,
                           obj = multi.GA_nG@fitnessValue,
                           k = k,
                           AIC = aic,
                           AICc = AICc,
                           R2m = fit.stats[[1]],
                           R2c = fit.stats[[2]],
                           LL = LL)
    
    colnames(AICc.tab) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        "k",
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL"
      )
    
    out <- list(GA.summary = multi.GA_nG.o,
                MLPE.model = MLPE.model,
                MLPE.model_REML = fit.mod_REML,
                AICc.tab = AICc.tab,
                cd = cd.list,
                percent.contribution = p.cont,
                k = k.df)
    
    return(out)
  }
}