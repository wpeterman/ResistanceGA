#' Simultaneous optimization of multiple resistance surfaces with kernel smoothing
#'
#' Optimize multiple resistance surfaces simultaneously using genetic algorithms and kernel smoothing
#'
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @return This function optimizes multiple resistance surfaces, returning a Genetic Algorithm (GA) object with summary information. Diagnostic plots of model fit are output to the "Results/Plots" folder that is automatically generated within the folder containing the optimized ASCII files. A text summary of the optimization settings and results is printed to the results folder.
#' @usage MS_optim.scale(CS.inputs, gdist.inputs, GA.inputs)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
MS_optim.scale <- function(CS.inputs = NULL,
                     gdist.inputs = NULL,
                     GA.inputs) {
  
  if (is.null(GA.inputs$scale)) {
    stop(
      "This function should only be used if you intend to apply kernel smoothing to your resistance surfaces"
    )
  }
  
  k.value <- GA.inputs$k.value
  
  if (!is.null(CS.inputs)) {
    if (GA.inputs$parallel != FALSE) {
      warning(
        "\n CIRCUITSCAPE cannot be optimized in parallel. \n Ignoring parallel arguement. \n If you want to optimize in parallel, use least cost paths or commute time with gdistance.",
        immediate. = TRUE
      )
    }
    t1 <- proc.time()[3]
    
    multi.GA_nG <- ga(
      type = "real-valued",
      fitness = Resistance.Opt_multi.scale,
      population = GA.inputs$population,
      selection = GA.inputs$selection,
      mutation = GA.inputs$mutation,
      pcrossover = GA.inputs$pcrossover,
      crossover = GA.inputs$crossover,
      pmutation = GA.inputs$pmutation,
      Min.Max = GA.inputs$Min.Max,
      GA.inputs = GA.inputs,
      CS.inputs = CS.inputs,
      min = GA.inputs$ga.min,
      max = GA.inputs$ga.max,
      popSize = GA.inputs$pop.size,
      maxiter = GA.inputs$maxiter,
      run = GA.inputs$run,
      parallel = FALSE,
      keepBest = GA.inputs$keepBest,
      seed = GA.inputs$seed,
      suggestions = GA.inputs$SUGGESTS,
      quiet = GA.inputs$quiet
    )
    rt <- proc.time()[3] - t1
    
    Opt.parm <- GA.opt <- multi.GA_nG@solution
    for (i in 1:GA.inputs$n.layers) {
    #   if (GA.inputs$surface.type[i] == "cat") {
    #     ga.p <-
    #       GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
    #     parm <- ga.p / min(ga.p)
    #     Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
    #                                                                    1])] <- parm
    #     
    #   } else {
        parm <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
      # }
    }
    multi.GA_nG@solution <- Opt.parm
    
    
    RAST <-
      Combine_Surfaces.scale(
        PARM = multi.GA_nG@solution,
        CS.inputs = CS.inputs,
        GA.inputs = GA.inputs,
        rescale = TRUE
      )
    NAME <- paste(GA.inputs$parm.type$name, collapse = ".")
    names(RAST) <- NAME
    Run_CS(
      CS.inputs,
      GA.inputs,
      r = RAST,
      CurrentMap = FALSE,
      EXPORT.dir = GA.inputs$Results.dir
    )
    
    # ifelse(length(unique(RAST)) > 15,
    #        type <- "continuous",
    #        type <- "categorical")
    
    type <- "continuous"
    
    Diagnostic.Plots(
      resistance.mat = paste0(GA.inputs$Results.dir, NAME, "_resistances.out"),
      genetic.dist = CS.inputs$response,
      plot.dir = GA.inputs$Plots.dir,
      type = type,
      ID = CS.inputs$ID,
      ZZ = CS.inputs$ZZ
    )
    
    fit.stats <-
      r.squaredGLMM(
        MLPE.lmm(
          resistance = paste0(GA.inputs$Results.dir, NAME, "_resistances.out"),
          pairwise.genetic = CS.inputs$response,
          REML = F,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
      )
    
    aic <-
      AIC(
        MLPE.lmm(
          resistance = paste0(GA.inputs$Results.dir, NAME, "_resistances.out"),
          pairwise.genetic = CS.inputs$response,
          REML = F,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
      )
    
    LL <-
      logLik(
        MLPE.lmm(
          resistance = paste0(GA.inputs$Results.dir, NAME, "_resistances.out"),
          pairwise.genetic = CS.inputs$response,
          REML = F,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
      )
    
    MLPE.model <-
      MLPE.lmm(
        resistance = paste0(GA.inputs$Results.dir, NAME, "_resistances.out"),
        pairwise.genetic = CS.inputs$response,
        REML = F,
        ID = CS.inputs$ID,
        ZZ = CS.inputs$ZZ
      )
    
    cd <- (read.table(paste0(
      GA.inputs$Results.dir,
      NAME,
      "_resistances.out"))[-1, -1])
    
    if (k.value == 1) {
      k <- 2
    } else if (k.value == 2) {
      k <-
        sum(GA.inputs$parm.type$n.parm) - sum(GA.inputs$parm.type == "cont") + 1
    } else {
      k <-
        sum(GA.inputs$parm.type$n.parm) - sum(GA.inputs$parm.type == "cont") + nrow(GA.inputs$parm.type) + 1
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
    
    file.remove(list.files(GA.inputs$Write.dir, full.names = TRUE))
    
    k.df <- data.frame(surface = NAME, k = k)
    
    cd.list <- list(as.matrix(cd))
    names(cd.list) <- NAME
    
    out <- list(GA.summary = multi.GA_nG,
                MLPE.model = MLPE.model,
                cd = cd.list,
                k = k.df)
    return(out)
  }
  
  #### Optimize using gdistance ####
  if (!is.null(gdist.inputs)) {
    t1 <- proc.time()[3]
    multi.GA_nG <- ga(
      type = "real-valued",
      fitness = Resistance.Opt_multi.scale,
      population = GA.inputs$population,
      selection = GA.inputs$selection,
      mutation = GA.inputs$mutation,
      pcrossover = GA.inputs$pcrossover,
      crossover = GA.inputs$crossover,
      pmutation = GA.inputs$pmutation,
      Min.Max = GA.inputs$Min.Max,
      GA.inputs = GA.inputs,
      gdist.inputs = gdist.inputs,
      min = GA.inputs$ga.min,
      max = GA.inputs$ga.max,
      popSize = GA.inputs$pop.size,
      maxiter = GA.inputs$maxiter,
      parallel = GA.inputs$parallel,
      run = GA.inputs$run,
      keepBest = GA.inputs$keepBest,
      seed = GA.inputs$seed,
      suggestions = GA.inputs$SUGGESTS,
      quiet = GA.inputs$quiet
    )
    rt <- proc.time()[3] - t1
    
    Opt.parm <- GA.opt <- multi.GA_nG@solution
    for (i in 1:GA.inputs$n.layers) {
      # if (GA.inputs$surface.type[i] == "cat") {
      #   ga.p <-
      #     GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
      #   parm <- ga.p / min(ga.p)
      #   Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
      #                                                                  1])] <- parm
      #   
      # } else {
        parm <-
          GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
        Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                       1])] <- parm
      # }
    }
    multi.GA_nG@solution <- Opt.parm
    
    RAST <-
      Combine_Surfaces.scale(
        PARM = multi.GA_nG@solution,
        gdist.inputs = gdist.inputs,
        GA.inputs = GA.inputs,
        rescale = TRUE
      )
    NAME <- paste(GA.inputs$parm.type$name, collapse = ".")
    names(RAST) <- NAME
    cd <- Run_gdistance(gdist.inputs, RAST)
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
    
    # ifelse(length(unique(RAST)) > 15,
    #        type <- "continuous",
    #        type <- "categorical")
    
    type <- "continuous"
    

    Diagnostic.Plots(
      resistance.mat = cd,
      genetic.dist = gdist.inputs$response,
      plot.dir = GA.inputs$Plots.dir,
      type = type,
      name = NAME,
      ID = gdist.inputs$ID,
      ZZ = gdist.inputs$ZZ
    )
    
    # Get parameter estimates
    MLPE.results <- MLPE.lmm_coef(
      resistance = GA.inputs$Results.dir,
      genetic.dist = gdist.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "gd",
      ID = gdist.inputs$ID,
      ZZ = gdist.inputs$ZZ
    )
    
    fit.stats <- r.squaredGLMM(
      MLPE.lmm(
        resistance = cd,
        pairwise.genetic = gdist.inputs$response,
        REML = F,
        ID = gdist.inputs$ID,
        ZZ = gdist.inputs$ZZ
      )
    )
    
    aic <- AIC(
      MLPE.lmm(
        resistance = cd,
        pairwise.genetic = gdist.inputs$response,
        REML = F,
        ID = gdist.inputs$ID,
        ZZ = gdist.inputs$ZZ
      )
    )
    
    LL <- logLik(
      MLPE.lmm(
        resistance = cd,
        pairwise.genetic = gdist.inputs$response,
        REML = F,
        ID = gdist.inputs$ID,
        ZZ = gdist.inputs$ZZ
      )
    )
    
    MLPE.model <- MLPE.lmm(
      resistance = cd,
      pairwise.genetic = gdist.inputs$response,
      REML = F,
      ID = gdist.inputs$ID,
      ZZ = gdist.inputs$ZZ
    )
    
    if (k.value == 1) {
      k <- 2
    } else if (k.value == 2) {
      k <-
        sum(GA.inputs$parm.type$n.parm) - sum(GA.inputs$parm.type == "cont") + 1
    } else {
      k <-
        sum(GA.inputs$parm.type$n.parm) - sum(GA.inputs$parm.type == "cont") + nrow(GA.inputs$parm.type) + 1
    }
    
    n <- gdist.inputs$n.Pops
    AICc <- (-2 * LL) + (2 * k) + ((2 * k) * (k + 1)) / (n - k - 1)
    
    Result.txt(
      GA.results = multi.GA_nG,
      GA.inputs = GA.inputs,
      method = gdist.inputs$method,
      Run.Time = rt,
      fit.stats = fit.stats,
      optim = GA.inputs$method,
      k = k,
      aic = aic,
      AICc = AICc,
      LL = LL[[1]]
    )
    
    file.remove(list.files(GA.inputs$Write.dir, full.names = TRUE))
    
    k.df <- data.frame(surface = NAME, k = k)
    
    cd.list <- list(as.matrix(cd))
    names(cd.list) <- NAME
    
    out <- list(GA.summary = multi.GA_nG,
                MLPE.model = MLPE.model,
                cd = cd.list,
                k = k.df)
    return(out)
  }
  
  #####  RUN BRENT OPTIMIZATION ####
  # ##  Run second optimization to determine if maximum resistance values should be adjusted
  #   Parm.multiplier <- optim(par=1,
  #                            fn = Max.optim_Brent,
  #                            method = "Brent",
  #                            lower = 0,
  #                            upper = 25,
  #                            GA.inputs = GA.inputs,
  #                            CS.inputs = CS.inputs,
  #                            GA.opt = multi.GA_nG@solution)
  # PARM<-Parm.multiplier$par
  # ###########################################
  # Opt.parm <- GA.opt <- multi.GA_nG@solution
  # for(i in 1:GA.inputs$n.layers){
  #     if(GA.inputs$surface.type[i]=="cat"){
  #       ga.p <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
  #       parm <- ga.p/min(ga.p)#*PARM[1])+1
  # #       parm <- parm/min(parm)
  #       Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
  #
  #     } else {
  #       parm <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
  #       Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
  #     }
  # }
  # # ####
  # multi.GA_nG@solution <- Opt.parm
  # multi.GA_nG@fitnessValue <- Parm.multiplier$value
  # ##################################
  
  #   RAST<-Combine_Surfaces(PARM=multi.GA_nG@solution,CS.inputs=CS.inputs,GA.inputs=GA.inputs)
  #   NAME<-paste(GA.inputs$parm.type$name,collapse=".")
  #   names(RAST)<-NAME
  #   Run_CS(CS.inputs,GA.inputs,r=RAST,CurrentMap=FALSE,EXPORT.dir=GA.inputs$Results.dir)
  #
  #   Diagnostic.Plots(resistance.mat=paste0(GA.inputs$Results.dir,NAME,"_resistances.out"),genetic.dist=CS.inputs$response,plot.dir=GA.inputs$Plots.dir,type="continuous")
  #
  #   # Get parameter estimates
  #   MLPE.results<-MLPE.lmm_coef(resistance=GA.inputs$Results.dir,genetic.dist=CS.inputs$response,out.dir=GA.inputs$Results.dir)
  #
  #   Result.txt(GA.results=multi.GA_nG,GA.inputs=GA.inputs, CS.inputs=CS.inputs)
  #   return(multi.GA_nG)
}