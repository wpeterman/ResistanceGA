#' Single surface optimization
#'
#' Optimize all surfaces contained in a directory using a genetic algorithm executed with the \code{\link[GA]{ga}} function in the Genetic Algorithms package \pkg{GA}
#'
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param jl.inputs Object created from running \code{\link[ResistanceGA]{jl.prep}} function. Defined if optimizing using CIRCUITSCAPE run in Julia
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param dist_mod Logical, if TRUE, a Distance model will be calculated and added to the output table (default = TRUE)
#' @param null_mod Logical, if TRUE, an intercept-only model will be calculated and added to the output table (default = TRUE)
#' @return This function optimizes resistance surfaces in isolation. Following optimization of all surfaces, several summary objects are created.\cr
#' \enumerate{
#' \item Diagnostic plots of model fit are output to the "Results/Plots" directory that is automatically generated within the folder containing the optimized ASCII files.
#' \item A .csv file with the Maximum Likelihood Population Effects mixed effects model coefficient estimates (MLPE_coeff_Table.csv)
#' \item Three summary .csv files are generated: CategoricalResults.csv, ContinuousResults.csv, & All_Results_AICc.csv. These tables contain AICc values and optimization summaries for each surface.
#' }
#' All results tables are also summarized in a named list ($ContinuousResults, $CategoricalResults, $AICc, $MLPE, $MLPE.list, $cd, $k)\cr
#' The \code{lmer} model objects stored $MLPE.list are fit using Restricted Maximum Likelihood \cr
#' $cd is a list of the optimized cost pairwise distance matrices and $k is a table of the surface names and number of parameters used to calculate AICc. These two objects can be passed to \code{\link[ResistanceGA]{Resist.boot}} to conduct a bootstrap analysis.
#' @usage SS_optim(CS.inputs, 
#'  gdist.inputs, 
#'  jl.inputs,
#'  GA.inputs,
#'  dist_mod, 
#'  null_mod)
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @export
#' 
#' @examples  
#' ## Not run:
#' ## *** TO BE COMPLETED *** ##
#' 
#' ## End (Not run)

SS_optim <- function(CS.inputs = NULL,
                     gdist.inputs = NULL,
                     jl.inputs = NULL,
                     GA.inputs,
                     dist_mod = TRUE,
                     null_mod = TRUE) {
  
  if (!is.null(GA.inputs$scale)) {
    stop(
      "This function should NOT be used if you intend to apply kernel smoothing to your resistance surfaces"
    )
  }
  t1 <- proc.time()[3]
  RESULTS.cat <- list() # List to store categorical results within
  RESULTS.cont <- list() # List to store continuous results within
  cnt1 <- 0
  cnt2 <- 0
  k.value <- GA.inputs$k.value
  MLPE.list <- list()
  cd.list <- list()
  k.list <- list()
  ga.list <- list()
  
  wd <- getwd()
  
  # Optimize each surface in turn
  for (i in 1:GA.inputs$n.layers) {
    r <- GA.inputs$Resistance.stack[[i]]
    names(r) <- GA.inputs$layer.names[i]
    
    # >>> CIRCUITSCAPE <<< ------------------------------------------------------------
    
    
    # * Categorical -----------------------------------------------------------
    
    # Processing of categorical surfaces
    if (!is.null(CS.inputs)) {
      # if (GA.inputs$parallel != FALSE) {
      #   warning(
      #     "\n CIRCUITSCAPE cannot be optimized in parallel. \n Ignoring parallel arguement. \n If you want to optimize in parallel, use least cost paths and gdistance.",
      #     immediate. = TRUE
      #   )
      # }
      if (GA.inputs$surface.type[i] == 'cat') {
        cnt1 <- cnt1 + 1
        names(r) <- GA.inputs$layer.names[i]
        
        single.GA <- ga(
          type = "real-valued",
          fitness = Resistance.Opt_single,
          Resistance = r,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          pcrossover = GA.inputs$pcrossover,
          pmutation = GA.inputs$pmutation,
          crossover = GA.inputs$crossover,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.inputs,
          CS.inputs = CS.inputs,
          lower = GA.inputs$min.list[[i]],
          upper = GA.inputs$max.list[[i]],
          parallel = GA.inputs$parallel,
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          elitism = GA.inputs$percent.elite,
          mutation = GA.inputs$mutation,
          # suggestions = GA.inputs$SUGGESTS,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          iter = i,
          quiet = GA.inputs$quiet
        )
        
        if(dim(single.GA@solution)[1] > 1) {
          single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
        }
        
        single.GA@solution <-
          single.GA@solution / min(single.GA@solution)
        df <- data.frame(id = unique(r), t(single.GA@solution))
        r <- subs(r, df)
        NAME <- GA.inputs$layer.names[i]
        names(r) <- NAME
        
        cd <- Run_CS(CS.inputs, 
                     r, 
                     full.mat = TRUE,
                     EXPORT.dir = GA.inputs$Results.dir)
        
        write.table(
          cd,
          file = paste0(GA.inputs$Results.dir, NAME, "_csResistMat.csv"),
          sep = ",",
          row.names = F,
          col.names = F
        )
        writeRaster(r,
                    paste0(GA.inputs$Results.dir, NAME, ".asc"),
                    overwrite = TRUE)
        
        Diagnostic.Plots(
          resistance.mat = lower(cd),
          genetic.dist = CS.inputs$response,
          plot.dir = GA.inputs$Plots.dir,
          type = "categorical",
          name = NAME,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
        
        fit.stats <-
          r.squaredGLMM(
            MLPE.lmm(
              resistance = lower(cd),
              pairwise.genetic = CS.inputs$response,
              REML = F,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        aic <-
          AIC(
            MLPE.lmm(
              resistance = lower(cd),
              pairwise.genetic = CS.inputs$response,
              REML = F,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        LL <-
          logLik(
            MLPE.lmm(
              resistance = lower(cd),
              pairwise.genetic = CS.inputs$response,
              REML = F,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        if (k.value == 1) {
          k <- 2
        } else if (k.value == 2) {
          k <- GA.inputs$parm.type$n.parm[i] + 1
        } else if (k.value == 3) {
          k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
        } else {
          k <- length(GA.inputs$layer.names[i]) + 1
        }
        
        n <- CS.inputs$n.Pops
        AICc <-
          # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
          (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
        
        
        RS <- data.frame(
          GA.inputs$layer.names[i],
          single.GA@fitnessValue,
          k,
          aic,
          AICc,
          fit.stats[[1]],
          fit.stats[[2]],
          LL[[1]],
          single.GA@solution
        )
        
        k <- GA.inputs$parm.type$n.parm[i]
        
        Features <- matrix()
        for (z in 1:(k)) {
          feature <- paste0("Feature", z)
          Features[z] <- feature
        }
        
        colnames(RS) <-
          c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            "k",
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL",
            Features
          )
        
        RESULTS.cat[[cnt1]] <- RS
        
        MLPE.list[[i]] <- MLPE.lmm(
          resistance = lower(cd),
          pairwise.genetic = CS.inputs$response,
          REML = TRUE,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
        
        cd.list[[i]] <- cd
        
        names(MLPE.list)[i] <- GA.inputs$layer.names[i]
        names(cd.list)[i] <- GA.inputs$layer.names[i]
        
        # * Continuous -----------------------------------------------------------
      } else {
        # Processing of continuous surfaces
        cnt2 <- cnt2 + 1
        r <- SCALE(r, 0, 10)
        names(r) <- GA.inputs$layer.names[i]
        
        single.GA <- ga(
          type = "real-valued",
          fitness = Resistance.Opt_single,
          Resistance = r,
          population = GA.inputs$population,
          selection = GA.inputs$selection,
          pcrossover = GA.inputs$pcrossover,
          pmutation = GA.inputs$pmutation,
          crossover = GA.inputs$crossover,
          Min.Max = GA.inputs$Min.Max,
          GA.inputs = GA.inputs,
          CS.inputs = CS.inputs,
          lower = GA.inputs$min.list[[i]],
          upper = GA.inputs$max.list[[i]],
          optim = GA.inputs$optim,
          optimArgs = GA.inputs$optimArgs,
          parallel = GA.inputs$parallel,
          popSize = GA.inputs$pop.size,
          maxiter = GA.inputs$maxiter,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          elitism = GA.inputs$percent.elite,
          mutation = GA.inputs$mutation,
          # suggestions = GA.inputs$SUGGESTS,
          seed = GA.inputs$seed,
          monitor = GA.inputs$monitor,
          iter = i,
          quiet = GA.inputs$quiet
        )
        
        if(dim(single.GA@solution)[1] > 1) {
          single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
        }
        
        ## Deleted nlm function here
        
        if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
          EQ <- get.EQ(9)
          c.names <- dimnames(single.GA@solution)
          single.GA@solution <- t(as.matrix(rep(9, 3)))
          dimnames(single.GA@solution) <- c.names
          
        } else {
          EQ <- get.EQ(single.GA@solution[1])
        }
        
        r.tran <-
          Resistance.tran(
            transformation = single.GA@solution[1],
            shape = single.GA@solution[2],
            max = single.GA@solution[3],
            r = r
          )
        NAME <- GA.inputs$layer.names[i]
        names(r.tran) <- GA.inputs$layer.names[i]
        
        cd <- Run_CS(CS.inputs,
                     r.tran,
                     full.mat = TRUE,
                     EXPORT.dir = GA.inputs$Results.dir)
        
        write.table(
          cd,
          file = paste0(GA.inputs$Results.dir, NAME, "_csResistMat.csv"),
          sep = ",",
          row.names = F,
          col.names = F
        )
        writeRaster(r.tran,
                    paste0(GA.inputs$Results.dir, NAME, ".asc"),
                    overwrite = TRUE)
        
        Diagnostic.Plots(
          resistance.mat = lower(cd),
          genetic.dist = CS.inputs$response,
          plot.dir = GA.inputs$Plots.dir,
          type = "continuous",
          name = NAME,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
        
        Plot.trans(
          PARM = single.GA@solution[-1],
          Resistance = GA.inputs$Resistance.stack[[i]],
          transformation = EQ,
          print.dir = GA.inputs$Plots.dir
        )
        
        fit.stats <- r.squaredGLMM(
          MLPE.lmm(
            REML = F,
            resistance = lower(cd),
            pairwise.genetic = CS.inputs$response,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
        )
        aic <- AIC(
          MLPE.lmm(
            REML = F,
            resistance = lower(cd),
            pairwise.genetic = CS.inputs$response,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
        )
        
        LL <-
          logLik(
            MLPE.lmm(
              resistance = lower(cd),
              pairwise.genetic = CS.inputs$response,
              REML = F,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        if (k.value == 1) {
          k <- 2
        } else if (k.value == 2) {
          k <- GA.inputs$parm.type$n.parm[i] + 1
        } else if (k.value == 3) {
          k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
        } else {
          k <- length(GA.inputs$layer.names[i]) + 1
        }
        
        k.list[[i]] <- k
        names(k.list)[i] <- GA.inputs$layer.names[i]
        
        n <- CS.inputs$n.Pops
        AICc <-
          # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
          (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
        
        
        RS <- data.frame(
          GA.inputs$layer.names[i],
          single.GA@fitnessValue,
          k,
          aic,
          AICc,
          fit.stats[[1]],
          fit.stats[[2]],
          LL[[1]],
          get.EQ(single.GA@solution[1]),
          single.GA@solution[2],
          single.GA@solution[3]
        )
        
        colnames(RS) <-
          c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            "k",
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL",
            "Equation",
            "shape",
            "max"
          )
        RESULTS.cont[[cnt2]] <- RS
        
        MLPE.list[[i]] <- MLPE.lmm(
          resistance = lower(cd),
          pairwise.genetic = CS.inputs$response,
          REML = TRUE,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
        
        cd.list[[i]] <- cd
        
        names(MLPE.list)[i] <- GA.inputs$layer.names[i]
        names(cd.list)[i] <- GA.inputs$layer.names[i]
        
        
        
      } # Close if-else
      if (dist_mod == TRUE) {
        # r <- reclassify(r, c(-Inf, Inf, 1))
        r <- (r * 0) + 1
        names(r) <- "dist"
        cd <- Run_CS(CS.inputs, r, full.mat = T)
        
        write.table(
          cd,
          file = paste0(GA.inputs$Results.dir, "Distance", "_csResistMat.csv"),
          sep = ",",
          row.names = F,
          col.names = F
        )
        
        Dist.AIC <-
          AIC(
            MLPE.lmm(
              resistance = lower(cd),
              pairwise.genetic = CS.inputs$response,
              REML = FALSE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        fit.stats <-
          r.squaredGLMM(
            MLPE.lmm(
              resistance = lower(cd),
              pairwise.genetic = CS.inputs$response,
              REML = FALSE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        LL <-
          logLik(
            MLPE.lmm(
              resistance = lower(cd),
              pairwise.genetic = CS.inputs$response,
              REML = FALSE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        MLPE.list[[i + 1]] <- MLPE.lmm(
          resistance = lower(cd),
          pairwise.genetic = CS.inputs$response,
          REML = TRUE,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
        
        cd.list[[i + 1]] <- cd
        # (read.table(paste0(GA.inputs$Write.dir, "dist_resistances.out"))[-1, -1])
        
        names(cd.list)[i + 1] <- 'Distance'
        
        names(MLPE.list)[i + 1] <- "Distance"
        
        if (GA.inputs$method == "AIC") {
          dist.obj <- Dist.AIC
        } else if (GA.inputs$method == "R2") {
          dist.obj <- fit.stats[[1]]
        } else {
          dist.obj <- LL[[1]]
        }
        
        k <- 2
        k.list[[i + 1]] <- k
        names(k.list)[i + 1] <- 'Distance'
        
        n <- CS.inputs$n.Pops
        AICc <-
          # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
          (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
        
        
        Dist.AICc <-
          data.frame("Distance",
                     dist.obj,
                     k,
                     Dist.AIC,
                     AICc,
                     fit.stats[[1]],
                     fit.stats[[2]],
                     LL[[1]])
        colnames(Dist.AICc) <-
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
      }
      
      if (null_mod == TRUE) {
        response = CS.inputs$response
        
        dat <- data.frame(CS.inputs$ID, response = CS.inputs$response)
        colnames(dat) <- c("pop1", "pop2", "response")
        
        # Fit model
        mod <-
          lFormula(response ~ 1 + (1 | pop1),
                   data = dat,
                   REML = FALSE)
        mod$reTrms$Zt <- CS.inputs$ZZ
        dfun <- do.call(mkLmerDevfun, mod)
        opt <- optimizeLmer(dfun)
        aic <- Null.AIC <-
          AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
        fit.stats <-
          r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
        LL <-
          logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
        
        if (GA.inputs$method == "AIC") {
          null.obj <- Null.AIC
        } else if (GA.inputs$method == "R2") {
          null.obj <- fit.stats[[1]]
        } else {
          null.obj <- LL[[1]]
        }
        k <- 1
        n <- CS.inputs$n.Pops
        AICc <-
          # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
          (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
        
        Null.AICc <-
          data.frame("Null",
                     null.obj,
                     k,
                     Null.AIC,
                     AICc,
                     fit.stats[[1]],
                     fit.stats[[2]],
                     LL[[1]])
        colnames(Null.AICc) <-
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
      }
      
    }
    
    
    # >>> gdistance <<< -------------------------------------------------
    
    if (!is.null(gdist.inputs)) {
      
      # MLPE with Covariates ----------------------------------------------------
      
      if(!is.null(gdist.inputs$covariates)) { 
        
        # * Island GA ---------------------------------------------------------------
        
        
        if(isTRUE(GA.inputs$gaisl)) {
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              single.GA@solution / min(single.GA@solution)
            df <- data.frame(id = unique(r), t(single.GA@solution))
            r <- subs(r, df)
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            scale.cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "categorical",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
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
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
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
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
            
            # *-* Continuous -----------------------------------------------------------
            
          } else {
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
            fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.mod_REML <- mlpe_rga(formula = gdist.inputs$formula,
                                     data = dat,
                                     ZZ = gdist.inputs$ZZ,
                                     REML = TRUE)
            
            aic <- Dist.AIC <- suppressWarnings(AIC(
              fit.mod
            ))
            
            fit.stats <- r.squaredGLMM(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            MLPE.list[[i + 1]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
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
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close gaisl cat-cont if else
        } else { # * Standard GA -------------------------------------------------------------
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              parallel = GA.inputs$parallel,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = unique(r), t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = unique(GA.inputs$Resistance.stack[[i]]), t(single.GA@solution))
            
            r <- subs(GA.inputs$Resistance.stack[[i]], df) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "categorical",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # if (k.value == 1) {
            #   k <- 2
            # } else if (k.value == 2) {
            #   k <- GA.inputs$parm.type$n.parm[i] + 1
            # } else if (k.value == 3) {
            #   k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            # } else {
            #   k <- length(GA.inputs$layer.names[i]) + 1
            # }
            
            fit.mod <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.mod_REML <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                     data = dat,
                                     ZZ = gdist.inputs$ZZ,
                                     REML = TRUE)
            
            fit.stats <- suppressWarnings(r.squaredGLMM(
              fit.mod
            ))
            
            aic <- AIC(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
          } 
          else { # *-* Continuous ----------
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
            # fit.stats <-
            #   r.squaredGLMM(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # 
            # aic <-
            #   AIC(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # LL <-
            #   logLik(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            
            fit.mod <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.mod_REML <- mlpe_rga(formula = gd ~ cd + (1|pop),
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
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close cat-cont if else
          
        } # Close GA vs. gaisl if-else
        
        if (dist_mod == TRUE) {
          rd <- reclassify(r, c(-Inf, Inf, 1))
          names(rd) <- "dist"
          cd <- Run_gdistance(gdist.inputs, rd)
          
          dat <- gdist.inputs$df
          dat$cd <- scale(c(cd))
          
          write.table(
            as.matrix(cd),
            file = paste0(GA.inputs$Results.dir, 'Distance', "_", gdist.inputs$method, "_distMat.csv"),
            sep = ",",
            row.names = F,
            col.names = F
          )
          
          fit.mod <- mlpe_rga(formula = gd ~ cd + (1|pop),
                              data = dat,
                              ZZ = gdist.inputs$ZZ,
                              REML = FALSE)
          
          fit.mod_REML <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                   data = dat,
                                   ZZ = gdist.inputs$ZZ,
                                   REML = TRUE)
          
          aic <- Dist.AIC <- suppressWarnings(AIC(
            fit.mod
          ))
          
          fit.stats <- r.squaredGLMM(
            fit.mod
          )
          
          LL <- logLik(
            fit.mod
          )[[1]]
          
          MLPE.list[[i + 1]] <-  fit.mod_REML
          
          cd.list[[i + 1]] <- as.matrix(cd)
          names(cd.list)[i + 1] <- "Distance"
          
          names(MLPE.list)[i + 1] <- "Distance"
          
          ROW <- nrow(gdist.inputs$ID)
          k <- 2
          
          k.list[[i + 1]] <- k
          names(k.list)[i + 1] <- 'Distance'
          
          if (GA.inputs$method == "AIC") {
            dist.obj <- -Dist.AIC
          } else if (GA.inputs$method == "R2") {
            dist.obj <- fit.stats[[1]]
          } else {
            dist.obj <- LL[[1]]
          }
          
          n <- gdist.inputs$n.Pops
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          Dist.AICc <- data.frame("Distance",
                                  dist.obj,
                                  k,
                                  Dist.AIC,
                                  AICc,
                                  fit.stats[[1]],
                                  fit.stats[[2]],
                                  LL[[1]])
          colnames(Dist.AICc) <- c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            'k',
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL"
          )
        }
        
        if (null_mod == TRUE) {
          dat <- data.frame(gdist.inputs$ID, response = gdist.inputs$response)
          colnames(dat) <- c("pop1", "pop2", "response")
          
          # Fit model
          mod <-
            lFormula(response ~ 1 + (1 | pop1),
                     data = dat,
                     REML = FALSE)
          mod$reTrms$Zt <- gdist.inputs$ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          aic <- Null.AIC <-
            AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          fit.stats <-
            r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          LL <-
            logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          ROW <- nrow(gdist.inputs$ID)
          k <- 1
          
          if (GA.inputs$method == "AIC") {
            null.obj <- -Null.AIC
          } else if (GA.inputs$method == "R2") {
            null.obj <- fit.stats[[1]]
          } else {
            null.obj <- LL[[1]]
          }
          n <- gdist.inputs$n.Pops
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          Null.AICc <-
            data.frame("Null",
                       null.obj,
                       k,
                       Null.AIC,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Null.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
      } # End Covariate
      
      # MLPE no Covariates ------------------------------------------------------
      if(is.null(gdist.inputs$covariates)) {
        
        # * Island GA ---------------------------------------------------------------
        
        
        if(isTRUE(GA.inputs$gaisl)) {
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              single.GA@solution / min(single.GA@solution)
            df <- data.frame(id = unique(r), t(single.GA@solution))
            r <- subs(r, df)
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "categorical",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # if (k.value == 1) {
            #   k <- 2
            # } else if (k.value == 2) {
            #   k <- GA.inputs$parm.type$n.parm[i] + 1
            # } else if (k.value == 3) {
            #   k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            # } else {
            #   k <- length(GA.inputs$layer.names[i]) + 1
            # }
            
            fit.mod <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.stats <- suppressWarnings(r.squaredGLMM(
              fit.mod
            ))
            
            aic <- AIC(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = TRUE,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
            
            # *-* Continuous -----------------------------------------------------------
            
          } else {
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
            # fit.stats <-
            #   r.squaredGLMM(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # 
            # aic <-
            #   AIC(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # LL <-
            #   logLik(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            
            fit.mod <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.stats <- suppressWarnings(r.squaredGLMM(
              fit.mod
            ))
            
            aic <- AIC(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            MLPE.list[[i]] <-  MLPE.lmm2(
              resistance = cd,
              response = gdist.inputs$response,
              REML = TRUE,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # if (k.value == 1) {
            #   k <- 2
            # } else if (k.value == 2) {
            #   k <- GA.inputs$parm.type$n.parm[i] + 1
            # } else if (k.value == 3) {
            #   k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            # } else {
            #   k <- length(GA.inputs$layer.names[i]) + 1
            # }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close gaisl cat-cont if else
        } else { # * Standard GA -------------------------------------------------------------
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              parallel = GA.inputs$parallel,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = unique(r), t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = unique(GA.inputs$Resistance.stack[[i]]), t(single.GA@solution))
            r <- subs(GA.inputs$Resistance.stack[[i]], df) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "categorical",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # if (k.value == 1) {
            #   k <- 2
            # } else if (k.value == 2) {
            #   k <- GA.inputs$parm.type$n.parm[i] + 1
            # } else if (k.value == 3) {
            #   k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            # } else {
            #   k <- length(GA.inputs$layer.names[i]) + 1
            # }
            
            fit.mod <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.mod_REML <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                     data = dat,
                                     ZZ = gdist.inputs$ZZ,
                                     REML = TRUE)
            
            fit.stats <- suppressWarnings(r.squaredGLMM(
              fit.mod
            ))
            
            aic <- AIC(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
          } 
          else { # *-* Continuous ----------
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
            # fit.stats <-
            #   r.squaredGLMM(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # 
            # aic <-
            #   AIC(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            # 
            # LL <-
            #   logLik(
            #     MLPE.lmm(
            #       resistance = cd,
            #       pairwise.genetic = gdist.inputs$response,
            #       REML = F,
            #       ID = gdist.inputs$ID,
            #       ZZ = gdist.inputs$ZZ
            #     )
            #   )
            
            fit.mod <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.mod_REML <- mlpe_rga(formula = gd ~ cd + (1|pop),
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
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close cat-cont if else
          
        } # Close GA vs. gaisl if-else
        
        if (dist_mod == TRUE) {
          r <- reclassify(r, c(-Inf, Inf, 1))
          names(r) <- "dist"
          cd <- Run_gdistance(gdist.inputs, r)
          
          dat <- gdist.inputs$df
          dat$cd <- scale(c(cd))
          
          write.table(
            as.matrix(cd),
            file = paste0(GA.inputs$Results.dir, 'Distance', "_", gdist.inputs$method, "_distMat.csv"),
            sep = ",",
            row.names = F,
            col.names = F
          )
          
          fit.mod <- mlpe_rga(formula = gd ~ cd + (1|pop),
                              data = dat,
                              ZZ = gdist.inputs$ZZ,
                              REML = FALSE)
          
          fit.mod_REML <- mlpe_rga(formula = gd ~ cd + (1|pop),
                                   data = dat,
                                   ZZ = gdist.inputs$ZZ,
                                   REML = TRUE)
          
          aic <- Dist.AIC <- suppressWarnings(AIC(
            fit.mod
          ))
          
          fit.stats <- r.squaredGLMM(
            fit.mod
          )
          
          LL <- logLik(
            fit.mod
          )[[1]]
          
          MLPE.list[[i + 1]] <-  fit.mod_REML
          
          cd.list[[i + 1]] <- as.matrix(cd)
          names(cd.list)[i + 1] <- "Distance"
          
          names(MLPE.list)[i + 1] <- "Distance"
          
          ROW <- nrow(gdist.inputs$ID)
          k <- 2
          
          k.list[[i + 1]] <- k
          names(k.list)[i + 1] <- 'Distance'
          
          if (GA.inputs$method == "AIC") {
            dist.obj <- -Dist.AIC
          } else if (GA.inputs$method == "R2") {
            dist.obj <- fit.stats[[1]]
          } else {
            dist.obj <- LL[[1]]
          }
          
          n <- gdist.inputs$n.Pops
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          Dist.AICc <- data.frame("Distance",
                                  dist.obj,
                                  k,
                                  Dist.AIC,
                                  AICc,
                                  fit.stats[[1]],
                                  fit.stats[[2]],
                                  LL[[1]])
          colnames(Dist.AICc) <- c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            'k',
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL"
          )
        }
        
        if (null_mod == TRUE) {
          dat <- data.frame(gdist.inputs$ID, response = gdist.inputs$response)
          colnames(dat) <- c("pop1", "pop2", "response")
          
          # Fit model
          mod <-
            lFormula(response ~ 1 + (1 | pop1),
                     data = dat,
                     REML = FALSE)
          mod$reTrms$Zt <- gdist.inputs$ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          aic <- Null.AIC <-
            AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          fit.stats <-
            r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          LL <-
            logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          ROW <- nrow(gdist.inputs$ID)
          k <- 1
          
          if (GA.inputs$method == "AIC") {
            null.obj <- -Null.AIC
          } else if (GA.inputs$method == "R2") {
            null.obj <- fit.stats[[1]]
          } else {
            null.obj <- LL[[1]]
          }
          n <- gdist.inputs$n.Pops
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          Null.AICc <-
            data.frame("Null",
                       null.obj,
                       k,
                       Null.AIC,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Null.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
        
      } # End no covariate
    } # End gdistance
    
    # >>> Julia <<< -------------------------------------------------
    if (!is.null(jl.inputs)) {

      setwd(jl.inputs$JULIA_HOME)
      # MLPE with Covariates ----------------------------------------------------
      
      if(!is.null(jl.inputs$covariates)) { 
        
        # * Island GA ---------------------------------------------------------------
        
        
        if(isTRUE(GA.inputs$gaisl)) {
          stop("This feature is not currently supported")
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = unique(r), t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = unique(GA.inputs$Resistance.stack[[i]]), t(single.GA@solution))
            r <- subs(GA.inputs$Resistance.stack[[i]], df) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- Run_gdistance(gdist.inputs, r)
            scale.cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method,  "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "categorical",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
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
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = cd,
            #     response = gdist.inputs$response,
            #     REML = F,
            #     ID = gdist.inputs$ID,
            #     ZZ = gdist.inputs$ZZ
            #   )
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
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
            
            # *-* Continuous -----------------------------------------------------------
            
          } else {
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              gdist.inputs = gdist.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- Run_gdistance(gdist.inputs, r)
            dat <- gdist.inputs$df
            dat$cd <- scale(c(cd))
            
            write.table(
              as.matrix(cd),
              file = paste0(GA.inputs$Results.dir, NAME, "_", gdist.inputs$method, "_distMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = cd,
              genetic.dist = gdist.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = gdist.inputs$ID,
              ZZ = gdist.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
            fit.mod <- mlpe_rga(formula = gdist.inputs$formula,
                                data = dat,
                                ZZ = gdist.inputs$ZZ,
                                REML = FALSE)
            
            fit.mod_REML <- mlpe_rga(formula = gdist.inputs$formula,
                                     data = dat,
                                     ZZ = gdist.inputs$ZZ,
                                     REML = TRUE)
            
            aic <- Dist.AIC <- suppressWarnings(AIC(
              fit.mod
            ))
            
            fit.stats <- r.squaredGLMM(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            MLPE.list[[i + 1]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
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
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- gdist.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close gaisl cat-cont if else
        } else { # * Standard GA -------------------------------------------------------------
          
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              parallel = GA.inputs$parallel,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = unique(r), t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = unique(GA.inputs$Resistance.stack[[i]]), t(single.GA@solution))
            r <- subs(GA.inputs$Resistance.stack[[i]], df) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            # save(cd, file = paste0(GA.inputs$Write.dir, NAME, ".rda"))
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = dat$cd,
              genetic.dist = jl.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "categorical",
              name = NAME,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ
            )
            
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
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <- fit.mod_REML
            
            cd.list[[i]] <- as.matrix(cd)
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            # rm(single.GA, r)
            gc()
          } 
          else { # *-* Continuous ----------
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single.cov,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            # save(single.GA, 
            #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = dat$cd,
              genetic.dist = jl.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
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
              k <- GA.inputs$parm.type$n.parm[i] + length(lme4::fixef(fit.mod)) - 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + length(lme4::fixef(fit.mod)) - 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + length(lme4::fixef(fit.mod)) - 1
            }
            
            MLPE.list[[i]] <- fit.mod_REML
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
            # rm(single.GA, r)
            gc()
          } # Close cat-cont if else
          
        } # Close GA vs. gaisl if-else
        
        if (dist_mod == TRUE) {
          r <- reclassify(r, c(-Inf, Inf, 1))
          names(r) <- "dist"
          cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
          cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
          dat <- jl.inputs$df
          dat$cd <- cd.l
          
          write.table(
            cd,
            file = paste0(GA.inputs$Results.dir, "Distance", "_jlResistMat.csv"),
            sep = ",",
            row.names = F,
            col.names = F
          )
          
          fit.mod <- mlpe_rga(formula = jl.inputs$formula,
                              data = dat,
                              ZZ = jl.inputs$ZZ,
                              REML = FALSE)
          fit.mod_REML <- mlpe_rga(formula = jl.inputs$formula,
                                   data = dat,
                                   ZZ = jl.inputs$ZZ,
                                   REML = TRUE)
          
          aic <- Dist.AIC <- suppressWarnings(AIC(
            fit.mod)   
          )
          
          fit.stats <- suppressWarnings(r.squaredGLMM(
            fit.mod
          ))
          
          LL <- logLik(
            fit.mod
          )[[1]]
          
          MLPE.list[[i + 1]] <-  fit.mod_REML
          
          cd.list[[i + 1]] <- cd
          names(cd.list)[i + 1] <- "Distance"
          
          names(MLPE.list)[i + 1] <- "Distance"
          
          ROW <- nrow(jl.inputs$ID)
          k <- length(lme4::fixef(fit.mod))
          
          k.list[[i + 1]] <- k
          names(k.list)[i + 1] <- 'Distance'
          
          if (GA.inputs$method == "AIC") {
            dist.obj <- -Dist.AIC
          } else if (GA.inputs$method == "R2") {
            dist.obj <- fit.stats[[1]]
          } else {
            dist.obj <- LL[[1]]
          }
          
          n <- jl.inputs$n.Pops
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          Dist.AICc <- data.frame("Distance",
                                  dist.obj,
                                  k,
                                  Dist.AIC,
                                  AICc,
                                  fit.stats[[1]],
                                  fit.stats[[2]],
                                  LL[[1]])
          colnames(Dist.AICc) <- c(
            "Surface",
            paste0("obj.func_", GA.inputs$method),
            'k',
            "AIC",
            "AICc",
            "R2m",
            "R2c",
            "LL"
          )
        }
        
        if (null_mod == TRUE) {
          dat <- data.frame(jl.inputs$ID, response = jl.inputs$response)
          colnames(dat) <- c("pop1", "pop2", "response")
          
          # Fit model
          mod <-
            lFormula(response ~ 1 + (1 | pop1),
                     data = dat,
                     REML = FALSE)
          mod$reTrms$Zt <- jl.inputs$ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          aic <- Null.AIC <-
            AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          fit.stats <-
            r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          LL <-
            logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          ROW <- nrow(jl.inputs$ID)
          k <- 1
          
          if (GA.inputs$method == "AIC") {
            null.obj <- -Null.AIC
          } else if (GA.inputs$method == "R2") {
            null.obj <- fit.stats[[1]]
          } else {
            null.obj <- LL[[1]]
          }
          n <- jl.inputs$n.Pops
          AICc <-
            # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
            (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
          
          Null.AICc <-
            data.frame("Null",
                       null.obj,
                       k,
                       Null.AIC,
                       AICc,
                       fit.stats[[1]],
                       fit.stats[[2]],
                       LL[[1]])
          colnames(Null.AICc) <-
            c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
        }
      } # End Covariate
      
      
      # MLPE no Covariates ------------------------------------------------------
      if(is.null(jl.inputs$covariates)) { 
        
        
        # * Island GA ---------------------------------------------------------------
        
        if(isTRUE(GA.inputs$gaisl)) {
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = unique(r), t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = unique(GA.inputs$Resistance.stack[[i]]), t(single.GA@solution))
            r <- subs(GA.inputs$Resistance.stack[[i]], df) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = dat$cd,
              genetic.dist = jl.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "categorical",
              name = NAME,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ
            )
            
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
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- cd
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
          } else { # *-* Continuous ---------------
            
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- gaisl(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              numIslands = GA.inputs$numIslands,
              migrationRate = GA.inputs$migrationRate,
              migrationInterval = GA.inputs$migrationInterval,
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              # keepBest = GA.inputs$keepBest,
              # suggestions = GA.inputs$SUGGESTS,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, GA.inputs$layer.names[i], "_full.rds"))
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = dat$cd,
              genetic.dist = jl.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
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
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- cd
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
          } # Close if-else
          
          if (dist_mod == TRUE) {
            r <- reclassify(r, c(-Inf, Inf, 1))
            names(r) <- "dist"
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, "Distance", "_jlResistMat.csv"),
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            fit.mod <- mlpe_rga(formula = jl.inputs$formula,
                                data = dat,
                                ZZ = jl.inputs$ZZ,
                                REML = FALSE)
            fit.mod_REML <- mlpe_rga(formula = jl.inputs$formula,
                                     data = dat,
                                     ZZ = jl.inputs$ZZ,
                                     REML = TRUE)
            
            aic <- Dist.AIC <- suppressWarnings(AIC(
              fit.mod
            ))
            
            fit.stats <- r.squaredGLMM(
              fit.mod
            )
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            MLPE.list[[i + 1]] <-  fit.mod_REML
            
            cd.list[[i + 1]] <- cd
            names(cd.list)[i + 1] <- "Distance"
            
            names(MLPE.list)[i + 1] <- "Distance"
            
            ROW <- nrow(jl.inputs$ID)
            k <- 2
            
            k.list[[i + 1]] <- k
            names(k.list)[i + 1] <- 'Distance'
            
            if (GA.inputs$method == "AIC") {
              dist.obj <- -Dist.AIC
            } else if (GA.inputs$method == "R2") {
              dist.obj <- fit.stats[[1]]
            } else {
              dist.obj <- LL[[1]]
            }
            
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            Dist.AICc <- data.frame("Distance",
                                    dist.obj,
                                    k,
                                    Dist.AIC,
                                    AICc,
                                    fit.stats[[1]],
                                    fit.stats[[2]],
                                    LL[[1]])
            colnames(Dist.AICc) <- c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
          }
          
          if (null_mod == TRUE) {
            dat <- data.frame(jl.inputs$ID, response = jl.inputs$response)
            colnames(dat) <- c("pop1", "pop2", "response")
            
            # Fit model
            mod <-
              lFormula(response ~ 1 + (1 | pop1),
                       data = dat,
                       REML = FALSE)
            mod$reTrms$Zt <- jl.inputs$ZZ
            dfun <- do.call(mkLmerDevfun, mod)
            opt <- optimizeLmer(dfun)
            aic <- Null.AIC <-
              AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
            fit.stats <-
              r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
            LL <-
              logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
            ROW <- nrow(jl.inputs$ID)
            k <- 1
            
            if (GA.inputs$method == "AIC") {
              null.obj <- -Null.AIC
            } else if (GA.inputs$method == "R2") {
              null.obj <- fit.stats[[1]]
            } else {
              null.obj <- LL[[1]]
            }
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            Null.AICc <-
              data.frame("Null",
                         null.obj,
                         k,
                         Null.AIC,
                         AICc,
                         fit.stats[[1]],
                         fit.stats[[2]],
                         LL[[1]])
            colnames(Null.AICc) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL"
              )
          }
          
          
        } else { # * Standard GA -----------------
          # *-* Categorical -----------------------------------------------------------
          
          if (GA.inputs$surface.type[i] == 'cat') {
            cnt1 <- cnt1 + 1
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            if(dim(single.GA@solution)[1] > 1) {
              single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
            }
            
            single.GA@solution <-
              single.GA@solution / min(single.GA@solution)
            # df <- data.frame(id = unique(r), t(single.GA@solution))
            # r <- subs(r, df)
            df <- data.frame(id = unique(GA.inputs$Resistance.stack[[i]]), t(single.GA@solution))
            r <- subs(GA.inputs$Resistance.stack[[i]], df) ## Modified 2019-03-26
            NAME <- GA.inputs$layer.names[i]
            names(r) <- NAME
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              sep = ",",
              row.names = F,
              col.names = F
            )
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = dat$cd,
              genetic.dist = jl.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "categorical",
              name = NAME,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ
            )
            
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
            
            # fit.stats <- r.squaredGLMM(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            # 
            # aic <- AIC(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            # 
            # LL <- logLik(
            #   MLPE.lmm2(
            #     resistance = dat$cd,
            #     response = jl.inputs$response,
            #     REML = F,
            #     ID = jl.inputs$ID,
            #     ZZ = jl.inputs$ZZ
            #   )
            # )
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              single.GA@solution
            )
            
            k <- GA.inputs$parm.type$n.parm[i]
            
            Features <- matrix()
            for (z in 1:(k)) {
              feature <- paste0("Feature", z)
              Features[z] <- feature
            }
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                "k",
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                Features
              )
            
            RESULTS.cat[[cnt1]] <- RS
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- cd
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
          } else { # *-* Continuous ---------------
            
            # Processing of continuous surfaces
            cnt2 <- cnt2 + 1
            r <- SCALE(r, 0, 10)
            names(r) <- GA.inputs$layer.names[i]
            
            single.GA <- ga(
              type = "real-valued",
              fitness = Resistance.Opt_single,
              Resistance = r,
              population = GA.inputs$population,
              selection = GA.inputs$selection,
              pcrossover = GA.inputs$pcrossover,
              pmutation = GA.inputs$pmutation,
              crossover = GA.inputs$crossover,
              Min.Max = GA.inputs$Min.Max,
              GA.inputs = GA.inputs,
              jl.inputs = jl.inputs,
              lower = GA.inputs$min.list[[i]],
              upper = GA.inputs$max.list[[i]],
              optim = GA.inputs$optim,
              optimArgs = GA.inputs$optimArgs,
              parallel = GA.inputs$parallel,
              popSize = GA.inputs$pop.size,
              maxiter = GA.inputs$maxiter,
              run = GA.inputs$run,
              keepBest = GA.inputs$keepBest,
              elitism = GA.inputs$percent.elite,
              mutation = GA.inputs$mutation,
              # suggestions = GA.inputs$SUGGESTS,
              seed = GA.inputs$seed,
              monitor = GA.inputs$monitor,
              iter = i,
              quiet = GA.inputs$quiet
            )
            
            ## Deleted nlm function here
            
            if(single.GA@fitnessValue == -99999 | dim(single.GA@solution)[1] > 1) {
              EQ <- get.EQ(9)
              c.names <- dimnames(single.GA@solution)
              single.GA@solution <- t(as.matrix(rep(9, 3)))
              dimnames(single.GA@solution) <- c.names
              
            } else {
              EQ <- get.EQ(single.GA@solution[1])
            }
            
            r <-
              Resistance.tran(
                transformation = single.GA@solution[1],
                shape = single.GA@solution[2],
                max = single.GA@solution[3],
                r = r
              )
            names(r) <- GA.inputs$layer.names[i]
            NAME <- GA.inputs$layer.names[i]
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, NAME, "_jlResistMat.csv"),
              
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            writeRaster(r,
                        paste0(GA.inputs$Results.dir, NAME, ".asc"),
                        overwrite = TRUE)
            
            saveRDS(single.GA, 
                    file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
            
            ga.list[[i]] <- single.GA
            names(ga.list[i]) <- NAME
            
            Diagnostic.Plots(
              resistance.mat = dat$cd,
              genetic.dist = jl.inputs$response,
              plot.dir = GA.inputs$Plots.dir,
              type = "continuous",
              name = NAME,
              ID = jl.inputs$ID,
              ZZ = jl.inputs$ZZ
            )
            
            Plot.trans(
              PARM = single.GA@solution[-1],
              Resistance = GA.inputs$Resistance.stack[[i]],
              transformation = EQ,
              print.dir = GA.inputs$Plots.dir
            )
            
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
            
            # fit.stats <-
            #   r.squaredGLMM(
            #     MLPE.lmm(
            #       resistance = dat$cd,
            #       pairwise.genetic = jl.inputs$response,
            #       REML = F,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ
            #     )
            #   )
            # 
            # 
            # aic <-
            #   AIC(
            #     MLPE.lmm(
            #       resistance = dat$cd,
            #       pairwise.genetic = jl.inputs$response,
            #       REML = F,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ
            #     )
            #   )
            # 
            # LL <-
            #   logLik(
            #     MLPE.lmm(
            #       resistance = dat$cd,
            #       pairwise.genetic = jl.inputs$response,
            #       REML = F,
            #       ID = jl.inputs$ID,
            #       ZZ = jl.inputs$ZZ
            #     )
            #   )
            
            MLPE.list[[i]] <-  fit.mod_REML
            
            cd.list[[i]] <- cd
            names(cd.list)[i] <- GA.inputs$layer.names[i]
            
            names(MLPE.list)[i] <- GA.inputs$layer.names[i]
            
            if (k.value == 1) {
              k <- 2
            } else if (k.value == 2) {
              k <- GA.inputs$parm.type$n.parm[i] + 1
            } else if (k.value == 3) {
              k <- GA.inputs$parm.type$n.parm[i] + length(GA.inputs$layer.names) + 1
            } else {
              k <- length(GA.inputs$layer.names[i]) + 1
            }
            
            k.list[[i]] <- k
            names(k.list)[i] <- GA.inputs$layer.names[i]
            
            n <- jl.inputs$n.Pops
            
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            
            RS <- data.frame(
              GA.inputs$layer.names[i],
              single.GA@fitnessValue,
              k,
              aic,
              AICc,
              fit.stats[[1]],
              fit.stats[[2]],
              LL[[1]],
              get.EQ(single.GA@solution[1]),
              single.GA@solution[2],
              single.GA@solution[3]
            )
            
            colnames(RS) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL",
                "Equation",
                "shape",
                "max"
              )
            RESULTS.cont[[cnt2]] <- RS
            
          } # Close  cat-cont ifelse
          
          if (dist_mod == TRUE) {
            r <- reclassify(r, c(-Inf, Inf, 1))
            names(r) <- "dist"
            
            cd <- suppressWarnings(Run_CS.jl(jl.inputs, r, full.mat = TRUE))
            cd.l <- scale(lower(cd)[jl.inputs$keep == 1])
            dat <- jl.inputs$df
            dat$cd <- cd.l
            
            write.table(
              cd,
              file = paste0(GA.inputs$Results.dir, "Distance", "_jlResistMat.csv"),
              sep = ",",
              row.names = F,
              col.names = F
            )
            
            fit.mod <- mlpe_rga(formula = jl.inputs$formula,
                                data = dat,
                                ZZ = jl.inputs$ZZ,
                                REML = FALSE)
            fit.mod_REML <- mlpe_rga(formula = jl.inputs$formula,
                                     data = dat,
                                     ZZ = jl.inputs$ZZ,
                                     REML = TRUE)
            
            aic <- Dist.AIC <- suppressWarnings(AIC(
              fit.mod)   
            )
            
            fit.stats <- suppressWarnings(r.squaredGLMM(
              fit.mod
            ))
            
            LL <- logLik(
              fit.mod
            )[[1]]
            
            MLPE.list[[i + 1]] <-  fit.mod_REML
            
            cd.list[[i + 1]] <- cd
            names(cd.list)[i + 1] <- "Distance"
            
            names(MLPE.list)[i + 1] <- "Distance"
            
            ROW <- nrow(jl.inputs$ID)
            k <- 2
            
            k.list[[i + 1]] <- k
            names(k.list)[i + 1] <- 'Distance'
            
            if (GA.inputs$method == "AIC") {
              dist.obj <- -Dist.AIC
            } else if (GA.inputs$method == "R2") {
              dist.obj <- fit.stats[[1]]
            } else {
              dist.obj <- LL[[1]]
            }
            
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            Dist.AICc <- data.frame("Distance",
                                    dist.obj,
                                    k,
                                    Dist.AIC,
                                    AICc,
                                    fit.stats[[1]],
                                    fit.stats[[2]],
                                    LL[[1]])
            colnames(Dist.AICc) <- c(
              "Surface",
              paste0("obj.func_", GA.inputs$method),
              'k',
              "AIC",
              "AICc",
              "R2m",
              "R2c",
              "LL"
            )
          }
          
          if (null_mod == TRUE) {
            dat <- data.frame(jl.inputs$ID, response = jl.inputs$response)
            colnames(dat) <- c("pop1", "pop2", "response")
            
            # Fit model
            mod <-
              lFormula(response ~ 1 + (1 | pop1),
                       data = dat,
                       REML = FALSE)
            mod$reTrms$Zt <- jl.inputs$ZZ
            dfun <- do.call(mkLmerDevfun, mod)
            opt <- optimizeLmer(dfun)
            aic <- Null.AIC <-
              AIC(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
            fit.stats <-
              r.squaredGLMM(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
            LL <-
              logLik(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
            ROW <- nrow(jl.inputs$ID)
            k <- 1
            
            if (GA.inputs$method == "AIC") {
              null.obj <- -Null.AIC
            } else if (GA.inputs$method == "R2") {
              null.obj <- fit.stats[[1]]
            } else {
              null.obj <- LL[[1]]
            }
            n <- jl.inputs$n.Pops
            AICc <-
              # (-2 * LL) + (((2 * k) * (k + 1)) / (n - k - 1))
              (aic) + (((2 * k) * (k + 1)) / max((n - k - 1), 1))
            
            Null.AICc <-
              data.frame("Null",
                         null.obj,
                         k,
                         Null.AIC,
                         AICc,
                         fit.stats[[1]],
                         fit.stats[[2]],
                         LL[[1]])
            colnames(Null.AICc) <-
              c(
                "Surface",
                paste0("obj.func_", GA.inputs$method),
                'k',
                "AIC",
                "AICc",
                "R2m",
                "R2c",
                "LL"
              )
          }
          
        } # End Island-Standard ifelse
      } # End no covariates
    } # End julia
  } # Close ascii loop
  
  
  # Optimization summary ----------------------------------------------------
  
  # Make results data frame
  Results.cat <- data.frame()
  Results.cont <- data.frame()
  # cnt1<-0
  # cnt2<-0
  for (i in 1:GA.inputs$n.layers) {
    if (GA.inputs$surface.type[i] == 'cat') {
      #     cnt1 <- cnt1+1
      #     RS <- data.frame(GA.inputs$layer.names[i], -(RESULTS.cat[[i]]@fitnessValue),RESULTS[[i]]@solution)
      Results.cat <- do.call(rbind.fill, RESULTS.cat)
    } else {
      #   cnt2 <-cnt2+1
      #   RS <- data.frame(GA.inputs$layer.names[i], -(RESULTS.cont[[i]]@fitnessValue), Cont.Param(RESULTS[[i]]@solution))
      Results.cont <- do.call(rbind, RESULTS.cont)
    }
  }
  
  
  # Compile results into tables
  cat("\n")
  cat("\n")
  if (nrow(Results.cat) > 0) {
    Features <- array()
    for (i in 1:ncol(Results.cat) - 8) {
      feature <- paste0("Feature", i)
      Features[i] <- feature
    }
    colnames(Results.cat) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        'k',
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL",
        Features
      )
    Results.cat <-  Results.cat[order(Results.cat$AICc), ]
    write.table(
      Results.cat,
      paste0(GA.inputs$Results.dir, "CategoricalResults.csv"),
      sep = ",",
      col.names = T,
      row.names = F
    )
  }
  
  if (ncol(Results.cont) > 0) {
    colnames(Results.cont) <-
      c(
        "Surface",
        paste0("obj.func_", GA.inputs$method),
        'k',
        "AIC",
        "AICc",
        "R2m",
        "R2c",
        "LL",
        "Equation",
        "shape",
        "max"
      )
    Results.cont <- Results.cont[order(Results.cont$AICc), ]
    write.table(
      Results.cont,
      paste0(GA.inputs$Results.dir, "ContinuousResults.csv"),
      sep = ",",
      col.names = T,
      row.names = F
    )
  }
  
  # Full Results
  if (nrow(Results.cat) > 0 & nrow(Results.cont) > 0) {
    Results.All <- rbind(Results.cat[, c(1:8)], Results.cont[, c(1:8)])
  } else if (nrow(Results.cat) < 1 & nrow(Results.cont) > 0) {
    Results.All <- (Results.cont[, c(1:8)])
  } else {
    Results.All <- (Results.cat[, c(1:8)])
  }
  
  if (dist_mod == TRUE)
    Results.All <- rbind(Results.All, Dist.AICc)
  if (null_mod == TRUE)
    Results.All <- rbind(Results.All, Null.AICc)
  
  Results.All <- Results.All[order(Results.All$AICc), ]
  
  cat("\n")
  cat("\n")
  write.table(
    Results.All,
    paste0(GA.inputs$Results.dir, "All_Results_Table.csv"),
    
    sep = ",",
    col.names = T,
    row.names = F
  )
  
  # if(!is.null(gdist.inputs$covariates)) { 
  #   MLPE.results <- NULL
  # } else {
  # Get parameter estimates
  if (!is.null(CS.inputs)) {
    MLPE.results <- MLPE.lmm_coef(
      resistance = GA.inputs$Results.dir,
      genetic.dist = CS.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "cs",
      ID = CS.inputs$ID,
      ZZ = CS.inputs$ZZ
    )
    
  } else if(!is.null(jl.inputs)) {
    MLPE.results <- MLPE.lmm_coef(
      formula = jl.inputs$formula,
      inputs = jl.inputs$df,
      resistance = GA.inputs$Results.dir,
      genetic.dist = jl.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "jl",
      ID = jl.inputs$ID,
      ZZ = jl.inputs$ZZ
    )
  } else {
    MLPE.results <- MLPE.lmm_coef(
      resistance = GA.inputs$Results.dir,
      genetic.dist = gdist.inputs$response,
      out.dir = GA.inputs$Results.dir,
      method = "gd",
      ID = gdist.inputs$ID,
      ZZ = gdist.inputs$ZZ
    )
  } 
  # } ## End covariate if-else
  
  k.list <- plyr::ldply(k.list)
  colnames(k.list) <- c("surface", "k")
  
  rt <- proc.time()[3] - t1
  # Full Results
  if (nrow(Results.cat) > 0 & nrow(Results.cont) > 0) {
    RESULTS <-
      list(
        ContinuousResults = Results.cont,
        CategoricalResults = Results.cat,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list,
        ga = ga.list
      )
    
  } else if (nrow(Results.cat) < 1 & nrow(Results.cont) > 0) {
    RESULTS <-
      list(
        ContinuousResults = Results.cont,
        CategoricalResults = NULL,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list,
        ga = ga.list
      )
    
  } else if (nrow(Results.cat) > 0 & nrow(Results.cont) < 1) {
    RESULTS <-
      list(
        ContinuousResults = NULL,
        CategoricalResults = Results.cat,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list,
        ga = ga.list
      )
  } else {
    RESULTS <-
      list(
        ContinuousResults = NULL,
        CategoricalResults = NULL,
        AICc = Results.All,
        MLPE = MLPE.results,
        Run.Time = rt,
        MLPE.list = MLPE.list,
        cd = cd.list,
        k = k.list,
        ga = ga.list
      )
  }
  # rm(single.GA, r)
  setwd(wd)
  gc()
  return(RESULTS)
}