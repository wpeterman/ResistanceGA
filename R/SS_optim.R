#' Single surface optimization
#'
#' Optimize all surfaces contained in a directory using a genetic algorithm executed with the \code{\link[GA]{ga}} function in the Genetic Algorithms package \pkg{GA}
#'
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param nlm Logical, if TRUE, the final step of optimization will use nlm to fine-tune parameter estimates. This may lead to overfitting in some cases. Default = FALSE.
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
#' @usage SS_optim(CS.inputs, gdist.inputs, GA.inputs, nlm, dist_mod, null_mod)
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @export
SS_optim <- function(CS.inputs = NULL,
                     gdist.inputs = NULL,
                     GA.inputs,
                     nlm = FALSE,
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

  # Optimize each surface in turn
  for (i in 1:GA.inputs$n.layers) {
    r <- GA.inputs$Resistance.stack[[i]]
    names(r) <- GA.inputs$layer.names[i]

    # CIRCUITSCAPE ------------------------------------------------------------
    
    
    # * Categorical -----------------------------------------------------------
    
    # Processing of categorical surfaces
    if (!is.null(CS.inputs)) {
      if (GA.inputs$parallel != FALSE) {
        warning(
          "\n CIRCUITSCAPE cannot be optimized in parallel. \n Ignoring parallel arguement. \n If you want to optimize in parallel, use least cost paths and gdistance.",
          immediate. = TRUE
        )
      }
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
          popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]),
          maxiter = GA.inputs$maxiter,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          elitism = GA.inputs$percent.elite,
          mutation = GA.inputs$mutation,
          seed = GA.inputs$seed,
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
        names(r) <- GA.inputs$layer.names[i]
        
        Run_CS(CS.inputs, GA.inputs, r, EXPORT.dir = GA.inputs$Results.dir)
        
        Diagnostic.Plots(
          resistance.mat = paste0(
            GA.inputs$Results.dir,
            GA.inputs$layer.names[i],
            "_resistances.out"
          ),
          genetic.dist = CS.inputs$response,
          plot.dir = GA.inputs$Plots.dir,
          type = "categorical",
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
        
        fit.stats <-
          r.squaredGLMM(
            MLPE.lmm(
              resistance = paste0(
                GA.inputs$Results.dir,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              pairwise.genetic = CS.inputs$response,
              REML = F,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        aic <-
          AIC(
            MLPE.lmm(
              resistance = paste0(
                GA.inputs$Results.dir,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              pairwise.genetic = CS.inputs$response,
              REML = F,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        LL <-
          logLik(
            MLPE.lmm(
              resistance = paste0(
                GA.inputs$Results.dir,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
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
          (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
        # AICc <- (aic)+(((2*k)*(k+1))/(CS.inputs$n.Pops-k-1))
        
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
          resistance = paste0(
            GA.inputs$Results.dir,
            GA.inputs$layer.names[i],
            "_resistances.out"
          ),
          pairwise.genetic = CS.inputs$response,
          REML = TRUE,
          ID = CS.inputs$ID,
          ZZ = CS.inputs$ZZ
        )
        
        cd.list[[i]] <- (read.table(paste0(
          GA.inputs$Results.dir,
          GA.inputs$layer.names[i],
          "_resistances.out"))[-1, -1])
        
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
          popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]),
          maxiter = GA.inputs$maxiter,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          elitism = GA.inputs$percent.elite,
          mutation = GA.inputs$mutation,
          seed = GA.inputs$seed,
          iter = i,
          quiet = GA.inputs$quiet
        )
        
        if(dim(single.GA@solution)[1] > 1) {
          single.GA@solution <- t(as.matrix(single.GA@solution[1,]))
        }
        
        # Using GA results, optimize with nlm
        start.vals <- single.GA@solution[-1]
        
        if (nlm == TRUE) {
          # Informed start values; these are the optimized values from the single parameter optimization
          EQ <- get.EQ(single.GA@solution[1])
          Optim.nlm <-
            nlm(
              Resistance.Optimization_cont.nlm,
              log(start.vals),
              Resistance = r,
              equation = single.GA@solution[1],
              get.best = FALSE,
              CS.inputs = CS.inputs,
              Min.Max = 'min',
              write.dir = GA.inputs$Write.dir
            )
          
          OPTIM <-
            Resistance.Optimization_cont.nlm(
              PARM = (Optim.nlm$estimate),
              Resistance = r,
              equation = single.GA@solution[1],
              get.best = TRUE,
              CS.inputs = CS.inputs,
              Min.Max = 'min',
              write.dir = GA.inputs$Results.dir
            )
          
          Diagnostic.Plots(
            resistance.mat = paste0(
              GA.inputs$Results.dir,
              GA.inputs$layer.names[i],
              "_resistances.out"
            ),
            genetic.dist = CS.inputs$response,
            plot.dir = GA.inputs$Plots.dir,
            type = "continuous",
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
          
          Plot.trans(
            PARM = exp(Optim.nlm$estimate),
            Resistance = GA.inputs$Resistance.stack[[i]],
            transformation = EQ,
            print.dir = GA.inputs$Plots.dir,
            Name = GA.inputs$layer.names[i]
          )
          
          RS <-
            data.frame(GA.inputs$layer.names[i],
                       Optim.nlm$minimum,
                       EQ,
                       Cont.Param(exp(Optim.nlm$estimate)))
          colnames(RS) <-
            c("Surface",
              "obj.func",
              "R2m",
              "R2c",
              "Equation",
              "shape",
              "max")
          RESULTS.cont[[cnt2]] <- RS
          
          MLPE.list[[i]] <- MLPE.lmm(
            resistance = paste0(
              GA.inputs$Results.dir,
              GA.inputs$layer.names[i],
              "_resistances.out"
            ),
            pairwise.genetic = CS.inputs$response,
            REML = TRUE,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
          
          cd.list[[i]] <- (read.table(paste0(
            GA.inputs$Results.dir,
            GA.inputs$layer.names[i],
            "_resistances.out"))[-1, -1])
          
          names(MLPE.list)[i] <- GA.inputs$layer.names[i]
          names(cd.list)[i] <- GA.inputs$layer.names[i]          

        } else {
          
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
          names(r.tran) <- GA.inputs$layer.names[i]
          
          Run_CS(CS.inputs,
                 GA.inputs,
                 r.tran,
                 EXPORT.dir = GA.inputs$Results.dir)
          
          Diagnostic.Plots(
            resistance.mat = paste0(
              GA.inputs$Results.dir,
              GA.inputs$layer.names[i],
              "_resistances.out"
            ),
            genetic.dist = CS.inputs$response,
            plot.dir = GA.inputs$Plots.dir,
            type = "continuous",
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
              resistance = paste0(
                GA.inputs$Results.dir,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              pairwise.genetic = CS.inputs$response,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
          aic <- AIC(
            MLPE.lmm(
              REML = F,
              resistance = paste0(
                GA.inputs$Results.dir,
                GA.inputs$layer.names[i],
                "_resistances.out"
              ),
              pairwise.genetic = CS.inputs$response,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
          
          LL <-
            logLik(
              MLPE.lmm(
                resistance = paste0(
                  GA.inputs$Results.dir,
                  GA.inputs$layer.names[i],
                  "_resistances.out"
                ),
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
            (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
          # AICc <- (aic)+(((2*k)*(k+1))/(CS.inputs$n.Pops-k-1))
          
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
            resistance = paste0(
              GA.inputs$Results.dir,
              GA.inputs$layer.names[i],
              "_resistances.out"
            ),
            pairwise.genetic = CS.inputs$response,
            REML = TRUE,
            ID = CS.inputs$ID,
            ZZ = CS.inputs$ZZ
          )
          
          cd.list[[i]] <- (read.table(paste0(
            GA.inputs$Results.dir,
            GA.inputs$layer.names[i],
            "_resistances.out"))[-1, -1])
          
          names(MLPE.list)[i] <- GA.inputs$layer.names[i]
          names(cd.list)[i] <- GA.inputs$layer.names[i]
          

        }
      } # Close if-else
      if (dist_mod == TRUE) {
        r <- reclassify(r, c(-Inf, Inf, 1))
        names(r) <- "dist"
        cd <- Run_CS(CS.inputs, GA.inputs, r, full.mat = T)

        Dist.AIC <-
          AIC(
            MLPE.lmm(
              resistance = paste0(GA.inputs$Write.dir, "dist_resistances.out"),
              pairwise.genetic = CS.inputs$response,
              REML = FALSE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        fit.stats <-
          r.squaredGLMM(
            MLPE.lmm(
              resistance = paste0(GA.inputs$Write.dir, "dist_resistances.out"),
              pairwise.genetic = CS.inputs$response,
              REML = FALSE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        LL <-
          logLik(
            MLPE.lmm(
              resistance = paste0(GA.inputs$Write.dir, "dist_resistances.out"),
              pairwise.genetic = CS.inputs$response,
              REML = FALSE,
              ID = CS.inputs$ID,
              ZZ = CS.inputs$ZZ
            )
          )
        
        MLPE.list[[i + 1]] <- MLPE.lmm(
          resistance = paste0(GA.inputs$Write.dir, "dist_resistances.out"),
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
          (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
        # (Dist.AIC)+(((2*k)*(k+1))/((CS.inputs$n.Pops)-k-1))
        
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
        Null.AIC <-
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
          (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
        # AICc <- (Null.AIC)+(((2*k)*(k+1))/((CS.inputs$n.Pops)-k-1))
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
    

    # Optimize with gdistance -------------------------------------------------
    
    if (!is.null(gdist.inputs)) {
      # * Categorical -----------------------------------------------------------

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
          popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]),
          maxiter = GA.inputs$maxiter,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          elitism = GA.inputs$percent.elite,
          mutation = GA.inputs$mutation,
          seed = GA.inputs$seed,
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
        
        cd <- Run_gdistance(gdist.inputs, r)
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
        
        Diagnostic.Plots(
          resistance.mat = cd,
          genetic.dist = gdist.inputs$response,
          plot.dir = GA.inputs$Plots.dir,
          type = "categorical",
          name = NAME,
          ID = gdist.inputs$ID,
          ZZ = gdist.inputs$ZZ
        )
        
        fit.stats <- r.squaredGLMM(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            REML = F,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ
          )
        )
        
        aic <- AIC(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            REML = F,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ
          )
        )
        
        LL <- logLik(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            REML = F,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ
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

        n <- gdist.inputs$n.Pops
        AICc <-
          (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
        # AICc <- (aic)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
        
        
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
          gdist.inputs = gdist.inputs,
          lower = GA.inputs$min.list[[i]],
          upper = GA.inputs$max.list[[i]],
          parallel = GA.inputs$parallel,
          popSize = GA.inputs$pop.mult * length(GA.inputs$max.list[[i]]),
          maxiter = GA.inputs$maxiter,
          run = GA.inputs$run,
          keepBest = GA.inputs$keepBest,
          elitism = GA.inputs$percent.elite,
          mutation = GA.inputs$mutation,
          seed = GA.inputs$seed,
          iter = i,
          quiet = GA.inputs$quiet
        )

# ** Second optim ---------------------------------------------------------

        
        # Using GA results, optimize with nlm
        start.vals <- single.GA@solution[-1]
        
        if (nlm == TRUE) {
          # Informed start values; these are the optimized values from the single parameter optimization
          EQ <- get.EQ(single.GA@solution[1])
          Optim.nlm <-
            nlm(
              Resistance.Optimization_cont.nlm,
              log(start.vals),
              Resistance = r,
              equation = single.GA@solution[1],
              get.best = FALSE,
              gdist.inputs = gdist.inputs,
              Min.Max = 'min',
              write.dir = GA.inputs$Write.dir
            )
          
          r <-
            Resistance.tran(
              transformation = EQ,
              shape = Optim.nlm$estimate[1],
              max = Optim.nlm$estimate[2],
              r = r,
              out = GA.inputs$Results.dir
            )
          
          names(r) <- GA.inputs$layer.names[i]
          NAME <- GA.inputs$layer.names[i]
          
          cd <- Run_gdistance(gdist.inputs, r)
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
            PARM = exp(Optim.nlm$estimate),
            Resistance = GA.inputs$Resistance.stack[[i]],
            transformation = EQ,
            print.dir = GA.inputs$Plots.dir,
            Name = GA.inputs$layer.names[i]
          )
          
          fit.stats <-
            r.squaredGLMM(
              MLPE.lmm(
                resistance = cd,
                pairwise.genetic = gdist.inputs$response,
                REML = F,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
            )
          
          
          aic <-
            AIC(
              MLPE.lmm(
                resistance = cd,
                pairwise.genetic = gdist.inputs$response,
                REML = F,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
            )
          
          RS <- data.frame(
            GA.inputs$layer.names[i],-single.GA@fitnessValue,
            fit.stats[[1]],
            fit.stats[[2]],
            single.GA@solution
          )
          
          # RS<-data.frame(GA.inputs$layer.names[i],Optim.nlm$minimum,EQ,Cont.Param(exp(Optim.nlm$estimate)))
          colnames(RS) <-
            c("Surface",
              "AICc",
              "R2m",
              "R2c",
              "Equation",
              "shape",
              "max")
          RESULTS.cont[[cnt2]] <- RS
          
          # * Continuous -----------------------------------------------------------

        } else {
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
          
          fit.stats <-
            r.squaredGLMM(
              MLPE.lmm(
                resistance = cd,
                pairwise.genetic = gdist.inputs$response,
                REML = F,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
            )
          
          
          aic <-
            AIC(
              MLPE.lmm(
                resistance = cd,
                pairwise.genetic = gdist.inputs$response,
                REML = F,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
            )
          
          LL <-
            logLik(
              MLPE.lmm(
                resistance = cd,
                pairwise.genetic = gdist.inputs$response,
                REML = F,
                ID = gdist.inputs$ID,
                ZZ = gdist.inputs$ZZ
              )
            )
          
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
            (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
          # AICc <- (aic)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
          
          
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
        }
      } # Close if-else
      
      if (dist_mod == TRUE) {
        r <- reclassify(r, c(-Inf, Inf, 1))
        names(r) <- "dist"
        cd <- Run_gdistance(gdist.inputs, r)
        
        write.table(
          as.matrix(cd),
          file = paste0(GA.inputs$Results.dir, 'Distance', "_", gdist.inputs$method, "_distMat.csv"),
          sep = ",",
          row.names = F,
          col.names = F
        )
        
        Dist.AIC <- suppressWarnings(AIC(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ,
            REML = FALSE
          )
        ))
        
        fit.stats <- r.squaredGLMM(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ,
            REML = FALSE
          )
        )
        
        LL <- logLik(
          MLPE.lmm2(
            resistance = cd,
            response = gdist.inputs$response,
            ID = gdist.inputs$ID,
            ZZ = gdist.inputs$ZZ,
            REML = FALSE
          )
        )
        
        MLPE.list[[i + 1]] <-  MLPE.lmm2(
          resistance = cd,
          response = gdist.inputs$response,
          REML = TRUE,
          ID = gdist.inputs$ID,
          ZZ = gdist.inputs$ZZ
        )
        
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
          (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
        # AICc <- (Dist.AIC)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
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
        Null.AIC <-
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
          (-2 * LL) + (2 * k) + (((2 * k) * (k + 1)) / (n - k - 1))
        # AICc <- (Null.AIC)+(((2*k)*(k+1))/(gdist.inputs$n.Pops-k-1))
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
    }
  } # Close ascii loop
  
  ####################################################
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
  ##################################
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
    paste0(GA.inputs$Results.dir, "All_Results_Table_", gdist.inputs$method,".csv"),

    sep = ",",
    col.names = T,
    row.names = F
  )
  
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
        k = k.list
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
        k = k.list
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
        k = k.list
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
        k = k.list
      )
  }
  
  unlink(GA.inputs$Write.dir, recursive = T, force = T)
  return(RESULTS)
  ###############################################################################################################
}