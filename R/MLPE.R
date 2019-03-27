# Run Mixed effects models, recover parameter estimates
#' Run maximum likelihood population effects mixed effects model (MLPE)
#'
#' Runs MLPE as detailed by Clarke et al. (2002). This function will run the model and return lmer object
#'
#' @param resistance Path to pairwise resistance distance matrix (resistances.out) from CS results. Alternatively, provide the pairwise resistances created from optimizing with `gdistance` (result of Run_gdistance).
#' @param pairwise.genetic Lower half of pairwise genetic distance matrix
#' @param REML Logical. If TRUE, mixed effects model will be fit using restricted maximum likelihood. Default = FALSE
#' @param ID The to_from ID list for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @param ZZ The sparse matrix object for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @param scale Specify whether the pairwise distance values be scaled and centered (Default = TRUE)
#' @return A lmer object from the fitted model
#' @details An AIC value will only be returned if \code{REML = FALSE}

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage MLPE.lmm(resistance, 
#'                 pairwise.genetic, 
#'                 REML = FALSE, 
#'                 ID = NULL, 
#'                 ZZ = NULL,
#'                 scale = TRUE)
#' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

MLPE.lmm <-
  function(resistance,
           pairwise.genetic,
           REML = FALSE,
           ID = NULL,
           ZZ = NULL,
           scale = TRUE) {
    response = pairwise.genetic
    
    if (class(resistance)[[1]] == 'dist') {
      mm <- as.vector(resistance)
      m <- attr(resistance, "Size")
      mm <- mm[which(mm != -1)]
      
      if (is.null(ID)) {
        ID <- To.From.ID(POPS = m)
      }
      if (is.null(ZZ)) {
        ZZ <- ZZ.mat(ID = ID)
      }
      
      if(scale == T) {
        cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
      } else {
        cs.matrix <- mm
      }
      
    } else if(!is.character(resistance)) {
      if(is.vector(resistance)) {
        mm <- resistance
        m <- 0.5 * (sqrt((8 * length(mm)) + 1) + 1)
        mm <- mm[which(mm != -1)]
      } else {
        mm <- resistance
        m <- nrow(mm)
        mm <- lower(mm)
        mm <- mm[which(mm != -1)]
      }
      
      if (is.null(ID)) {
        ID <- To.From.ID(POPS = m)
      }
      if (is.null(ZZ)) {
        ZZ <- ZZ.mat(ID = ID)
      }
      
      if(scale == T) {
        cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
      } else {
        cs.matrix <- mm
      }  
      
    } else {
      mm <- (read.table(resistance)[-1, -1])
      m <- nrow(mm)
      mm <- lower(mm)
      mm <- mm[which(mm != -1)]
      
      if (is.null(ID)) {
        ID <- To.From.ID(POPS = m)
      }
      if (is.null(ZZ)) {
        ZZ <- ZZ.mat(ID = ID)
      }
      
      if(scale == T) {
        cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
      } else {
        cs.matrix <- mm
      }
    }
    
    dat <- data.frame(ID, resistance = cs.matrix, response = response)
    colnames(dat) <- c("pop1", "pop2", "resistance", "response")
    
    # Assign value to layer
    #     LAYER<-assign("Resist",value=dat$cs.matrix)
    
    # Fit model
    mod <-
      lFormula(response ~ resistance + (1 | pop1),
               data = dat,
               REML = REML)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun, mod)
    opt <- optimizeLmer(dfun)
    MOD <-
      (mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    return(MOD)
  }


MLPE.lmm2 <- function(resistance, response, REML = FALSE, ID, ZZ) {
  if (class(resistance)[1] != 'dist') {
    res <- resistance[which(resistance != -1)]
    
    dat <- data.frame(ID, resistance = res, response = response)
    colnames(dat) <- c("pop1", "pop2", "resistance", "response")
  } else {
    resistance <- as.vector(resistance)
    resistance <- resistance[which(resistance != -1)]
    
    dat <- data.frame(ID, resistance = resistance, response = response)
    colnames(dat) <- c("pop1", "pop2", "resistance", "response")
    
  }
  # Assign value to layer
  #   LAYER<-assign("Resist",value=dat$resistance)
  
  # Fit model
  mod <-
    lFormula(response ~ scale(resistance) + (1 | pop1),
             data = dat,
             REML = REML)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun, mod)
  opt <- optimizeLmer(dfun)
  MOD <-
    (mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
  return(MOD)
}


MLPE.lmm_coef <-
  function(resistance,
           genetic.dist,
           out.dir = NULL,
           method,
           dat = NULL,
           formula = NULL,
           inputs = NULL,
           ID = NULL,
           ZZ = NULL) {
    
    if(!is.null(formula)) {
      if (method == "cs") {
        response = genetic.dist
        resist.mat <-
          list.files(resistance, pattern = "*_csResistMat.csv", full.names = TRUE)
        resist.names <-
          gsub(pattern = "_csResistMat.csv",
               "",
               x = list.files(resistance, pattern = "*_csResistMat.csv"))
        COEF.Table <- array()
        for (i in 1:length(resist.mat)) {
          cd <- read.csv(resist.mat[i], header = F)
          cd.l <- lower(cd)
          cd.l <- cd.l[cd.l != -1]

          scale.cd <- scale(cd.l, center = TRUE, scale = TRUE)
          cs.unscale <- cd.l

          # Assign value to layer
          LAYER <- assign(resist.names[i], value = dat$cs.matrix)
          
          # Fit model
          mod <- mlpe_rga(formula = formula,
                          data = inputs$df,
                          ZZ = inputs$ZZ,
                          REML = TRUE)
          
          Mod.Summary <- summary(mod)
          COEF <- Mod.Summary$coefficients
          row.names(COEF) <- c("Intercept", resist.names[i])
          COEF.Table <- rbind(COEF.Table, COEF)
        }
      } else if(method == "jl") {
        response <- genetic.dist
        resist.mat <-
          list.files(resistance, pattern = "*_jlResistMat.csv", full.names = TRUE)
        resist.names <-
          gsub(pattern = "_jlResistMat.csv",
               "",
               x = list.files(resistance, pattern = "*_jlResistMat.csv"))
        
        
        COEF.Table <- array()
        for (i in 1:length(resist.mat)) {
          cd.mat <- read.csv(resist.mat[i], header = F)
          cd.l <- lower(cd.mat)
          cd.l <- cd.l[cd.l != -1]
          
          scale.cd <- scale(cd.l, center = TRUE, scale = TRUE)
          cs.unscale <- cd.l
          
          # scale.cd <- scale(cd.l[which(cd.l != -1)], center = TRUE, scale = TRUE)
          # cs.unscale <- cd.l[which(cd.l != -1)]
          
          dat <- inputs
          dat$cd <- scale.cd
          
          # # Assign value to layer
          # cd <- assign(resist.names[i], value = inputs$cd)
          
          # Fit model
          mod <- mlpe_rga(formula = formula,
                          data = dat,
                          REML = TRUE,
                          ZZ = ZZ)
            
            
          #   lFormula(formula,
          #            data = inputs,
          #            REML = TRUE)
          # mod$reTrms$Zt <- ZZ
          # dfun <- do.call(mkLmerDevfun, mod)
          # opt <- optimizeLmer(dfun)
          Mod.Summary <-
            summary(mod)
          COEF <- Mod.Summary$coefficients
          if(nrow(COEF) > 2) {
            row.names(COEF)[length(row.names(COEF))] <- resist.names[i]
            
          } else {
            row.names(COEF) <- c("Intercept", resist.names[i])
            
          }
          COEF.Table <- rbind(COEF.Table, COEF)
        }
      } else {
        response <- genetic.dist
        resist.mat <-
          list.files(resistance, pattern = "*_distMat.csv", full.names = TRUE)
        resist.names <-
          gsub(pattern = "_distMat.csv",
               "",
               x = list.files(resistance, pattern = "*_distMat.csv"))
        
        if(length(agrep("_commuteDistance", resist.names)) > 0){
          resist.names <- plyr::ldply(strsplit(resist.names, "_commuteDistance"))[,1]
          
        } else {
          resist.names <- plyr::ldply(strsplit(resist.names, "_costDistance"))[,1]
        }
        
        
        COEF.Table <- array()
        for (i in 1:length(resist.mat)) {
          cd <- read.csv(resist.mat[i], header = F)
          mm <- lower(cd)
          # mm <- lower(cd)
          m <- dim(cd)[1]
          ID <- To.From.ID(POPS = m)
          ZZ <- ZZ.mat(ID = ID)
          cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
          cs.unscale <- mm
          dat <- cbind(ID, cs.matrix, response)
          
          # Assign value to layer
          LAYER <- assign(resist.names[i], value = dat$cs.matrix)
          
          # Fit model
          mod <-
            lFormula(response ~ LAYER + (1 | pop1),
                     data = dat,
                     REML = TRUE)
          mod$reTrms$Zt <- ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          Mod.Summary <-
            summary(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          COEF <- Mod.Summary$coefficients
          row.names(COEF) <- c("Intercept", resist.names[i])
          COEF.Table <- rbind(COEF.Table, COEF)
        }
      }
    } else {
      if (method == "cs") {
        response = genetic.dist
        resist.mat <-
          list.files(resistance, pattern = "*_csResistMat.csv", full.names = TRUE)
        resist.names <-
          gsub(pattern = "_csResistMat.csv",
               "",
               x = list.files(resistance, pattern = "*_csResistMat.csv"))
        COEF.Table <- array()
        for (i in 1:length(resist.mat)) {
          cd <- read.csv(resist.mat[i], header = F)
          mm <- lower(cd)
          # mm <- lower(cd)
          m <- dim(cd)[1]
          ID <- To.From.ID(POPS = m)
          ZZ <- ZZ.mat(ID = ID)
          cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
          cs.unscale <- mm
          dat <- cbind(ID, cs.matrix, response)
          
          # Assign value to layer
          LAYER <- assign(resist.names[i], value = dat$cs.matrix)
          
          # Fit model
          mod <-
            lFormula(response ~ LAYER + (1 | pop1),
                     data = dat,
                     REML = TRUE)
          mod$reTrms$Zt <- ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          Mod.Summary <-
            summary(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          COEF <- Mod.Summary$coefficients
          row.names(COEF) <- c("Intercept", resist.names[i])
          COEF.Table <- rbind(COEF.Table, COEF)
        }
      } else if(method == "jl") {
        response <- genetic.dist
        resist.mat <-
          list.files(resistance, pattern = "*_jlResistMat.csv", full.names = TRUE)
        resist.names <-
          gsub(pattern = "_jlResistMat.csv",
               "",
               x = list.files(resistance, pattern = "*_jlResistMat.csv"))
        
        
        COEF.Table <- array()
        for (i in 1:length(resist.mat)) {
          cd <- read.csv(resist.mat[i], header = F)
          mm <- lower(cd)
          # mm <- lower(cd)
          m <- dim(cd)[1]
          ID <- To.From.ID(POPS = m)
          ZZ <- ZZ.mat(ID = ID)
          cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
          cs.unscale <- mm
          dat <- cbind(ID, cs.matrix, response)
          
          # Assign value to layer
          LAYER <- assign(resist.names[i], value = dat$cs.matrix)
          
          # Fit model
          mod <-
            lFormula(response ~ LAYER + (1 | pop1),
                     data = dat,
                     REML = TRUE)
          mod$reTrms$Zt <- ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          Mod.Summary <-
            summary(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          COEF <- Mod.Summary$coefficients
          row.names(COEF) <- c("Intercept", resist.names[i])
          COEF.Table <- rbind(COEF.Table, COEF)
        }
      } else {
        response <- genetic.dist
        resist.mat <-
          list.files(resistance, pattern = "*_distMat.csv", full.names = TRUE)
        resist.names <-
          gsub(pattern = "_distMat.csv",
               "",
               x = list.files(resistance, pattern = "*_distMat.csv"))
        
        if(length(agrep("_commuteDistance", resist.names)) > 0){
          resist.names <- plyr::ldply(strsplit(resist.names, "_commuteDistance"))[,1]
          
        } else {
          resist.names <- plyr::ldply(strsplit(resist.names, "_costDistance"))[,1]
        }
        
        
        COEF.Table <- array()
        for (i in 1:length(resist.mat)) {
          cd <- read.csv(resist.mat[i], header = F)
          mm <- lower(cd)
          # mm <- lower(cd)
          m <- dim(cd)[1]
          ID <- To.From.ID(POPS = m)
          ZZ <- ZZ.mat(ID = ID)
          cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
          cs.unscale <- mm
          dat <- cbind(ID, cs.matrix, response)
          
          # Assign value to layer
          LAYER <- assign(resist.names[i], value = dat$cs.matrix)
          
          # Fit model
          mod <-
            lFormula(response ~ LAYER + (1 | pop1),
                     data = dat,
                     REML = TRUE)
          mod$reTrms$Zt <- ZZ
          dfun <- do.call(mkLmerDevfun, mod)
          opt <- optimizeLmer(dfun)
          Mod.Summary <-
            summary(mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
          COEF <- Mod.Summary$coefficients
          row.names(COEF) <- c("Intercept", resist.names[i])
          COEF.Table <- rbind(COEF.Table, COEF)
        }
      }
    }
    
    
    
    
    if (is.null(out.dir)) {
      COEF.Table <- (COEF.Table[-1, ])
    } else {
      COEF.Table <- COEF.Table[-1, ]
      write.table(
        COEF.Table,
        file = paste0(out.dir, "MLPE_coeff_Table.csv"),
        sep = ",",
        row.names = T,
        col.names = NA
      )
      return(COEF.Table)
    }
  }