##################################
#' Create diagnostic plots
#'
#' This function will generate mixed effect model diagnostic plots following optimization
#'
#' @param resistance.mat Path to CIRCUITSCAPE "_resistances.out" file, or costDistance object created from running gdistance
#' @param genetic.dist Vector of pairwise genetic distances (lower half of pairwise matrix). Can be input as CS.inputs$response
#' @param XLAB Label for x-axis (Defaults to "Estimated resistance")
#' @param YLAB Label for y-axis (Defaults to "Genetic distance")
#' @param plot.dir Directory to output TIFF of diagnostic plots
#' @param type Specify whether the optimized surface is "continuous" or "categorical"
#' @param name The name to be attached to the output file. Must be specified when getting diagnostic plots for gdistance models
#' @param ID The to_from ID list for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @param ZZ The sparse matrix object for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @return A multipanel panel .tif including histogram of residuals and qqplot of fitted mixed effects model

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage Diagnostic.Plots(resistance.mat, genetic.dist, XLAB,YLAB, plot.dir, type, name, ID, ZZ)

Diagnostic.Plots <-
  function(resistance.mat,
           genetic.dist,
           XLAB = "Estimated resistance",
           YLAB = "Genetic distance",
           plot.dir,
           type = "categorical",
           name = NULL,
           ID = NULL,
           ZZ = NULL) {
    if (length(resistance.mat) == 1) {
      response = genetic.dist
      if (is.null(name)) {
        NAME <-
          gsub(pattern = "*_resistances.out", "", x = (basename(resistance.mat)))
      }
      mm <- read.table(resistance.mat)[-1, -1]
      m <- length(mm)
      mm <- lower(mm)
      mm <- mm[which(mm != -1)]
      
      if (is.null(ID)) {
        ID <- To.From.ID(POPS = m)
      }
      if (is.null(ZZ)) {
        ZZ <- ZZ.mat(ID = ID)
      }
      
      cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
      cs.unscale <- mm
      dat <- data.frame(ID, cs.matrix = cs.matrix, response = response)
      colnames(dat) <- c("pop1", "pop2", "cs.matrix", "response")
      
      # Assign value to layer
      LAYER <- assign("LAYER", value = dat$cs.matrix)
      
      # Fit model
      mod <- lFormula(response ~ LAYER + (1 | pop1),
                      data = dat,
                      REML = TRUE)
      mod$reTrms$Zt <- ZZ
      dfun <- do.call(mkLmerDevfun, mod)
      opt <- optimizeLmer(dfun)
      Mod <-
        (mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    }
    
    
    
    if (length(resistance.mat) > 1) {
      response = genetic.dist
      if (is.null(name)) {
        stop("Output file 'name' must be specified!!!")
      }
      NAME <- name
      mm <- lower(as.matrix(resistance.mat))
      m <- attr(resistance.mat, "Size")
      mm <- mm[which(mm != -1)]
      
      if (is.null(ID)) {
        ID <- To.From.ID(POPS = m)
      }
      if (is.null(ZZ)) {
        ZZ <- ZZ.mat(ID = ID)
      }
      cs.matrix <- scale(mm, center = TRUE, scale = TRUE)
      cs.unscale <- mm
      dat <- data.frame(ID, cs.matrix = cs.matrix, response = response)
      colnames(dat) <- c("pop1", "pop2", "cs.matrix", "response")
      
      # Assign value to layer
      LAYER <- assign("LAYER", value = dat$cs.matrix)
      
      # Fit model
      mod <- lFormula(response ~ LAYER + (1 | pop1),
                      data = dat,
                      REML = TRUE)
      mod$reTrms$Zt <- ZZ
      dfun <- do.call(mkLmerDevfun, mod)
      opt <- optimizeLmer(dfun)
      Mod <-
        (mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    }
    #######
    # Make diagnostic plots
    if (type == "categorical") {
      tiff(
        filename = paste0(plot.dir, NAME, "_DiagnosticPlots.tif"),
        width = 279,
        height = 215,
        units = "mm",
        compression = c("lzw"),
        bg = "white",
        res = 300
      )
      par(
        mfrow = c(2, 1),
        oma = c(0, 4, 0, 0) + 0.1,
        mar = c(4, 4, 1, 1) + 0.1
      )
      hist(residuals(Mod), xlab = "Residuals", main = "")
      qqnorm(resid(Mod), main = "")
      qqline(resid(Mod))
      dev.off()
      par(mfrow = c(1, 1))
    } else {
      tiff(
        filename = paste0(plot.dir, NAME, "_DiagnosticPlots.tif"),
        width = 279,
        height = 215,
        units = "mm",
        compression = c("lzw"),
        bg = "white",
        res = 300
      )
      par(
        mfrow = c(2, 2),
        oma = c(0, 4, 0, 0) + 0.1,
        mar = c(4, 4, 1, 1) + 0.1
      )
      plot(dat$response ~ cs.unscale,
           xlab = XLAB,
           ylab = YLAB)
      abline(lm(dat$response ~ cs.unscale))
      plot(residuals(Mod) ~ cs.unscale,
           xlab = XLAB,
           ylab = "Residuals")
      abline(lm(residuals(Mod) ~ cs.unscale))
      hist(residuals(Mod), xlab = "Residuals", main = "")
      qqnorm(resid(Mod), main = "")
      qqline(resid(Mod))
      dev.off()
      par(mfrow = c(1, 1))
    }
  }