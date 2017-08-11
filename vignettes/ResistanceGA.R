## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo = FALSE-------------------------------------------------------
library(ggplot2, warn.conflicts = F, quietly = T)
library(raster, warn.conflicts = F, quietly = T)

## ----install.package, eval=FALSE-----------------------------------------
#  # Install 'devtools' package, if needed
#  if(!("devtools" %in% list.files(.libPaths()))) {
#      install.packages("devtools", repo = "http://cran.rstudio.com", dep = TRUE)
#  }
#  
#  devtools::install_github("wpeterman/ResistanceGA", build_vignettes = TRUE) # Download package

## ----results='hide',message=FALSE, warning=FALSE-------------------------
library(ResistanceGA)
rm(list = ls())

## ----Plot.trans.demo, fig.height = 4, fig.width = 4----------------------
Ricker.plot <- Plot.trans(PARM = c(1.5, 200),
                          Resistance = c(1,10),
                          transformation="Ricker")

# Change title of plot
Ricker.plot$labels$title <- "Ricker Tansformation"
Ricker.plot

# Find original data value that now has maximum resistance
Ricker.plot$data$original[which(Ricker.plot$data$transformed==max(Ricker.plot$data$transformed))]

## ----warning=FALSE, results='hide',message=FALSE, eval=FALSE-------------
#  
#  if("ResistanceGA_Examples"%in%dir("C:/")==FALSE)
#    dir.create(file.path("C:/", "ResistanceGA_Examples"))
#  
#  # Create a subdirectory for the first example
#  dir.create(file.path("C:/ResistanceGA_Examples/","SingleSurface"))
#  
#  write.dir <- "C:/ResistanceGA_Examples/SingleSurface/"      # Directory to write .asc files and results
#  
#  # Give path to CIRCUITSCAPE .exe file
#  # Default = '"C:/Program Files/Circuitscape/cs_run.exe"'
#  CS.program <- paste('"C:/Program Files/Circuitscape/cs_run.exe"')

## ----load.data, echo = FALSE, message = FALSE, warning=FALSE-------------
data(resistance_surfaces)
continuous <- resistance_surfaces[[2]]
data(samples)
sample.locales <- SpatialPoints(samples[,c(2,3)])

## ----eval=FALSE----------------------------------------------------------
#  data(resistance_surfaces)
#  continuous <- resistance_surfaces[[2]]
#  writeRaster(continuous,
#              filename = paste0(write.dir,"cont.asc"),
#              overwrite = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  data(samples)
#  write.table(samples,file=paste0(write.dir,"samples.txt"),sep="\t",col.names=F,row.names=F)
#  
#  # Create a spatial points object for plotting
#  sample.locales <- SpatialPoints(samples[,c(2,3)])

## ----single.surface.plot, fig.height = 4, fig.width = 4------------------
plot(continuous)
plot(sample.locales, pch = 16, col = "blue", add = TRUE) # Add points

## ----eval=FALSE----------------------------------------------------------
#  # Set the random number seed to reproduce the results presented
#  GA.inputs <- GA.prep(ASCII.dir = write.dir,
#                       max.cat = 500,
#                       max.cont = 500,
#                       select.trans = "M",
#                       method = "LL"
#                       seed = 555)
#  
#  CS.inputs <- CS.prep(n.Pops = length(sample.locales),
#                     CS_Point.File = paste0(write.dir,"samples.txt"),
#                     CS.program = CS.program)

## ----monomolec.plot, eval = FALSE----------------------------------------
#  r.tran <- Resistance.tran(transformation = "Monomolecular",
#                            shape = 2,
#                            max = 275,
#                            r = continuous)
#  
#  plot.t <- Plot.trans(PARM = c(2, 275),
#                       Resistance = continuous,
#                       transformation = "Monomolecular")

## ----eval=FALSE----------------------------------------------------------
#  # Create the true resistance/response surface
#  CS.response <- Run_CS(CS.inputs = CS.inputs,
#                        GA.inputs = GA.inputs,
#                        r = r.tran)

## ----eval=FALSE----------------------------------------------------------
#  CS.inputs <- CS.prep(n.Pops = length(sample.locales),
#                     response = CS.response,
#                     CS_Point.File = paste0(write.dir,"samples.txt"),
#                     CS.program = CS.program)

## ----eval=FALSE----------------------------------------------------------
#  SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
#                         GA.inputs=GA.inputs)

## ----eval=FALSE----------------------------------------------------------
#  GA.inputs <- GA.prep(ASCII.dir = write.dir)
#  
#  CS.inputs <- CS.prep(n.Pops = length(sample.locales),
#                       response = CS.response,
#                       CS_Point.File = paste0(write.dir,"samples.txt"),
#                       CS.program = CS.program)
#  
#  SS_RESULTS <- SS_optim(CS.inputs = CS.inputs,
#                         GA.inputs = GA.inputs)

## ----gdistance, eval=FALSE-----------------------------------------------
#  # Import data
#  data(resistance_surfaces)
#  continuous <- resistance_surfaces[[2]]
#  
#  data(samples)
#  sample.locales <- SpatialPoints(samples[,c(2,3)])
#  
#  # Set the random number seed to reproduce the results presented
#  # Run in parallel on 4 cores
#  GA.inputs <- GA.prep(ASCII.dir=continuous,
#                       Results.dir=write.dir,
#                       select.trans = "M",
#                       max.cat=500,
#                       max.cont=500,
#                       seed = 555,
#                       parallel = 4)
#  
#  
#  gdist.inputs <- gdist.prep(length(sample.locales),
#                             samples = sample.locales,
#                             response = CS.response,
#                             method = 'commuteDistance') ## Optimize using commute distance
#  
#  # Run optimization
#  SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
#                               GA.inputs = GA.inputs)

## ----eval=FALSE----------------------------------------------------------
#  Grid.Results <- Grid.Search(shape = seq(1, 3, by = 0.025),
#                              max = seq(125, 425, by = 50),
#                              transformation = "Monomolecular",
#                              Resistance = continuous,
#                              gdist.inputs = gdist.inputs,
#                              GA.inputs = GA.inputs)

## ----recreate.grid, eval=FALSE-------------------------------------------
#  filled.contour(Grid.Results$Plot.data,
#                 col = rainbow(14),
#                 xlab = "Shape parameter",
#                 ylab = "Maximum value parameter")

## ----warning=FALSE, results='hide',message=FALSE, eval=FALSE-------------
#  if("ResistanceGA_Examples"%in%dir("C:/")==FALSE)
#    dir.create(file.path("C:/", "ResistanceGA_Examples"))
#  
#  # Create a subdirectory for the second example
#  dir.create(file.path("C:/ResistanceGA_Examples/","MultipleSurfaces"))
#  
#  write.dir <- "C:/ResistanceGA_Examples/MultipleSurfaces/"      # Directory to write .asc files and results

## ----multi_surface.sim, warning=FALSE, message=FALSE,results='hide'------
data(resistance_surfaces)
data(samples)
sample.locales <- SpatialPoints(samples[ ,c(2,3)])

## ----feature.sim, warning=FALSE,message=FALSE, fig.height = 4, fig.width = 4----
plot(resistance_surfaces[[1]],main = resistance_surfaces[[1]]@data@names)
plot(sample.locales, pch=16, col="blue", add=TRUE)
plot(resistance_surfaces[[2]],main = resistance_surfaces[[2]]@data@names)
plot(sample.locales, pch=16, col="blue", add=TRUE)
plot(resistance_surfaces[[3]],main = resistance_surfaces[[3]]@data@names)
plot(sample.locales, pch=16, col="blue", add=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  ## Note that the `resistance_surfaces` is already a RasterStack object.
#  ## The code below for demonstration of how to make a stack.
#  r.stack <- stack(resistance_surfaces$categorical,
#                   resistance_surfaces$continuous,
#                   resistance_surfaces$feature)
#  
#  GA.inputs <- GA.prep(ASCII.dir = r.stack,
#                       Results.dir = write.dir,
#                       method = "LL",
#                       max.cat = 500,
#                       max.cont = 500,
#                       seed = 555,
#                       parallel = 4)
#  
#  gdist.inputs <- gdist.prep(length(sample.locales),
#                             samples = sample.locales,
#                             method = 'commuteDistance') ## Optimize using commute distance
#  

## ----IR_Mono, eval=FALSE-------------------------------------------------
#  plot.t <- Plot.trans(PARM = c(3.5, 150),
#                       Resistance = continuous,
#                       transformation = "Inverse-Reverse Monomolecular")

## ----combine.surfaces, eval=FALSE----------------------------------------
#  PARM <- c(1, 250, 75, 1, 3.5, 150, 1, 350)
#  
#  # PARM<- c(1,   # First feature of categorical
#  #          250, # Second feature of categorical
#  #          75,  # Third feature of categorical
#  #          1,   # Transformation equation for continuous surface = Inverse-Reverse Monomolecular
#  #          3.5,   # Shape parameter
#  #          150, # Scale parameter
#  #          1,   # First feature of feature surface
#  #          350) # Second feature of feature surface
#  
#  # Combine resistance surfaces
#  Resist <- Combine_Surfaces(PARM = PARM,
#                             gdist.inputs = gdist.inputs,
#                             GA.inputs = GA.inputs,
#                             out = NULL,
#                             rescale = TRUE)
#  
#  # View combined surface
#  plot(Resist,  main = "scaled composite resistance")

## ----combine.cs, eval=FALSE----------------------------------------------
#  # Create the true resistance/response surface
#  gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
#                                  r = Resist)
#  
#  gdist.inputs <- gdist.prep(n.Pops = length(sample.locales),
#                             samples = sample.locales,
#                             response = as.vector(gdist.response),
#                             method = 'commuteDistance')

## ----eval=FALSE----------------------------------------------------------
#  Multi.Surface_optim <- MS_optim(gdist.inputs = gdist.inputs,
#                                  GA.inputs = GA.inputs)

## ----eval=FALSE----------------------------------------------------------
#  Summary.table <- data.frame(PARM,round(t(Multi.Surface_optim$GA.summary@solution),2))
#  colnames(Summary.table)<-c("Truth", "Optimized")
#  row.names(Summary.table)<-c("Category1", "Category2", "Category3", "Transformation", "Shape", "Max", "Feature1", "Feature2")

## ----combined.plots,fig.width=12,fig.height=8, eval=FALSE----------------
#  # Make combined, optimized resistance surface.
#  optim.resist <- Combine_Surfaces(PARM = Multi.Surface_optim$GA.summary@solution,
#                                   gdist.inputs = gdist.inputs,
#                                   GA.inputs =  GA.inputs,
#                                   rescale = TRUE)
#  ms.stack <- stack(Resist, optim.resist)
#  names(ms.stack) <- c("Truth", "Optimized")
#  plot(ms.stack)
#  
#  # Correlation between the two surfaces
#  pairs(ms.stack)

## ----CS.maps, eval=FALSE-------------------------------------------------
#  ## Note: You must run `CS.prep` to generate the CS.inputs object for doing this.
#  CS.inputs <- CS.prep(n.Pops = length(sample.locales),
#                       response = gdist.response,
#                       CS_Point.File = paste0(write.dir,"samples.txt"),
#                       CS.program = CS.program)
#  
#  Resist.true <- Run_CS(CS.inputs = CS.inputs,
#                        GA.inputs = GA.inputs,
#                        r = Resist,
#                        CurrentMap = TRUE,
#                        output = "raster")
#  
#  Resist.opt <- Run_CS(CS.inputs = CS.inputs,
#                       GA.inputs = GA.inputs,
#                       r = optim.resist,
#                       CurrentMap = TRUE,
#                       output = "raster")
#  
#  ## We can confirm that, like the resistance surfaces above,
#  ## the CIRCUITSCAPE current maps are also correlated
#  cs.stack <- stack(Resist.true, Resist.opt)
#  names(cs.stack) <- c("Truth", "Optimized")
#  pairs(cs.stack)

## ----eval = FALSE--------------------------------------------------------
#  data(resistance_surfaces)
#  cat <- resistance_surfaces[[1]]
#  cat[cat < 2] <- 0
#  
#  ## Make categorical surface binary
#  cat[cat == 2] <- 1
#  
#  ## Smooth and visualize
#  ## The `SCALE` parameter re-scales the surface to 0-10
#  cat.smooth <- k.smooth(raster = cat,
#                         sigma = 1,
#                         SCALE = TRUE)
#  par(mfrow = c(1,2))
#  plot(cat, main = "Original")
#  plot(cat.smooth, main = "Smoothed, sigma = 1")
#  par(mfrow = c(1,1))

## ----k_smooth, eval=FALSE------------------------------------------------
#  data(samples)
#  sample.locales <- SpatialPoints(samples[,c(2,3)])
#  
#  ## Set the random number seed to reproduce the results presented
#  ## Run in parallel on 4 cores
#  ## NOTE: `scale = TRUE` to indicate optimization of scaling/smoothing parameter
#  GA.inputs <- GA.prep(ASCII.dir = cat,
#                       Results.dir = write.dir,
#                       select.trans = "M",
#                       scale = TRUE,
#                       max.cat = 500,
#                       max.cont = 500,
#                       seed = 321,
#                       run = 35,
#                       parallel = 4)
#  
#  ## Optimize using commute distance
#  gdist.inputs <- gdist.prep(n.Pops = length(sample.locales),
#                             samples = sample.locales,
#                             method = 'commuteDistance')
#  
#  # Transform resistance surface
#  r.tran_smooth <- Resistance.tran(transformation = "Monomolecular",
#                                   shape = 2,
#                                   max = 275,
#                                   r = cat.smooth)
#  
#  # Create the true resistance/response surface
#  gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
#                                  r = r.tran_smooth)
#  
#  # Rerun `gdist.prep` to include response
#  gdist.inputs <- gdist.prep(n.Pops = length(sample.locales),
#                             response = as.vector(gdist.response),
#                             samples = sample.locales,
#                             method = 'commuteDistance')
#  
#  # Run optimization: NOTE use of `SS_optim.scale` to optimize the smoothing parameter
#  SS_RESULTS.gdist_scale <- SS_optim.scale(gdist.inputs = gdist.inputs,
#                                           GA.inputs = GA.inputs)

## ----analysis1, eval=FALSE-----------------------------------------------
#  data(samples)
#  data("resistance_surfaces")
#  
#  # Create a spatial points object
#  sample.locales <- SpatialPoints(samples[, c(2, 3)])
#  
#  # Set output directory
#  write.dir <-
#    "C:/ResistanceGA_Examples/run1/"      # Directory to write .asc files and results
#  
#  
#  # Run `gdist.prep` & GA.prep
#  gdist.inputs <-  gdist.prep(n.Pops = length(sample.locales),
#                              samples = sample.locales,
#                              method = 'commuteDistance')
#  
#  ## This will be used again later
#  GA.inputs_NoFeature <- GA.prep(method = "LL",
#                                 ASCII.dir = resistance_surfaces[[-3]],
#                                 Results.dir = "C:/ResistanceGA_Examples/run2/",
#                                 max.cat = 500,
#                                 max.cont = 500,
#                                 seed = 555,
#                                 parallel = 4)
#  
#  ## The 'true' resistance surface will be the composite surface
#  ## Combine resistance surfaces, omitting the feature surface
#  ## Use an Inverse Ricker transformation of the continuous surface
#  ## Inverse Ricker  = 8
#  PARM <- c(1, 250, 75, 8, 4, 150)
#  
#  Resist <- Combine_Surfaces(PARM = PARM,
#                             gdist.inputs = gdist.inputs,
#                             GA.inputs = GA.inputs_NoFeature,
#                             out = NULL,
#                             rescale = TRUE,
#                             p.contribution = T)
#  
#  ## Assess contribution of each surface
#  Resist$percent.contribution
#  

## ----analysis2, eval=FALSE-----------------------------------------------
#  # Create the true response
#  gd.true <- Run_gdistance(gdist.inputs = gdist.inputs,
#                               r = Resist$combined.surface)
#  gd.true <- as.vector(gd.true)
#  
#  ## Add some noise to response
#  set.seed(321)
#  gd.response <- gd.true + rnorm(length(gd.true), 0, 5)
#  
#  plot(gd.response ~ gd.true)
#  ecodist::mantel(gd.response ~ gd.true) # Mantel r = 0.58
#  
#  ## Re-run `gdist.prep` function
#  gdist.inputs <- gdist.prep(n.Pops = length(sample.locales),
#                             response = gd.response,
#                             samples = sample.locales)
#  
#  
#  ## Re-run GA.prep to include all surfaces
#  GA.inputs_All <- GA.prep(method = "LL",
#                                 ASCII.dir = resistance_surfaces,
#                                 Results.dir = write.dir,
#                                 max.cat = 500,
#                                 max.cont = 500,
#                                 seed = 555,
#                                 parallel = 4)
#  
#  ## First run all single surfaces, Multi-surface is response variable
#  SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
#                                GA.inputs = GA.inputs_All)
#  
#  ## Run `MS_optim` with all surfaces
#  Multi.Surface_optim.gd <- MS_optim(gdist.inputs = gdist.inputs,
#                                     GA.inputs = GA.inputs_All)
#  
#  ## Run `MS_optim` with without Feature surface
#  Multi.Surface_optim.gd2 <- MS_optim(gdist.inputs = gdist.inputs,
#                                      GA.inputs = GA.inputs_NoFeature)

## ----analysis3, eval=FALSE-----------------------------------------------
#  Multi.Surface_optim.gd$percent.contribution
#  
#  Multi.Surface_optim.gd2$percent.contribution

## ----boot, eval = F------------------------------------------------------
#  ## Extract relevant components from optimization outputs
#  ## Make a list of cost/resistance distance matrices
#  mat.list <- c(SS_RESULTS.gdist$cd,
#                Multi.Surface_optim.gd$cd,
#                Multi.Surface_optim.gd2$cd)
#  
#  k <- rbind(SS_RESULTS.gdist$k,
#             Multi.Surface_optim.gd$k,
#             Multi.Surface_optim.gd2$k)
#  
#  ## Make square genetic distance matrix
#  g.mat <- matrix(rep(0, 25^2),nrow = 25)
#  g.mat[lower.tri(g.mat)] <- gd.response
#  
#  ## Run bootstrap
#  (AIC.boot <- Resist.boot(mod.names = names(mat.list),
#                           dist.mat = mat.list,
#                           n.parameters = k[,2],
#                           sample.prop = 0.75,
#                           iters = 1000,
#                           obs = 25,
#                           genetic.mat = gd.response
#  )
#  )
#  

## ---- eval = FALSE-------------------------------------------------------
#  Summary.table <- data.frame(PARM,round(t(Multi.Surface_optim.gd$GA.summary@solution),2))
#  colnames(Summary.table)<-c("Truth", "Optimized")
#  row.names(Summary.table)<-c("Category1", "Category2", "Category3", "Transformation", "Shape", "Max", "Feature1", "Feature2")

## ---- eval = FALSE-------------------------------------------------------
#  opt.r <- raster("C:/ResistanceGA_Examples/run1/categorical.continuous.feature.asc")
#  r.stack <- stack(Resist, opt.r)
#  names(r.stack) <- c("Truth", "Optimized")
#  
#  plot(r.stack)
#  
#  pairs(r.stack)

