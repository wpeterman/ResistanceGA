## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install.package, eval=FALSE-----------------------------------------
#  # Install 'devtools' package, if needed
#  if(!("devtools" %in% list.files(.libPaths()))) {
#      install.packages("devtools", repo = "http://cran.rstudio.com", dep = TRUE)
#  }
#  
#  library(devtools) # Loads devtools
#  
#  install_github("wpeterman/ResistanceGA") # Download package

## ----results='hide',message=FALSE, warning=FALSE-------------------------
library(ResistanceGA)
rm(list = ls())

## ----results='hide', eval=FALSE------------------------------------------
#  Ricker.plot <- Plot.trans(PARM=c(1.5, 200),Resistance=c(1,10),transformation="Ricker")

## ----results='hide', eval=FALSE------------------------------------------
#  # Change title of plot
#  Ricker.plot$labels$title<-"Ricker Transformation"
#  Ricker.plot

## ----max_resist, eval=FALSE----------------------------------------------
#  # Find original data value that now has maximum resistance
#  Ricker.plot$data$original[which(Ricker.plot$data$transformed==max(Ricker.plot$data$transformed))]
#  
#  ## [1] 2.356784

## ----warning=FALSE, results='hide',message=FALSE, eval=FALSE-------------
#  
#  if("ResistanceGA_Examples"%in%dir("C:/")==FALSE)
#    dir.create(file.path("C:/", "ResistanceGA_Examples"))
#  
#  # Create a subdirectory for the first example
#  dir.create(file.path("C:/ResistanceGA_Examples/","SingleSurface"))
#  
#  # Directory to write .asc files and results
#  write.dir <- "C:/ResistanceGA_Examples/SingleSurface/"
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
#  writeRaster(continuous,filename=paste0(write.dir,"cont.asc"),overwrite=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  data(samples)
#  write.table(samples,file=paste0(write.dir,"samples.txt"),sep="\t",col.names=F,row.names=F)
#  
#  # Create a spatial points object for plotting
#  sample.locales <- SpatialPoints(samples[,c(2,3)])

## ----single.surface.plot, eval = FALSE-----------------------------------
#  plot(continuous)
#  plot(sample.locales, pch=16, col="blue", add=TRUE) # Add points

## ----eval=FALSE----------------------------------------------------------
#  # Set the random number seed to reproduce the results presented
#  GA.inputs <- GA.prep(ASCII.dir = write.dir,
#                       max.cat = 500,
#                       max.cont = 500,
#                       select.trans = "A",
#                       seed = 555)
#  
#  CS.inputs <- CS.prep(n.Pops = length(sample.locales),
#                       CS_Point.File = paste0(write.dir,"samples.txt"),
#                       CS.program = CS.program)

## ----monomolec.plot, eval = FALSE----------------------------------------
#  r.tran <- Resistance.tran(transformation="Monomolecular", shape=2, max=275, r=continuous)
#  
#  plot.t <- Plot.trans(PARM=c(2,275), Resistance=continuous, transformation="Monomolecular")

## ----eval=FALSE----------------------------------------------------------
#  # Create the true resistance/response surface
#  CS.response <- Run_CS(CS.inputs=CS.inputs,GA.inputs=GA.inputs, r=r.tran)

## ----eval=FALSE----------------------------------------------------------
#  CS.inputs <- CS.prep(n.Pops = length(sample.locales),
#                       response = CS.response,
#                       CS_Point.File = paste0(write.dir,"samples.txt"),
#                       CS.program = CS.program)

## ----eval=FALSE----------------------------------------------------------
#  SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
#                         GA.inputs=GA.inputs)

## ----eval=FALSE----------------------------------------------------------
#  Grid.Results <- Grid.Search(shape=seq(1,4,by=0.1),
#                              max=seq(50,500,by=75),
#                              transformation="Monomolecular",
#                              Resistance=continuous,
#                              CS.inputs,
#                              GA.inputs=GA.inputs)

## ----recreate.grid, eval=FALSE-------------------------------------------
#  filled.contour(Grid.Results$Plot.data,
#                 col=rainbow(30),
#                 xlab="Shape parameter",
#                 ylab="Maximum value parameter")

## ----eval=FALSE----------------------------------------------------------
#  GA.inputs <- GA.prep(ASCII.dir=write.dir)
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
#                       max.cat=500,
#                       max.cont=500,
#                       seed = 555,
#                       parallel = 4)
#  
#  gdist.inputs <- gdist.prep(length(sample.locales),
#                             samples=sample.locales)
#  
#  # Transform resistance surface
#  r.tran <- Resistance.tran(transformation="Monomolecular", shape=2, max=275, r=continuous)
#  
#  # Create the true resistance/response surface
#  gdist.response <- Run_gdistance(gdist.inputs=gdist.inputs, r=r.tran)
#  
#  # Rerun `gdist.prep` to include response
#  gdist.inputs <- gdist.prep(n.Pops = length(sample.locales),
#                             response=lower(as.matrix(gdist.response)),
#                             samples=sample.locales)
#  
#  # Run optimization
#  SS_RESULTS.gdist <- SS_optim(gdist.inputs=gdist.inputs,
#                               GA.inputs=GA.inputs)
#  
#  # Grid search of response surface
#  Grid.Results.gdist <- Grid.Search(shape=seq(1,4,by=0.1),
#                                    max=seq(50,500,by=75),
#                                    transformation="Monomolecular",
#                                    Resistance=continuous,
#                                    gdist.inputs=gdist.inputs,
#                                    GA.inputs=GA.inputs)

## ----warning=FALSE, results='hide',message=FALSE, eval=FALSE-------------
#  if("ResistanceGA_Examples"%in%dir("C:/")==FALSE)
#    dir.create(file.path("C:/", "ResistanceGA_Examples"))
#  
#  # Create a subdirectory for the second example
#  dir.create(file.path("C:/ResistanceGA_Examples/","MultipleSurfaces"))
#  
#  # Directory to write .asc files and results
#  write.dir <- "C:/ResistanceGA_Examples/MultipleSurfaces/"

## ----multi_surface.sim, warning=FALSE, message=FALSE,results='hide'------
data(resistance_surfaces)
data(samples)
sample.locales <- SpatialPoints(samples[ ,c(2,3)])

## ----feature.sim, warning=FALSE,message=FALSE,eval=FALSE-----------------
#  plot(resistance_surfaces[[1]],main = resistance_surfaces[[1]]@data@names)
#  plot(sample.locales, pch=16, col="blue", add=TRUE)
#  plot(resistance_surfaces[[2]],main = resistance_surfaces[[2]]@data@names)
#  plot(sample.locales, pch=16, col="blue", add=TRUE)
#  plot(resistance_surfaces[[3]],main = resistance_surfaces[[3]]@data@names)
#  plot(sample.locales, pch=16, col="blue", add=TRUE)

## ----write_multi.ASCII, warning=FALSE, results='hide',message=FALSE, eval=FALSE----
#  writeRaster(categorical,filename=paste0(write.dir,"cat.asc"),overwrite=TRUE)
#  writeRaster(continuous,filename=paste0(write.dir,"cont.asc"),overwrite=TRUE)
#  writeRaster(feature,filename=paste0(write.dir,"feature.asc"),overwrite=TRUE)
#  
#  write.table(samples,file=paste0(write.dir,"samples.txt"),sep="\t",col.names=F,row.names=F)

## ----eval=FALSE----------------------------------------------------------
#  GA.inputs <- GA.prep(ASCII.dir=write.dir,
#                       method = "R2"
#                       max.cat=500,
#                       max.cont=500,
#                       seed = 555,
#                       quiet = TRUE)

## ----reverse.ricker, eval=FALSE------------------------------------------
#  plot.t <- Plot.trans(PARM=c(3.5,400),Resistance=continuous,transformation="Reverse Ricker")

## ----combine.surfaces, eval=FALSE----------------------------------------
#  PARM <- c(1, 250, 75, 6, 3.5, 150, 1, 350)
#  
#  # PARM<- c(1,   # First feature of categorical
#  #          250, # Second feature of categorical
#  #          75,  # Third feature of categorical
#  #          6,   # Transformation equation for continuous surface
#  #          3.5,   # Shape parameter
#  #          150, # Scale parameter
#  #          1,   # First feature of feature surface
#  #          350) # Second feature of feature surface
#  
#  # Combine resistance surfaces
#  Resist <- Combine_Surfaces(PARM=PARM,
#                             CS.inputs=CS.inputs,
#                             GA.inputs=GA.inputs,
#                             out=NULL,
#                             rescale = TRUE)
#  
#  # View combined surface
#  plot(Resist,  main = "scaled composite resistance")

## ----combine.cs, eval=FALSE----------------------------------------------
#  # Create the true resistance/response surface
#  CS.response <- Run_CS(CS.inputs=CS.inputs, GA.inputs=GA.inputs, r=Resist)
#  
#  CS.inputs<-CS.prep(n.Pops=length(sample.locales),
#                        response=CS.response,
#                        CS_Point.File=paste0(write.dir,"samples.txt"),
#                        CS.program=CS.program)

## ----eval=FALSE----------------------------------------------------------
#  Multi.Surface_optim <- MS_optim(CS.inputs=CS.inputs,
#                                 GA.inputs=GA.inputs)

## ----eval=FALSE----------------------------------------------------------
#  Summary.table <- data.frame(PARM,round(t(Multi.Surface_optim@solution),2))
#  colnames(Summary.table)<-c("Truth", "Optimized")
#  row.names(Summary.table)<-c("Category1",
#                              "Category2",
#                              "Category3",
#                              "Transformation",
#                              "Shape",
#                              "Max",
#                              "Feature1",
#                              "Feature2")

## ----combined.plots, eval=FALSE------------------------------------------
#  # Make combined, optimized resistance surface.
#  optim.resist <- Combine_Surfaces(PARM=Multi.Surface_optim@solution,
#                                   CS.inputs =  CS.inputs,
#                                   GA.inputs =  GA.inputs,
#                                   rescale = TRUE)
#  ms.stack <- stack(Resist, optim.resist)
#  names(ms.stack) <- c("Truth", "Optimized")
#  plot(ms.stack)
#  
#  # Correlation between the two surfaces
#  pairs(ms.stack)

## ----CS.maps, eval=FALSE-------------------------------------------------
#  Resist.true <- Run_CS(CS.inputs=CS.inputs,
#                        GA.inputs=GA.inputs,
#                        r=Resist,
#                        CurrentMap=TRUE,
#                        output="raster")
#  
#  Resist.opt <- Run_CS(CS.inputs=CS.inputs,
#                       GA.inputs=GA.inputs,
#                       r=optim.resist,
#                       CurrentMap=TRUE,
#                       output="raster")
#  
#  ## We can confirm that, like the resistance surfaces above,
#  ## the CIRCUITSCAPE current maps are also correlated
#  cs.stack <- stack(Resist.true, Resist.opt)
#  names(cs.stack) <- c("Truth", "Optimized")
#  pairs(cs.stack)

## ----multisurface_lcp, eval=FALSE----------------------------------------
#  # Run `gdist.prep`
#  gdist.inputs<-gdist.prep(n.Pops=length(sample.locales),
#                           samples=sample.locales)
#  
#  GA.inputs <- GA.prep(ASCII.dir=resistance_surfaces,
#                       Results.dir=write.dir,
#                       max.cat=500,
#                       max.cont=500,
#                       seed = 999,
#                       parallel = 4)
#  
#  # Combine resistance surfaces
#  Resist <- Combine_Surfaces(PARM=PARM,
#                             gdist.inputs=gdist.inputs,
#                             GA.inputs=GA.inputs,
#                             out=NULL,
#                             rescale=TRUE)
#  
#  # Create the true resistance/response surface
#  gd.response <- Run_gdistance(gdist.inputs=gdist.inputs,
#                               r=Resist)
#  
#  # Run `CS.prep` functions
#  gdist.inputs<-gdist.prep(n.Pops=length(sample.locales),
#                           response=lower(as.matrix(gd.response)),
#                           samples=sample.locales)
#  
#  # Run `MS_optim`
#  Multi.Surface_optim.gd <- MS_optim(gdist.inputs=gdist.inputs,
#                                     GA.inputs=GA.inputs)
#  
#  summary(Multi.Surface_optim.gd)

