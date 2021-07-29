## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("wpeterman/ResistanceGA", build_vignettes = TRUE)

## ----fig.width=6, fig.asp=0.65------------------------------------------------
library(ResistanceGA)

plot(raster_orig)

plot(raster_true)

plot(raster_true$multi_true, main = 'Multivariate Resistance')
plot(sample_pops$sample_multi, add = T, pch = 19, col = 'blue')

par(mfrow = c(2,1))
plot(lower(Dc_list$Dc_multi) ~
       c(dist(sample_pops$sample_multi@coords)),
     xlab = 'Euclidean distance', ylab = 'Chord distance')

plot(lower(Dc_list$Dc_multi) ~
       lower(resist_list$resist_multi),
     xlab = 'Resistance distance', ylab = 'Chord distance')
par(mfrow = c(1,1))


## ----eval = F-----------------------------------------------------------------
#  GA.prep(ASCII.dir,            # REQUIRED
#          Results.dir = NULL,   # REQUIRED
#          min.cat = NULL,
#          max.cat = 2500,
#          max.cont = 2500,
#          min.scale = NULL,
#          max.scale = NULL,
#          cont.shape = NULL,
#          select.trans = 'M',
#          method = "LL",
#          scale = FALSE,
#          scale.surfaces = NULL,
#          k.value = 2,
#          pop.mult = 15,
#          percent.elite = 0.05,
#          type = "real-valued",
#          pcrossover = 0.85,
#          pmutation = 0.125,
#          maxiter = 1000,
#          run = NULL,
#          keepBest = TRUE,
#          population = gaControl(type)$population,
#          selection = gaControl(type)$selection,
#          crossover = "gareal_blxCrossover",
#          mutation = gaControl(type)$mutation,
#          pop.size = NULL,
#          parallel = FALSE,
#          gaisl = FALSE,
#          island.pop = 20,
#          numIslands = NULL,
#          migrationRate = NULL,
#          migrationInterval = NULL,
#          optim = FALSE,
#          optim.method = "L-BFGS-B",
#          poptim = 0.0,
#          pressel = 1.00,
#          control = list(fnscale = -1, maxit = 100),
#          hessian = FALSE,
#          opt.digits = NULL,
#          seed = NULL,
#          monitor = TRUE,
#          quiet = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  gdist.prep(n.Pops,           # REQUIRED
#             response = NULL,  # REQUIRED
#             samples,          # REQUIRED
#             covariates = NULL,
#             formula = NULL,
#             transitionFunction = function(x)  1 / mean(x),
#             directions = 8,
#             longlat = FALSE,
#             method = 'commuteDistance',
#             min.max_dist = NULL,
#             keep = NULL)

## ----eval = FALSE-------------------------------------------------------------
#  jl.prep(n.Pops,             # REQUIRED
#          response = NULL,    # REQUIRED
#          CS_Point.File,      # REQUIRED
#          covariates = NULL,
#          formula = NULL,
#          JULIA_HOME = NULL,  # REQUIRED
#          Neighbor.Connect = 8,
#          pairs_to_include = NULL,
#          pop2ind = NULL,
#          nb = NULL,
#          parallel = FALSE,
#          cores = NULL,
#          cholmod = TRUE,
#          precision = FALSE,
#          run_test = TRUE,
#          write.files = NULL,
#          write.criteria = NULL,
#          silent = TRUE,
#          Julia_link = 'JuliaCall',
#          scratch = NULL,
#          rm.files = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  Plot.trans(PARM,
#             Resistance,
#             transformation,
#             scale,
#             print.dir,
#             marginal.plot,
#             marg.type,
#             Name)

## ----fig.width=6, fig.asp=0.65------------------------------------------------
(Plot.trans(PARM = c(14.5, 100), 
            Resistance = raster_orig$cont_orig, 
            transformation = 1))

## This creates the same plot
# (Plot.trans(PARM = c(2.5, 100), 
#             Resistance = raster_orig$cont_orig, 
#             transformation = "Reverse Monomolecular"))

(Plot.trans(PARM = c(2.5, 100), 
            Resistance = c(1,10), 
            transformation = 1))

## ----eval=F-------------------------------------------------------------------
#  gdist.inputs <- gdist.prep(n.Pops = length(sample_pops$sample_cont),
#                             samples = sample_pops$sample_cont)
#  
#  gdist.out <- Run_gdistance(gdist.inputs = gdist.inputs,
#                             r = raster_true$cont_true)

## ---- echo=F, eval=T----------------------------------------------------------
try(jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cont),
                         CS_Point.File = sample_pops$sample_cont,
                     JULIA_HOME = "C:/Users/peterman.73/AppData/Local/Programs/Julia-1.6.1/bin/"), silent = T)

## ----eval=F-------------------------------------------------------------------
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cont),
#                       CS_Point.File = sample_pops$sample_cont,
#                       JULIA_HOME = "C:/Users/peterman.73/AppData/Local/Programs/Julia-1.6.1/bin/")
#  jl.out <- Run_CS.jl(jl.inputs = jl.inputs,
#                      r = raster_true$cont_true)

## ----fig.width=6, fig.asp=0.65------------------------------------------------
jl.out <- Run_CS.jl(jl.inputs = jl.inputs,
                    r = raster_true$cont_true,
                    CurrentMap = TRUE,
                    output = 'raster',
                    EXPORT.dir = "C:/Rga/")
plot(jl.out)

## ----eval=F-------------------------------------------------------------------
#  # Get lower half of distance matrices
#  gd <- lower(Dc_list$Dc_cont) # Genetic distance
#  rd <- lower(resist_list$resist_cont) # Resistance distance
#  
#  ## Run `To.From.ID`
#  id <- To.From.ID(length(sample_pops$sample_cont))
#  
#  df <-data.frame(gd = gd,
#                  rd = rd,
#                  pop = id$pop1)
#  
#  mlpe.mod <- mlpe_rga(gd ~ rd + (1 | pop),
#                       data = df)

## ---- eval = F, echo=T--------------------------------------------------------
#  dir.create("C:/Rga/ExampleAnalysis/", recursive = T)
#  results_dir <- "C:/Rga/ExampleAnalysis/"
#  
#  ## For simplicity, pull out the genetic and spatial data from the package data objects
#  gen_dist <- lower(Dc_list$Dc_cont)
#  r_stack <- raster_orig
#  sample_locs <- sample_pops$sample_cont
#  
#  ## First, run GA.prep . In the interest of time, we'll reduce the total number of iterations (`maxiter`) to 10. Run in `parallel` on at least 2 cores, more if available.
#  ## Check # of cores available to use
#  ## DON'T USE ALL CORES!
#  parallel::detectCores()
#  
#  ##  ** [You fill in the function!] **
#  GA.inputs <- GA.prep(ASCII.dir = r_stack,
#                       Results.dir = results_dir,
#                       pop.size = 25,
#                       maxiter = 10,
#                       parallel = 5)
#  
#  
#  ## Next, run prep for your optimization method. We'll use `gdistance` here, as it should work seamlessly for all
#  ## Feel free to use `jl.prep` and Julia instead, if you have everything installed and working.
#  ##  ** [You fill in the function!] **
#  gdist.inputs <- gdist.prep(n.Pops = length(sample_locs),
#                             response = gen_dist,
#                             samples = sample_locs)
#  
#  ## Alternatively
#  # jl.inputs <- jl.prep()
#  
#  
#  ## Finally, optimize your surface(s)
#  ## Single Surface
#  ss_results <- SS_optim(gdist.inputs = gdist.inputs,
#                         GA.inputs = GA.inputs)
#  
#  ## Multiple Surfaces
#  ms_results <- MS_optim(gdist.inputs = gdist.inputs,
#                         GA.inputs = GA.inputs)
#  
#  ## If optimizing all surfaces together, it may make sense to use the `all_comb` to have a more comprehensive analysis. We'll need to adjust our `GA.prep` function `Results.dir` to accommodate this. Set the `iters` parameter to 100.
#  
#  GA.inputs2 <- GA.prep(ASCII.dir = r_stack,
#                       Results.dir = 'all_comb',
#                       pop.size = 25,
#                       maxiter = 10,
#                       parallel = 5)
#  
#  ac_results <- all_comb(gdist.inputs = gdist.inputs,
#                         GA.inputs = GA.inputs2,
#                         results.dir = results_dir,
#                         iters = 100)
#  
#  

## ----eval = F, echo=T---------------------------------------------------------
#  ## Path to your Julia installation
#  JULIA_HOME <- "C:/Users/peterman.73/AppData/Local/Programs/Julia-1.6.1/bin/"
#  
#  # Single Surface ----------------------------------------------------------
#  
#  # # >> Categorical ----------------------------------------------------------
#  dir.create("C:/Rga/single/cat5", recursive = T)
#  
#  GA.inputs <- GA.prep(ASCII.dir = raster_orig$cat_orig,
#                       Results.dir = "C:/Rga/single/cat5/",
#                       max.cat = 500,
#                       parallel = 5)
#  
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cat),
#                       response = lower(Dc_list$Dc_cat),
#                       CS_Point.File = sample_pops$sample_cat,
#                       JULIA_HOME = JULIA_HOME,
#                       cholmod = T)
#  
#  ss.cat <- SS_optim(jl.inputs = jl.inputs,
#                     GA.inputs = GA.inputs)
#  
#  # >> Continuous -----------------------------------------------------------
#  dir.create("C:/Rga/single/cont", recursive = T)
#  
#  GA.inputs <- GA.prep(ASCII.dir = raster_orig$cont_orig,
#                       Results.dir = "C:/Rga/single/cont/",
#                       max.cont = 500,
#                       parallel = 5)
#  
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cont),
#                       response = lower(Dc_list$Dc_cont),
#                       CS_Point.File = sample_pops$sample_cont,
#                       JULIA_HOME = JULIA_HOME,
#                       cholmod = T)
#  
#  ss.cont <- SS_optim(jl.inputs = jl.inputs,
#                      GA.inputs = GA.inputs)
#  
#  
#  
#  
#  # >> Continuous: Keep -----------------------------------------------------
#  # Sometimes we have samples that are very close together spatially and we don't necessarily expect there to be a prominent or meaningful effect of the landscape between between them. In these instances, we can identify which pairwise distances we want to retain, and omit the others from inclusion in our MLPE model.
#  
#  dir.create("C:/Rga/single/cont_keep", recursive = T)
#  
#  GA.inputs <- GA.prep(ASCII.dir = raster_orig$cont_orig,
#                       Results.dir = "C:/Rga/single/cont_keep/",
#                       max.cont = 500,
#                       parallel = 5)
#  
#  ## Keep only pairs that are > 4 distance units apart
#  keep <- c(dist(sample_pops$sample_cont@coords))
#  keep[keep < 4] <- 0
#  keep[keep > 0] <- 1
#  
#  length(keep)
#  sum(keep)
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cont),
#                       response = lower(Dc_list$Dc_cont),
#                       CS_Point.File = sample_pops$sample_cont,
#                       JULIA_HOME = JULIA_HOME,
#                       pairs_to_include = keep,
#                       cholmod = T)
#  
#  ss.cont_keep <- SS_optim(jl.inputs = jl.inputs,
#                           GA.inputs = GA.inputs)
#  
#  
#  
#  # >> Categorical: Thematic Resolution -------------------------------------
#  ## Another application of ResistanceGA is identification of thematic resolution of categorical land cover surfaces.In these analyses, we try to determine the best way to characterize our land cover surface. For this example, we reclassify category 5 from the original categorical raster to be part of category 1.
#  
#  dir.create("C:/Rga/single/cat4", recursive = T)
#  
#  cat4 <- subs(raster_orig$cat_orig, data.frame(c(1,2,3,4,5),
#                                                c(1,2,3,4,1)))
#  names(cat4) <- 'cat4'
#  cat_stack <- stack(cat4,
#                     raster_orig$cat_orig)
#  
#  GA.inputs <- GA.prep(ASCII.dir = cat_stack,
#                       Results.dir = "C:/Rga/single/cat4/",
#                       max.cat = 500,
#                       parallel = 8)
#  
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cat),
#                       response = lower(Dc_list$Dc_cat),
#                       CS_Point.File = sample_pops$sample_cat,
#                       JULIA_HOME = JULIA_HOME,
#                       cholmod = T)
#  
#  ss.cat4 <- SS_optim(jl.inputs = jl.inputs,
#                      GA.inputs = GA.inputs)
#  

## ----eval = F, echo=T---------------------------------------------------------
#  # Multisurface ---------------------------------------------------------
#  
#  
#  # >> Multi ----------------------------------------------------------------
#  dir.create("C:/Rga/multi/multi", recursive = T)
#  
#  GA.inputs <- GA.prep(ASCII.dir = raster_orig,
#                       Results.dir = "C:/Rga/multi/multi/",
#                       max.cat = 500,
#                       max.cont = 500,
#                       parallel = 5)
#  
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_multi),
#                       response = lower(Dc_list$Dc_multi),
#                       CS_Point.File = sample_pops$sample_multi,
#                       JULIA_HOME = JULIA_HOME,
#                       cholmod = T)
#  
#  ms <- MS_optim(jl.inputs = jl.inputs,
#                 GA.inputs = GA.inputs)
#  
#  
#  # All combination ---------------------------------------------------------
#  
#  ## Sometimes we want to go beyond just comparing single surfaces or running specific multi-surface combinations.In these instances, you may want to use the `all_comb` function to run all combinations of your surfaces.
#  
#  # >> cat ------------------------------------------------------------------
#  
#  
#  dir.create("C:/Rga/multi/all_comb_cat", recursive = T)
#  
#  GA.inputs <- GA.prep(ASCII.dir = raster_orig,
#                       Results.dir = 'all_comb',
#                       max.cat = 500,
#                       max.cont = 500,
#                       parallel = 5)
#  
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cat),
#                       response = lower(Dc_list$Dc_cat),
#                       CS_Point.File = sample_pops$sample_cat,
#                       JULIA_HOME = JULIA_HOME,
#                       cholmod = T)
#  
#  ac_cat <- all_comb(jl.inputs = jl.inputs,
#                     GA.inputs = GA.inputs,
#                     results.dir = "C:/Rga/multi/all_comb_cat/",
#                     replicate = 3)
#  
#  
#  
#  # >> cont -----------------------------------------------------------------
#  
#  
#  dir.create("C:/Rga/multi/all_comb_cont", recursive = T)
#  
#  GA.inputs <- GA.prep(ASCII.dir = raster_orig,
#                       Results.dir = 'all_comb',
#                       max.cat = 500,
#                       max.cont = 500,
#                       parallel = 5)
#  
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cont),
#                       response = lower(Dc_list$Dc_cont),
#                       CS_Point.File = sample_pops$sample_cont,
#                       JULIA_HOME = JULIA_HOME,
#                       cholmod = T)
#  
#  ac_cont <- all_comb(jl.inputs = jl.inputs,
#                      GA.inputs = GA.inputs,
#                      results.dir = "C:/Rga/multi/all_comb_cont/",
#                      replicate = 3)
#  
#  
#  # >> multi ----------------------------------------------------------------
#  
#  
#  dir.create("C:/Rga/multi/all_comb_multi", recursive = T)
#  
#  GA.inputs <- GA.prep(ASCII.dir = raster_orig,
#                       Results.dir = 'all_comb',
#                       max.cat = 500,
#                       max.cont = 500,
#                       parallel = 5)
#  
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_multi),
#                       response = lower(Dc_list$Dc_multi),
#                       CS_Point.File = sample_pops$sample_multi,
#                       JULIA_HOME = JULIA_HOME,
#                       cholmod = T)
#  
#  ac_multi <- all_comb(jl.inputs = jl.inputs,
#                       GA.inputs = GA.inputs,
#                       results.dir = "C:/Rga/multi/all_comb_multi/",
#                       replicate = 3)
#  
#  # >> Continuous linear -----------------------------------------------------
#  dir.create("C:/Rga/single/cont_linear", recursive = T)
#  
#  Plot.trans(PARM = c(15,100),
#             Resistance = c(1,10),
#             transformation = 1)
#  
#  GA.inputs <- GA.prep(ASCII.dir = raster_orig$cont_orig,
#                       Results.dir = "C:/Rga/single/cont_linear/",
#                       max.cont = 500,
#                       shape.min = 14.9,
#                       shape.max = 15.1,
#                       parallel = 5)
#  
#  
#  jl.inputs <- jl.prep(n.Pops = length(sample_pops$sample_cont),
#                       response = lower(Dc_list$Dc_cont),
#                       CS_Point.File = sample_pops$sample_cont,
#                       JULIA_HOME = JULIA_HOME,
#                       ss.cont_lin <- SS_optim(jl.inputs = jl.inputs,
#                                               GA.inputs = GA.inputs)

