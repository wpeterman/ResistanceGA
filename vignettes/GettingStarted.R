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

