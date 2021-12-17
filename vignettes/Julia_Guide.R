## ----Rga install, eval=FALSE--------------------------------------------------
#  devtools::install_github("wpeterman/ResistanceGA", build_vignettes = TRUE)

## ----Use Julia----------------------------------------------------------------
library(ResistanceGA)

## Sample locations
sp.dat <- sample_pops$sample_cont 

## Continuous landscape surface
cont.rast <- raster_orig$cont_orig

## Genetic distance measured between sample locations (chord distance)
gen.dist <- Dc_list$Dc_cont

plot(cont.rast)
plot(sp.dat, add = T, pch = 19)

## ----load julia, eval=FALSE---------------------------------------------------
#  ## Specify path the to `bin` directory
#  ## This path may vary depending upon Julia version or operating system
#  JULIA_HOME <- "C:/Users/peterman.73/AppData/Local/Programs/Julia 1.5.2/bin/"
#  JuliaCall::julia_setup(JULIA_HOME)

## ----prep function, eval = FALSE----------------------------------------------
#  GA.inputs <- GA.prep(ASCII.dir = cont.rast,
#                       Results.dir = "C:/Rga_examples/",
#                       parallel = 4)
#  jl.inputs <- jl.prep(n.Pops = length(sp.dat),
#                       response = lower(gen.dist),
#                       CS_Point.File = sp.dat,
#                       JULIA_HOME = JULIA_HOME)

## ----optimize,eval=FALSE------------------------------------------------------
#  jl.optim <- SS_optim(jl.inputs = jl.inputs,
#                       GA.inputs = GA.inputs)

