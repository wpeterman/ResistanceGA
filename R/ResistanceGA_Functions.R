#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# CONDUCT GRID SEARCH OF PARAMETER SPACE TO VIEW response SURFACE #!##!##!### 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#' Conduct grid search of response surface 
#' 
#' Visualize the AICc response surface
#' 
#' @param shape A vector of values for the shape parameter
#' @param max A vector of values for the maximum value parameter
#' @param transformation Transformation to apply. Can be either numeric or character of transformation name
#' @param Resistance An R Raster object, or path to a .asc file
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param write.dir Directory for writing intermeidate CIRCUITSCAPE files
#' @usage Grid.Search(shape, max, transformation, Resistance, CS.inputs, gdist.inputs, write.dir)
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @return This function will return a \code{filled.contour} plot. Additionally, an object with values that can be plotted with \code{filled.contour} to visualize the response surface
#' @details This function will perform a full factorial grid search of the values provided in the shape and max.scale vectors. Depending on the number of values provided for each, and the time it takes to run each iteration, this process may take a while to complete. \cr Suitable values for transformation:\cr
#' \tabular{ll}{
#'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#'    \tab 2 = "Inverse-Reverse Ricker"\cr
#'    \tab 3 = "Monomolecular"\cr
#'    \tab 4 = "Ricker"\cr
#'    \tab 5 = "Reverse Monomolecular"\cr
#'    \tab 6 = "Reverse Ricker"\cr
#'    \tab 7 = "Inverse Monomolecular"\cr
#'    \tab 8 = "Inverse Ricker"\cr
#'    \tab 9 = "Distance"\cr
#'    }



Grid.Search <- function(shape, max, transformation, Resistance, CS.inputs=NULL, gdist.inputs=NULL, write.dir) {
  if(class(Resistance)[1]!='RasterLayer') {  
    r <- raster(Resistance)
    r <- SCALE(r,0,10)
  } else {    
    r <-SCALE(Resistance,0,10)
  }
   
  GRID <- expand.grid(shape,max)
  RESULTS <- matrix(nrow=nrow(GRID),ncol=3); colnames(RESULTS)<-c("shape","max","AICc")
  EQ<-get.EQ(transformation)

if(!is.null(CS.inputs)){
for(i in 1:nrow(GRID)){
  AICc<-Resistance.Optimization_cont.nlm(PARM=log(c(t(GRID[i,]))),Resistance=r,equation=EQ, get.best=FALSE,CS.inputs,Min.Max='min',write.dir=write.dir)
  
  results<-as.matrix(cbind(GRID[i,],AICc))
  
  RESULTS[i,]<-results  
  }
} else {
  for(i in 1:nrow(GRID)){
    AICc<-Resistance.Optimization_cont.nlm(PARM=log(c(t(GRID[i,]))),Resistance=r,equation=EQ, get.best=FALSE,gdist.inputs=gdist.inputs,Min.Max='min',write.dir=write.dir)
    
    results<-as.matrix(cbind(GRID[i,],AICc))
    
    RESULTS[i,]<-results  
  }
}
RESULTS <- data.frame(RESULTS)
Results.mat <- interp(RESULTS$shape,RESULTS$max,RESULTS$AICc,duplicate='strip')
filled.contour(Results.mat,col=topo.colors(20),xlab="Shape parameter",ylab="Maximum value parameter")

AICc<-RESULTS
colnames(AICc)<-c("shape","max","AICc")
Results.mat<-list(Plot.data=Results.mat,AICc=AICc)

return(Results.mat)
}
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# Single command function to execute single surface optimization #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#' Single surface optimization
#' 
#' Optimize all surfaces contained in a directory using a genetic algorithm executed with the \code{\link[GA]{ga}} function in the Genetic Algorithms package \pkg{GA}
#' 
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param nlm Logical, if TRUE, the final step of optimization will use nlm to fine-tune parameter estimates. This may lead to overfitting in some cases. Default = FALSE.
#' @param dist_mod Logical, if TRUE, a Distance model will be calculated and added to the AICc output table (default = TRUE)
#' @param null_mod Logical, if TRUE, a Null model will be calculated and added to the AICc output table (default = TRUE)
#' @return This function optimizes resistance surfaces in isolation. Following optimization of all surfaces, several summary objects are created.\cr
#' \enumerate{
#' \item Diagnostic plots of model fit are output to the "Results/Plots" directory that is automatically generated within the folder containing the optimized ASCII files.
#' \item A .csv file with the Maximum Likelihood Population Effects mixed effects model coefficient estimates (MLPE_coeff_Table.csv)
#' \item Three summary .csv files are generated: CategoricalResults.csv, ContinuousResults.csv, & All_Results_AICc.csv. These tables contain AICc values and optimization summaries for each surface.
#' }
#' All results tables are also summarized in a named list ($ContinuousResults, $CategoricalResults, $AICc, $MLPE)
#' @usage SS_optim(CS.inputs, gdist.inputs, GA.inputs, nlm, dist_mod, null_mod)
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @export
SS_optim <- function(CS.inputs=NULL, gdist.inputs=NULL, GA.inputs, nlm=FALSE, dist_mod=TRUE, null_mod=TRUE){
  t1<-proc.time()[3]
  RESULTS.cat <- list() # List to store categorical results within
  RESULTS.cont <-list() # List to store continuous results within
  cnt1<-0
  cnt2<-0
  # Optimize each surface in turn
  for (i in 1:GA.inputs$n.layers){
    r<-GA.inputs$Resistance.stack[[i]]
    names(r)<-GA.inputs$layer.names[i]
    
    # Processing of categorical surfaces  
    if(!is.null(CS.inputs)){
      if(GA.inputs$parallel!=FALSE) {warning("\n CIRCUITSCAPE cannot be optimized in parallel. \n Ignoring parallel arguement. \n If you want to optimize in parallel, use least cost paths and gdistance.",immediate. = TRUE)}
      if(is.null(CS.inputs$sub)){
    if (GA.inputs$surface.type[i]=='cat'){
      cnt1 <- cnt1+1    
      names(r)<-GA.inputs$layer.names[i]
      
      single.GA <-ga(type= "real-valued",
                     fitness=Resistance.Opt_single,
                     Resistance=r, 
                     population = GA.inputs$population,
                     selection = GA.inputs$selection,
                     pcrossover=GA.inputs$pcrossover,
                     pmutation=GA.inputs$pmutation,
                     crossover=GA.inputs$crossover,
                     Min.Max=GA.inputs$Min.Max,
                     GA.inputs=GA.inputs,
                     CS.inputs=CS.inputs, 
                     min=GA.inputs$min.list[[i]],
                     max=GA.inputs$max.list[[i]],
                     popSize=GA.inputs$pop.mult*length(GA.inputs$max.list[[i]]),
                     maxiter=GA.inputs$maxiter,
                     run=GA.inputs$run,
                     keepBest=GA.inputs$keepBest,
                     elitism=GA.inputs$percent.elite, 
                     mutation = GA.inputs$mutation,
                     seed = GA.inputs$seed,
                     iter=i,
                     quiet = GA.inputs$quiet)
      
      single.GA@solution <- single.GA@solution/min(single.GA@solution)
      df <- data.frame(id=unique.rast(r),single.GA@solution) 
      r <-subs(r,df)
      names(r)<-GA.inputs$layer.names[i]
      
      Run_CS(CS.inputs,GA.inputs,r,EXPORT.dir=GA.inputs$Results.dir)
      
      Diagnostic.Plots(resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),
                       genetic.dist=CS.inputs$response,
                       plot.dir=GA.inputs$Plots.dir,
                       type="categorical",
                       ID = CS.inputs$ID, 
                       ZZ = CS.inputs$ZZ )
  
   
      RS <- data.frame(GA.inputs$layer.names[i], -single.GA@fitnessValue,single.GA@solution)
      k=GA.inputs$parm.type$n.parm[i]
      Features <- matrix()
      for(z in 1:(k)){
        feature <- paste0("Feature",z)
        Features[z]<-feature
      }
      
      colnames(RS)<-c("Surface","AICc", Features)
      
      RESULTS.cat[[cnt1]]<-RS
     
    } else {   # Processing of continuous surfaces   
      cnt2 <- cnt2+1    
      r<-SCALE(r,0,10)
      names(r)<-GA.inputs$layer.names[i]
      
      single.GA <-ga(type= "real-valued",
                     fitness=Resistance.Opt_single,
                     Resistance=r, 
                     population = GA.inputs$population,
                     selection = GA.inputs$selection,
                     pcrossover=GA.inputs$pcrossover,
                     pmutation=GA.inputs$pmutation,
                     crossover=GA.inputs$crossover,
                     Min.Max=GA.inputs$Min.Max,
                     GA.inputs=GA.inputs,
                     CS.inputs=CS.inputs, 
                     min=GA.inputs$min.list[[i]],
                     max=GA.inputs$max.list[[i]],
                     popSize=GA.inputs$pop.mult*length(GA.inputs$max.list[[i]]),
                     maxiter=GA.inputs$maxiter,
                     run=GA.inputs$run,
                     keepBest=GA.inputs$keepBest,
                     elitism=GA.inputs$percent.elite, 
                     mutation = GA.inputs$mutation,
                     seed = GA.inputs$seed,
                     iter=i,
                     quiet = GA.inputs$quiet) 
      
      # Using GA results, optimize with nlm  
      start.vals <- single.GA@solution[-1]
      
     if(nlm==TRUE){
      # Informed start values; these are the optimized values from the single parameter optimization
      EQ <-get.EQ(single.GA@solution[1])
      Optim.nlm <-nlm(Resistance.Optimization_cont.nlm, log(start.vals), Resistance=r, equation=single.GA@solution[1],get.best=FALSE,CS.inputs=CS.inputs,Min.Max='min',write.dir=GA.inputs$Write.dir)
      
      OPTIM <- Resistance.Optimization_cont.nlm(PARM=(Optim.nlm$estimate),Resistance=r, equation=single.GA@solution[1],get.best=TRUE,CS.inputs=CS.inputs,Min.Max='min',write.dir=GA.inputs$Results.dir)
      
      Diagnostic.Plots(resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),genetic.dist=CS.inputs$response,plot.dir=GA.inputs$Plots.dir,type="continuous",ID = CS.inputs$ID, ZZ = CS.inputs$ZZ)
      
      Plot.trans(PARM=exp(Optim.nlm$estimate), Resistance=GA.inputs$Resistance.stack[[i]], transformation=EQ, print.dir=GA.inputs$Plots.dir,Name=GA.inputs$layer.names[i])
      
      RS<-data.frame(GA.inputs$layer.names[i],Optim.nlm$minimum,EQ,Cont.Param(exp(Optim.nlm$estimate)))
      colnames(RS) <- c("Surface","AICc","Equation","shape","max")
      RESULTS.cont[[cnt2]] <- RS
     
      } else {
        EQ <-get.EQ(single.GA@solution[1])
        r.tran <- Resistance.tran(transformation=single.GA@solution[1],shape=single.GA@solution[2],max=single.GA@solution[3],r=r) 
        names(r.tran)<-GA.inputs$layer.names[i]
      
        Run_CS(CS.inputs,GA.inputs,r.tran,EXPORT.dir=GA.inputs$Results.dir)
      
        Diagnostic.Plots(resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),genetic.dist=CS.inputs$response,plot.dir=GA.inputs$Plots.dir,type="continuous",ID = CS.inputs$ID, ZZ = CS.inputs$ZZ)
        
        Plot.trans(PARM=single.GA@solution[-1], 
                   Resistance=GA.inputs$Resistance.stack[[i]], 
                   transformation=EQ, 
                   print.dir=GA.inputs$Plots.dir)
          
        RS <- data.frame(GA.inputs$layer.names[i], -single.GA@fitnessValue,get.EQ(single.GA@solution[1]),single.GA@solution[2],single.GA@solution[3])
      colnames(RS) <- c("Surface","AICc","Equation","shape","max")
      RESULTS.cont[[cnt2]] <- RS
      }     
    } # Close if-else  
    if (dist_mod==TRUE){
      r <- reclassify(r, c(-Inf,Inf, 1))
      names(r)<-"dist"
      Run_CS(CS.inputs,GA.inputs,r)
      Dist.AIC <- AIC(MLPE.lmm(resistance=paste0(GA.inputs$Write.dir,"dist_resistances.out"),
                               pairwise.genetic=CS.inputs$response,
                               REML=FALSE,
                               ID=CS.inputs$ID,
                               ZZ=CS.inputs$ZZ))
      k<-2
  AICc <- (Dist.AIC)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
  Dist.AICc<-data.frame("Distance", AICc); colnames(Dist.AICc)<-c("Surface","AICc")      
    }
  
  if(null_mod==TRUE){
    response=CS.inputs$response
        
    dat<-data.frame(CS.inputs$ID,response=CS.inputs$response)
    colnames(dat) <- c("pop1", "pop2","response")
        
    # Fit model
    mod <- lFormula(response ~ 1 + (1|pop1), data=dat,REML=FALSE)
    mod$reTrms$Zt <- CS.inputs$ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    Null.AIC <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr)) 
    k<-1
    AICc <- (Null.AIC)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
    Null.AICc<-data.frame("Null", AICc); colnames(Null.AICc)<-c("Surface","AICc")  
     }
  
      } # Close no-sublandscape
  
  #!#!#!#! SUB-LANDSCAPE IMPLEMENTATION !#!#!#!#!
  if(!is.null(CS.inputs$sub)){
    if (GA.inputs$surface.type[i]=='cat'){
      cnt1 <- cnt1+1    
      names(r)<-GA.inputs$layer.names[i]
      
      single.GA <-ga(type= "real-valued",
                     fitness=Resistance.Opt_single_sub,
                     Resistance=r, 
                     population = GA.inputs$population,
                     selection = GA.inputs$selection,
                     pcrossover=GA.inputs$pcrossover,
                     pmutation=GA.inputs$pmutation,
                     crossover=GA.inputs$crossover,
                     Min.Max=GA.inputs$Min.Max,
                     GA.inputs=GA.inputs,
                     CS.inputs=CS.inputs, 
                     min=GA.inputs$min.list[[i]],
                     max=GA.inputs$max.list[[i]],
                     popSize=GA.inputs$pop.mult*length(GA.inputs$max.list[[i]]),
                     maxiter=GA.inputs$maxiter,
                     run=GA.inputs$run,
                     keepBest=GA.inputs$keepBest,
                     elitism=GA.inputs$percent.elite, 
                     mutation = GA.inputs$mutation,
                     seed = GA.inputs$seed,
                     iter=i,
                     quiet = GA.inputs$quiet)
      
      single.GA@solution <- single.GA@solution/min(single.GA@solution)
      df <- data.frame(id=unique.rast(r),single.GA@solution) 
      r <-subs(r,df)
      names(r)<-GA.inputs$layer.names[i]
      
      Run_CS(CS.inputs,GA.inputs,r,EXPORT.dir=GA.inputs$Results.dir)
      
      Diagnostic.Plots(resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),
                       genetic.dist=CS.inputs$response,
                       plot.dir=GA.inputs$Plots.dir,
                       type="categorical",
                       ID = CS.inputs$ID, 
                       ZZ = CS.inputs$ZZ,
                       sublandscape = CS.inputs$sub)
      
      
      RS <- data.frame(GA.inputs$layer.names[i], -single.GA@fitnessValue,single.GA@solution)
      k=GA.inputs$parm.type$n.parm[i]
      Features <- matrix()
      for(z in 1:(k)){
        feature <- paste0("Feature",z)
        Features[z]<-feature
      }
      
      colnames(RS)<-c("Surface","AICc", Features)
      
      RESULTS.cat[[cnt1]]<-RS
      
    } else {   # Processing of continuous surfaces   
      cnt2 <- cnt2+1    
      r<-SCALE(r,0,10)
      names(r)<-GA.inputs$layer.names[i]
      
      single.GA <-ga(type= "real-valued",
                     fitness=Resistance.Opt_single_sub,
                     Resistance=r, 
                     population = GA.inputs$population,
                     selection = GA.inputs$selection,
                     pcrossover=GA.inputs$pcrossover,
                     pmutation=GA.inputs$pmutation,
                     crossover=GA.inputs$crossover,
                     Min.Max=GA.inputs$Min.Max,
                     GA.inputs=GA.inputs,
                     CS.inputs=CS.inputs, 
                     min=GA.inputs$min.list[[i]],
                     max=GA.inputs$max.list[[i]],
                     popSize=GA.inputs$pop.mult*length(GA.inputs$max.list[[i]]),
                     maxiter=GA.inputs$maxiter,
                     run=GA.inputs$run,
                     keepBest=GA.inputs$keepBest,
                     elitism=GA.inputs$percent.elite, 
                     mutation = GA.inputs$mutation,
                     seed = GA.inputs$seed,
                     iter=i,
                     quiet = GA.inputs$quiet) 
      
      # Using GA results, optimize with nlm  
      start.vals <- single.GA@solution[-1]
      
      if(nlm==TRUE)stop("nlm optimization cannot be used with sub-landscape analysis at this time. Please set nlm=FALSE and re-run anallysis.")
      
        EQ <-get.EQ(single.GA@solution[1])
        r.tran <- Resistance.tran(transformation=single.GA@solution[1],shape=single.GA@solution[2],max=single.GA@solution[3],r=r) 
        names(r.tran)<-GA.inputs$layer.names[i]
        
        Run_CS(CS.inputs,GA.inputs,r.tran,EXPORT.dir=GA.inputs$Results.dir)
        
        Diagnostic.Plots(resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),
                         genetic.dist=CS.inputs$response,
                         plot.dir=GA.inputs$Plots.dir,
                         type="continuous",
                         ID = CS.inputs$ID, 
                         ZZ = CS.inputs$ZZ,
                         sublandscape=CS.inputs$sub)
        
        Plot.trans(PARM=single.GA@solution[-1], 
                   Resistance=GA.inputs$Resistance.stack[[i]], 
                   transformation=EQ, 
                   print.dir=GA.inputs$Plots.dir)
        
        RS <- data.frame(GA.inputs$layer.names[i], -single.GA@fitnessValue,get.EQ(single.GA@solution[1]),single.GA@solution[2],single.GA@solution[3])
        colnames(RS) <- c("Surface","AICc","Equation","shape","max")
        RESULTS.cont[[cnt2]] <- RS
      }     
    } # Close if-else  
  
    if (dist_mod==TRUE){
      r <- reclassify(r, c(-Inf,Inf, 1))
      names(r)<-"dist"
      
      dist.resist <- Run_CS(CS.inputs,GA.inputs,r)
      
      
      Dist.AIC <- AIC(MLPE.lmm.sub(resistance=dist.resist,
                               response=CS.inputs$response,
                               REML=FALSE,
                               ID=CS.inputs$ID,
                               ZZ=CS.inputs$ZZ,
                               sub = CS.inputs$sub))
      k<-2
      AICc <- (Dist.AIC)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
      Dist.AICc<-data.frame("Distance", AICc); colnames(Dist.AICc)<-c("Surface","AICc")      
    }
    
    if(null_mod==TRUE){
      response=CS.inputs$response
      
      dat<-data.frame(CS.inputs$ID,response=CS.inputs$response,sub=CS.inputs$sub)
      colnames(dat) <- c("pop1", "pop2","response","sub")
      
      # Fit model
      mod <- lFormula(response ~ 1 + (1|sub/pop1), data=dat,REML=FALSE)
      mod$reTrms$Zt <- CS.inputs$ZZ
      dfun <- do.call(mkLmerDevfun,mod)
      opt <- optimizeLmer(dfun)
      Null.AIC <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr)) 
      k<-1
      AICc <- (Null.AIC)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
      Null.AICc<-data.frame("Null", AICc); colnames(Null.AICc)<-c("Surface","AICc")  
    }
    
  } # Close sublandscape
  
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!###
if(!is.null(gdist.inputs)){
    if (GA.inputs$surface.type[i]=='cat'){
      cnt1 <- cnt1+1    
      names(r)<-GA.inputs$layer.names[i]
      
      single.GA <-ga(type= "real-valued",
                     fitness=Resistance.Opt_single,
                     Resistance=r, 
                     population = GA.inputs$population,
                     selection = GA.inputs$selection,
                     pcrossover=GA.inputs$pcrossover,
                     pmutation=GA.inputs$pmutation,
                     crossover=GA.inputs$crossover,
                     Min.Max=GA.inputs$Min.Max,
                     GA.inputs=GA.inputs,
                     gdist.inputs=gdist.inputs, 
                     min=GA.inputs$min.list[[i]],
                     max=GA.inputs$max.list[[i]],
                     parallel = GA.inputs$parallel,
                     popSize=GA.inputs$pop.mult*length(GA.inputs$max.list[[i]]),
                     maxiter=GA.inputs$maxiter,
                     run=GA.inputs$run,
                     keepBest=GA.inputs$keepBest,
                     elitism=GA.inputs$percent.elite, 
                     mutation = GA.inputs$mutation,
                     seed = GA.inputs$seed,
                     iter=i,
                     quiet = GA.inputs$quiet)
      
      single.GA@solution <- single.GA@solution/min(single.GA@solution)
      df <- data.frame(id=unique.rast(r),single.GA@solution) 
      r <-subs(r,df)
      NAME<-GA.inputs$layer.names[i]
      names(r)<-NAME
      
      cd <- Run_gdistance(gdist.inputs,r)
      save(cd,file=paste0(GA.inputs$Write.dir,NAME,".rda"))
      writeRaster(r,paste0(GA.inputs$Results.dir,NAME,".asc"), overwrite=TRUE)      
      Diagnostic.Plots(resistance.mat=cd,genetic.dist=gdist.inputs$response,plot.dir=GA.inputs$Plots.dir,type="categorical", name=NAME,ID = gdist.inputs$ID, ZZ = gdist.inputs$ZZ)
      
      RS <- data.frame(GA.inputs$layer.names[i], -single.GA@fitnessValue,single.GA@solution)
      k=GA.inputs$parm.type$n.parm[i]
      Features <- matrix()
      for(z in 1:(k)){
        feature <- paste0("Feature",z)
        Features[z]<-feature
      }
      
      colnames(RS)<-c("Surface","AICc", Features)
      
      RESULTS.cat[[cnt1]]<-RS
      
    } else {   # Processing of continuous surfaces   
      cnt2 <- cnt2+1    
      r<-SCALE(r,0,10)
      names(r)<-GA.inputs$layer.names[i]
      
      single.GA <-ga(type= "real-valued",
                     fitness=Resistance.Opt_single,
                     Resistance=r, 
                     population = GA.inputs$population,
                     selection = GA.inputs$selection,
                     pcrossover=GA.inputs$pcrossover,
                     pmutation=GA.inputs$pmutation,
                     crossover=GA.inputs$crossover,
                     Min.Max=GA.inputs$Min.Max,
                     GA.inputs=GA.inputs,
                     gdist.inputs=gdist.inputs, 
                     min=GA.inputs$min.list[[i]],
                     max=GA.inputs$max.list[[i]],
                     parallel = GA.inputs$parallel,
                     popSize=GA.inputs$pop.mult*length(GA.inputs$max.list[[i]]),
                     maxiter=GA.inputs$maxiter,
                     run=GA.inputs$run,
                     keepBest=GA.inputs$keepBest,
                     elitism=GA.inputs$percent.elite, 
                     mutation = GA.inputs$mutation,
                     seed = GA.inputs$seed,
                     iter=i,
                     quiet = GA.inputs$quiet) 
      
      # Using GA results, optimize with nlm  
      start.vals <- single.GA@solution[-1]
      
      if(nlm==TRUE){
        # Informed start values; these are the optimized values from the single parameter optimization
        EQ <-get.EQ(single.GA@solution[1])
        Optim.nlm <-nlm(Resistance.Optimization_cont.nlm, log(start.vals), Resistance=r, equation=single.GA@solution[1],get.best=FALSE,gdist.inputs=gdist.inputs,Min.Max='min',write.dir=GA.inputs$Write.dir)
        
        r<-Resistance.tran(transformation=EQ,shape=Optim.nlm$estimate[1],max=Optim.nlm$estimate[2],r=r,out=GA.inputs$Results.dir)
        
        names(r)<-GA.inputs$layer.names[i]        
        NAME<-GA.inputs$layer.names[i]        
        
        cd <- Run_gdistance(gdist.inputs,r)
        save(cd,file=paste0(GA.inputs$Write.dir,NAME,".rda"))
        writeRaster(r,paste0(GA.inputs$Results.dir,NAME,".asc"), overwrite=TRUE)  
        
        Diagnostic.Plots(resistance.mat=cd,genetic.dist=gdist.inputs$response,plot.dir=GA.inputs$Plots.dir,type="continuous",name=NAME,ID = gdist.inputs$ID, ZZ = gdist.inputs$ZZ)
                
        Plot.trans(PARM=exp(Optim.nlm$estimate), Resistance=GA.inputs$Resistance.stack[[i]], transformation=EQ, print.dir=GA.inputs$Plots.dir,Name=GA.inputs$layer.names[i])
        
        RS<-data.frame(GA.inputs$layer.names[i],Optim.nlm$minimum,EQ,Cont.Param(exp(Optim.nlm$estimate)))
        colnames(RS) <- c("Surface","AICc","Equation","shape","max")
        RESULTS.cont[[cnt2]] <- RS
        
      } else {
        EQ <-get.EQ(single.GA@solution[1])
        r <- Resistance.tran(transformation=single.GA@solution[1],shape=single.GA@solution[2],max=single.GA@solution[3],r=r) 
        names(r)<-GA.inputs$layer.names[i]        
        NAME<-GA.inputs$layer.names[i]        
        
        cd <- Run_gdistance(gdist.inputs,r)
        save(cd,file=paste0(GA.inputs$Write.dir,NAME,".rda"))
        writeRaster(r,paste0(GA.inputs$Results.dir,NAME,".asc"), overwrite=TRUE)      
  
        
        Diagnostic.Plots(resistance.mat=cd,genetic.dist=gdist.inputs$response,plot.dir=GA.inputs$Plots.dir,type="continuous",name=NAME,ID = gdist.inputs$ID, ZZ = gdist.inputs$ZZ)
        
        Plot.trans(PARM=single.GA@solution[-1], 
                   Resistance=GA.inputs$Resistance.stack[[i]], 
                   transformation=EQ, 
                   print.dir=GA.inputs$Plots.dir)
        
        RS <- data.frame(GA.inputs$layer.names[i], -single.GA@fitnessValue,get.EQ(single.GA@solution[1]),single.GA@solution[2],single.GA@solution[3])
        colnames(RS) <- c("Surface","AICc","Equation","shape","max")
        RESULTS.cont[[cnt2]] <- RS
      }     
    } # Close if-else  

    if (dist_mod==TRUE){
      r <- reclassify(r, c(-Inf,Inf, 1))
      names(r)<-"dist"  
      cd <- Run_gdistance(gdist.inputs,r)
      
      Dist.AIC <- suppressWarnings(AIC(MLPE.lmm2(resistance=cd,
                                  response=gdist.inputs$response,
                                  ID=gdist.inputs$ID,
                                  ZZ=gdist.inputs$ZZ,
                                  REML=FALSE)))
      ROW <- nrow(gdist.inputs$ID)
      k<-2
      AICc <- (Dist.AIC)+(((2*k)*(k+1))/(ROW-k-1))
      Dist.AICc<-data.frame("Distance", AICc); colnames(Dist.AICc)<-c("Surface","AICc")      
    }
    
    if(null_mod==TRUE){          
      dat<-data.frame(gdist.inputs$ID,response=gdist.inputs$response)
      colnames(dat) <- c("pop1", "pop2","response")
      
      # Fit model
      mod <- lFormula(response ~ 1 + (1|pop1), data=dat,REML=FALSE)
      mod$reTrms$Zt <- gdist.inputs$ZZ
      dfun <- do.call(mkLmerDevfun,mod)
      opt <- optimizeLmer(dfun)
      Null.AIC <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr)) 
      ROW <- nrow(gdist.inputs$ID)
      k<-1
      AICc <- (Null.AIC)+(((2*k)*(k+1))/(ROW-k-1))
      Null.AICc<-data.frame("Null", AICc); colnames(Null.AICc)<-c("Surface","AICc")  
    }
  }
} # Close ascii loop
  
  #!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##
  # Make results data frame
  Results.cat<-data.frame()
  Results.cont<-data.frame()
  # cnt1<-0
  # cnt2<-0
  for (i in 1:GA.inputs$n.layers){
    if(GA.inputs$surface.type[i]=='cat'){
      #     cnt1 <- cnt1+1
      #     RS <- data.frame(GA.inputs$layer.names[i], -(RESULTS.cat[[i]]@fitnessValue),RESULTS[[i]]@solution)
      Results.cat <- do.call(rbind.fill,RESULTS.cat)
    } else {
      #   cnt2 <-cnt2+1
      #   RS <- data.frame(GA.inputs$layer.names[i], -(RESULTS.cont[[i]]@fitnessValue), Cont.Param(RESULTS[[i]]@solution))
      Results.cont <- do.call(rbind,RESULTS.cont)
    }
  }
  #!##!##!##!##!##!##!##!##!##!##!##
  # Compile results into tables
  cat("\n")
  cat("\n")
  if(nrow(Results.cat)>0){
  Features <- array()
  for(i in 1:ncol(Results.cat)-2){
    feature <- paste0("Feature",i)
    Features[i]<-feature
  }
  colnames(Results.cat)<-c("Surface","AICc", Features)
  Results.cat <-  Results.cat[order(Results.cat$AICc),]
  write.table(Results.cat,paste0(GA.inputs$Results.dir,"CategoricalResults.csv"),sep=",",col.names=T,row.names=F)
  }
  
  if(ncol(Results.cont)>0){    
    colnames(Results.cont)<-c("Surface","AICc","Equation","shape","max")
    Results.cont <- Results.cont[order(Results.cont$AICc),]
    write.table(Results.cont,paste0(GA.inputs$Results.dir,"ContinuousResults.csv"),sep=",",col.names=T,row.names=F)
  }
  
  # Full Results
  if(nrow(Results.cat)>0 & nrow(Results.cont)>0){
    Results.All<-rbind(Results.cat[,c(1,2)],Results.cont[,c(1,2)])
  } else if(nrow(Results.cat)<1 & nrow(Results.cont)>0){
    Results.All<-(Results.cont[,c(1,2)])  
  } else {
    Results.All<-(Results.cat[,c(1,2)])    
  }
  
  if(dist_mod==TRUE) Results.All<-rbind(Results.All,Dist.AICc)
  if(null_mod==TRUE) Results.All<-rbind(Results.All,Null.AICc)

  Results.All <- Results.All[order(Results.All$AICc),]

  cat("\n")
  cat("\n")
  write.table(Results.All,paste0(GA.inputs$Results.dir,"All_Results_AICc.csv"),sep=",",col.names=T,row.names=F)
  
  # Get parameter estimates
if(!is.null(CS.inputs)){  
  MLPE.results<-MLPE.lmm_coef(resistance=GA.inputs$Results.dir,
                              genetic.dist=CS.inputs$response,
                              out.dir=GA.inputs$Results.dir,
                              method="cs",
                              ID=CS.inputs$ID,
                              ZZ=CS.inputs$ZZ,
                              sub=CS.inputs$sub)

  
} else {  
  MLPE.results<-MLPE.lmm_coef(resistance=GA.inputs$Write.dir,
                              genetic.dist=gdist.inputs$response,
                              out.dir=GA.inputs$Results.dir,
                              method="gd",
                              ID=gdist.inputs$ID,
                              ZZ=gdist.inputs$ZZ)  
}

  rt<-proc.time()[3]-t1
  # Full Results
  if(nrow(Results.cat)>0 & nrow(Results.cont)>0){
    RESULTS<-list(ContinuousResults=Results.cont, CategoricalResults=Results.cat,AICc=Results.All,MLPE=MLPE.results, Run.Time=rt)
  } else if(nrow(Results.cat)<1 & nrow(Results.cont)>0){
    RESULTS<-list(ContinuousResults=Results.cont, CategoricalResults=NULL,AICc=Results.All,MLPE=MLPE.results, Run.Time=rt)
  } else if(nrow(Results.cat)>0 & nrow(Results.cont)<1){
    RESULTS<-list(ContinuousResults=NULL, CategoricalResults=Results.cat,AICc=Results.All,MLPE=MLPE.results, Run.Time=rt)
  } else {    
    RESULTS<-list(ContinuousResults=NULL, CategoricalResults=NULL,AICc=Results.All,MLPE=MLPE.results, Run.Time=rt)
  }

file.remove(list.files(GA.inputs$Write.dir,full.names=TRUE))
return(RESULTS)
  #!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#
}


#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# Single command function to execute mulit surface optimization #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#' Simultaneous optimization of multiple resistance surfaces
#' 
#' Optimize multiple resistance surfaces simultaneously using genetic algorithms
#' 
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @return This function optimizes multiple resistance surfaces, returning a Genetic Algorithm (GA) object with summary information. Diagnostic plots of model fit are output to the "Results/Plots" folder that is automatically generated within the folder containing the optimized ASCII files. A text summary of the optimization settings and results is printed to the results folder.
#' @usage MS_optim(CS.inputs, gdist.inputs, GA.inputs)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
MS_optim<-function(CS.inputs=NULL, gdist.inputs=NULL, GA.inputs){
  if(!is.null(CS.inputs)){
    if(GA.inputs$parallel!=FALSE) {warning("\n CIRCUITSCAPE cannot be optimized in parallel. \n Ignoring parallel arguement. \n If you want to optimize in parallel, use least cost paths and gdistance.",immediate. = TRUE)}
    if(is.null(CS.inputs$pairs_to_include)){
    t1<-proc.time()[3]
    multi.GA_nG <-ga(type= "real-valued",
                   fitness=Resistance.Opt_multi,
                   population=GA.inputs$population,
                   selection=GA.inputs$selection,
                   mutation=GA.inputs$mutation,
                   pcrossover=GA.inputs$pcrossover,
                   crossover=GA.inputs$crossover,
                   pmutation=GA.inputs$pmutation,
                   Min.Max=GA.inputs$Min.Max,
                   GA.inputs=GA.inputs,
                   CS.inputs=CS.inputs,
                   min=GA.inputs$ga.min,
                   max=GA.inputs$ga.max,
                   popSize=GA.inputs$pop.size,
                   maxiter=GA.inputs$maxiter,
                   run=GA.inputs$run,
                   parallel = FALSE,
                   keepBest=GA.inputs$keepBest,
                   seed = GA.inputs$seed,
                   suggestions=GA.inputs$SUGGESTS,
                   quiet = GA.inputs$quiet) 
      rt<-proc.time()[3]-t1
   
    Opt.parm <- GA.opt <- multi.GA_nG@solution
    for(i in 1:GA.inputs$n.layers){
      if(GA.inputs$surface.type[i]=="cat"){
        ga.p <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
        parm <- ga.p/min(ga.p)
        Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
              
      } else {
        parm <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]        
        Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
      }
    }
   } else { # Start multi-surface sub-landscape
     t1<-proc.time()[3]
     multi.GA_nG <-ga(type= "real-valued",
                      fitness=Resistance.Opt_multi_sub,
                      population=GA.inputs$population,
                      selection=GA.inputs$selection,
                      mutation=GA.inputs$mutation,
                      pcrossover=GA.inputs$pcrossover,
                      crossover=GA.inputs$crossover,
                      pmutation=GA.inputs$pmutation,
                      Min.Max=GA.inputs$Min.Max,
                      GA.inputs=GA.inputs,
                      CS.inputs=CS.inputs,
                      min=GA.inputs$ga.min,
                      max=GA.inputs$ga.max,
                      popSize=GA.inputs$pop.size,
                      maxiter=GA.inputs$maxiter,
                      run=GA.inputs$run,
                      parallel = FALSE,
                      keepBest=GA.inputs$keepBest,
                      seed = GA.inputs$seed,
                      suggestions=GA.inputs$SUGGESTS,
                      quiet = GA.inputs$quiet) 
     rt<-proc.time()[3]-t1
     
     Opt.parm <- GA.opt <- multi.GA_nG@solution
     for(i in 1:GA.inputs$n.layers){
       if(GA.inputs$surface.type[i]=="cat"){
         ga.p <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
         parm <- ga.p/min(ga.p)
         Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
         
       } else {
         parm <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]        
         Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
       } # End if-else
     } # End layer loop
     
   } # End multi-surface sub-landscape
   
  } # End multi-surface CS optimization
   
   
    multi.GA_nG@solution <- Opt.parm
    
  
  RAST<-Combine_Surfaces(PARM=multi.GA_nG@solution,CS.inputs=CS.inputs,GA.inputs=GA.inputs, rescale = TRUE)
  NAME<-paste(GA.inputs$parm.type$name,collapse=".")
  names(RAST)<-NAME
  Run_CS(CS.inputs,GA.inputs,r=RAST,CurrentMap=FALSE,EXPORT.dir=GA.inputs$Results.dir)
  
  ifelse(length(unique.rast(RAST))>15,type<-"continuous", type<-"categorical")
  
  if(is.null(CS.inputs$sub)){
  
  Diagnostic.Plots(resistance.mat=paste0(GA.inputs$Results.dir,NAME,"_resistances.out"),genetic.dist=CS.inputs$response,plot.dir=GA.inputs$Plots.dir,type=type,ID = CS.inputs$ID, ZZ = CS.inputs$ZZ)
  
  # Get parameter estimates
  MLPE.results<-MLPE.lmm_coef(resistance=GA.inputs$Results.dir,
                              genetic.dist=CS.inputs$response,
                              out.dir=GA.inputs$Results.dir,
                              method="cs",
                              ID = CS.inputs$ID, 
                              ZZ = CS.inputs$ZZ)  
  
  Result.txt(GA.results=multi.GA_nG,GA.inputs=GA.inputs, method="CIRCUITSCAPE", Run.Time=rt) 
  file.remove(list.files(GA.inputs$Write.dir,full.names=TRUE))
  return(multi.GA_nG)
    } else { # Diagnostic plots for sub-landscape
      Diagnostic.Plots(resistance.mat=paste0(GA.inputs$Results.dir, NAME,"_resistances.out"),
                       genetic.dist=CS.inputs$response,
                       plot.dir=GA.inputs$Plots.dir,
                       type=type,
                       ID = CS.inputs$ID, 
                       ZZ = CS.inputs$ZZ,
                       sublandscape=CS.inputs$sub)
      
      # Get parameter estimates
      MLPE.results<-MLPE.lmm_coef(resistance=GA.inputs$Results.dir,
                                  genetic.dist=CS.inputs$response,
                                  out.dir=GA.inputs$Results.dir,
                                  method="cs",
                                  ID = CS.inputs$ID, 
                                  ZZ = CS.inputs$ZZ)  
      
      Result.txt(GA.results=multi.GA_nG,GA.inputs=GA.inputs, method="CIRCUITSCAPE", Run.Time=rt) 
      file.remove(list.files(GA.inputs$Write.dir,full.names=TRUE))
      return(multi.GA_nG)
  }

 # GDIST code 
  if(!is.null(gdist.inputs)){
        t1<-proc.time()[3]
      multi.GA_nG <-ga(type= "real-valued",
                     fitness=Resistance.Opt_multi,
                     population=GA.inputs$population,
                     selection=GA.inputs$selection,
                     mutation=GA.inputs$mutation,
                     pcrossover=GA.inputs$pcrossover,
                     crossover=GA.inputs$crossover,
                     pmutation=GA.inputs$pmutation,
                     Min.Max=GA.inputs$Min.Max,
                     GA.inputs=GA.inputs,
                     gdist.inputs=gdist.inputs,
                     min=GA.inputs$ga.min,
                     max=GA.inputs$ga.max,
                     popSize=GA.inputs$pop.size,
                     maxiter=GA.inputs$maxiter,
                     parallel = GA.inputs$parallel,
                     run=GA.inputs$run,
                     keepBest=GA.inputs$keepBest,
                     seed = GA.inputs$seed,
                     suggestions=GA.inputs$SUGGESTS,
                     quiet = GA.inputs$quiet) 
       rt<-proc.time()[3]-t1
      
      Opt.parm <- GA.opt <- multi.GA_nG@solution
      for(i in 1:GA.inputs$n.layers){
        if(GA.inputs$surface.type[i]=="cat"){
          ga.p <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
          parm <- ga.p/min(ga.p)
          Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
          
        } else {
          parm <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]          
          Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
        }
      }
      multi.GA_nG@solution <- Opt.parm
      
    RAST<-Combine_Surfaces(PARM=multi.GA_nG@solution,gdist.inputs=gdist.inputs,GA.inputs=GA.inputs, rescale = TRUE)
    NAME<-paste(GA.inputs$parm.type$name,collapse=".")
    names(RAST)<-NAME
    cd <- Run_gdistance(gdist.inputs,RAST)
    save(cd,file=paste0(GA.inputs$Results.dir,NAME,".rda"))
    writeRaster(RAST,paste0(GA.inputs$Results.dir,NAME,".asc"), overwrite=TRUE)
    
    ifelse(length(unique.rast(RAST))>15,type<-"continuous", type<-"categorical")

#     Run_CS(CS.inputs,GA.inputs,r=RAST,CurrentMap=FALSE,EXPORT.dir=GA.inputs$Results.dir)
    
    Diagnostic.Plots(resistance.mat=cd,genetic.dist=gdist.inputs$response,plot.dir=GA.inputs$Plots.dir,type=type, name=NAME,ID = gdist.inputs$ID, ZZ = gdist.inputs$ZZ)
    
    # Get parameter estimates
    MLPE.results<-MLPE.lmm_coef(resistance=GA.inputs$Results.dir,
                                genetic.dist=gdist.inputs$response, 
                                out.dir=GA.inputs$Results.dir, 
                                method="gd",
                                ID=gdist.inputs$ID,
                                ZZ=gdist.inputs$ZZ)  
    
    Result.txt(GA.results=multi.GA_nG,GA.inputs=GA.inputs,method="cost distance", Run.Time=rt) 
file.remove(list.files(GA.inputs$Write.dir,full.names=TRUE))
return(multi.GA_nG)
  }

  #!###  RUN BRENT OPTIMIZATION #!##
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
# #!##!##!##!##!##!##!##!##!##!##!##!##!##!##
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
# # #!##
# multi.GA_nG@solution <- Opt.parm
# multi.GA_nG@fitnessValue <- Parm.multiplier$value
# #!##!##!##!##!##!##!##!##!##!##!##
  
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

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
# Max.optim_Brent <- function(PARM,CS.inputs,GA.inputs, Min.Max='min', quiet=FALSE, GA.opt){
#   t1<-Sys.time()
#   
#   ID<-CS.inputs$ID
#   ZZ<-CS.inputs$ZZ
#   response<-CS.inputs$response
#   Opt.parm <-vector(length=length(PARM))
# 
#   CS_Point.File<-CS.inputs$CS_Point.File
#   CS.program<-CS.inputs$CS.program
#   EXPORT.dir<-GA.inputs$Write.dir
#   #!##!#
#   r <- GA.inputs$Resistance.stack
#   
#   for(i in 1:GA.inputs$n.layers){
#     if(GA.inputs$surface.type[i]=="cat"){
#       ga.p <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
#       parm <- ((ga.p-min(ga.p))*PARM[1])+1
#       Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
#       parm<-parm/min(parm)
#       df <- data.frame(id=unique.rast(r[[i]]),parm) # Data frame with original raster values and replacement values
#       r[[i]] <-subs(r[[i]],df)
#       
#       r[[i]]<-r[[i]]-(cellStats(x=r[[i]],stat="min"))
#       
#     } else {
#       rast <-SCALE(data=r[[i]],MIN=0,MAX=10)
#       parm <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
#       mx<-parm[3]*PARM[1]
#       parm[3]<-mx
#       Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
#       
#       # Set equation for continuous surface
#       equation <- floor(parm[1]) # Parameter can range from 1-9.99
#       
#       # Read in resistance surface to be optimized
#       SHAPE <-  (parm[2])
#       Max.SCALE <- (parm[3])
#       
#       rick.eq<-(equation==2||equation==4||equation==6||equation==8)
#       if(rick.eq==TRUE & SHAPE>6){
#         equation<-9
#       }
#       
#       # Apply specified transformation
#       if(equation==1){
#         SIGN=-1 # Inverse
#         R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE))+SIGN # Monomolecular
#         R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
#         R.vec <- rev(R) # Reverse
#         rast.R <- setValues(R,values=R.vec)
#         r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
#         EQ <- "Inverse-Reverse Monomolecular"
#         
#       } else if(equation==5){
#         SIGN=1
#         R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE))+SIGN # Monomolecular
#         R.vec <- rev(R) # Reverse
#         rast.R <- setValues(R,values=R.vec)
#         r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
#         EQ <- "Reverse Monomolecular"        
#         
#       } else if(equation==3){
#         SIGN=1
#         r[[i]] <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE))+SIGN # Monomolecular    
#         r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
#         EQ <- "Monomolecular"
#         
#       } else if (equation==7) {
#         SIGN=-1 #Inverse
#         R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE))+SIGN # Monomolecular
#         r[[i]] <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
#         r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
#         EQ <- "Inverse Monomolecular"        
#         
#       } else if (equation==8) {
#         SIGN=-1 #Inverse
#         R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE))+SIGN # Ricker
#         r[[i]] <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
#         r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
#         EQ <- "Inverse Ricker"  
#         
#       } else if (equation==4) {
#         SIGN=1
#         r[[i]] <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE))+SIGN #  Ricker
#         r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
#         EQ <- "Ricker"
#         
#       } else if (equation==6) {
#         SIGN=1
#         R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE))+SIGN #  Ricker
#         R.vec <- rev(R)
#         rast.R <- setValues(R,values=R.vec)
#         r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
#         EQ <- "Reverse Ricker"        
#         
#       } else if (equation==2) {
#         SIGN=-1 # Inverse
#         R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE))+SIGN # Ricker
#         R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
#         R.vec <- rev(R) # Reverse
#         rast.R <- setValues(R,values=R.vec)
#         r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
#         EQ <- "Inverse-Reverse Ricker"
#         
#       } else {
#         r[[i]] <- (rast*0) #  Cancel layer...set to zero
#       } # End if-else  
#     } # Close parameter type if-else  
#   } # Close layer loop
#   
#   
#   File.name <- "resist_surface"
#   
#   multi_surface <- sum(r)+1 # Add all surfaces together (+1 for distance)
#   if(cellStats(multi_surface,"max")>5e5)  multi_surface<-SCALE(multi_surface,1,5e5) # Rescale surface in case resistance are too high
#   
#   writeRaster(x=multi_surface,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
#   
#   # Modify and write Circuitscape.ini file
#   #!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#  
#   BATCH<-paste0(EXPORT.dir,File.name,".ini")        
#   OUT<-paste0(paste0("output_file = ",EXPORT.dir), File.name,".out")
#   HABITAT<-paste0("habitat_file = ",paste0(EXPORT.dir,File.name,".asc"))
#   LOCATION.FILE <- paste0("point_file = ", CS.inputs$CS_Point.File)
#   ifelse(CS.inputs$Neighbor.Connect==4,connect<-"True",connect<-"False")
#   CONNECTION=paste0("connect_four_neighbors_only=",connect)
#   
#   #   if(CS.version=='3.5.8'){
#   #     write.CS_3.5.8(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,VERSION=VERSION)
#   #   } else {
#   write.CS_4.0(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION)    
#   #   }
#   
#   #!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#
#   # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
#   hidden = TRUE  
#   # Run Circuitscape
#   if(CS.inputs$platform=="pc"){
#     CS.ini <- paste0(EXPORT.dir,File.name,".ini")
#     CS.Run.output<-system(paste(CS.inputs$CS.program, CS.ini), hidden)
#   } else {
#     CS.ini <- paste0(EXPORT.dir,File.name,".ini")
#     CS.Run.output<-system(paste("python", CS.inputs$CS.program, CS.ini), hidden)
#   }
#   
#   
#   #!##!##!##!##!##!##!##!##!##!##!##!##!###
#   # Run mixed effect model on each Circuitscape effective resistance
#   
#   CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")
#   
#   # Get AIC statistic for transformed-scaled resistance surface
#   cs.matrix<-read.matrix(CS.results)
#   cs.matrix<-scale(cs.matrix,center=TRUE,scale=TRUE)
#   # cs.matrix2<-round(read.matrix(CS.results),4)
#   
#   data<-cbind(ID,cs.matrix,response)
#   
#   # Assign value to layer
#   LAYER<-assign("LAYER",value=data$cs.matrix)
#   
#   # Fit model
#   mod <- lFormula(response ~ LAYER + (1|pop1), data=data,REML=FALSE)
#   mod$reTrms$Zt <- ZZ
#   dfun <- do.call(mkLmerDevfun,mod)
#   opt <- optimizeLmer(dfun)
#   AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
#   #    summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
#   
#   k<-max(GA.inputs$parm.index)+1
#   AICc <- (AIC.stat)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
#   
#  
#   t2 <-Sys.time()
#   if(quiet==FALSE){
#   cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
#   cat(paste0("\t", "AICc = ",round(AICc,4)),"\n","\n")
#   }
#   
#   
#   OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
# }

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# Create continuous surface response figures #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!## 
# Response.Figs<- function(Optim.input){
#   Top.params <- read.csv(file=paste0(Optim.input$Results.dir,"/TopModel_Optimization_Parameters.csv"),header=T)
#   dir.create(file.path(Optim.input$Results.dir, "Plots"))
#   
#   PDF.dir <- paste0(Optim.input$Results.dir,"Plots/")
#   for (i in 1:nrow(Top.params)){
#     ASCII.file <- list.files(Optim.input$ASCII.dir,pattern=paste0(Top.params[i,1],".asc"),full.names=TRUE)
#     if(Top.params[i,5]<1e-5 | Top.params[i,6]<1e-5) {
#       cat("\n","\n",paste0("Plotting of ", Top.params[i,1]," could not be completed due to extremely small parameter estimates."),"\n")
#       next # Use 'next' command to avoid testing invalid combinations    
#     } else {
#       PLOT.response(PARM=Top.params[i,c(5,6)],Resistance=ASCII.file,equation=Top.params[i,2],AIC=Top.params[i,7], OutputFolder=PDF.dir)
#     }
#   }
# }

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# RUN CIRCUITSCAPE #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#' Run CIRCUITSCAPE in R
#' 
#' Execute CS from R
#' 
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param r Accepts two types of inputs. Provide either the path to the resistance surface file (.asc) or specify an R RasterLayer object
#' @param CurrentMap Logical. If TRUE, the cumulative current resistance map will be generated during the CS run (Default = FALSE)
#' @param EXPORT.dir Directory where CS results should be written (Default = GA.inputs$Write.dir, which is a temporary directory for reading/writing CS results). It is critical that there are NO SPACES in the directory, as this will cause the function to fail.
#' @param output Specifiy either "matrix" or "raster". "matrix" will return the lower half of the pairwise resistance matrix (default), while "raster" will return a \code{raster} object of the cumulative current map. The raster map can only be returned if \code{CurrentMap=TRUE}
#' @param hidden Logical. If TRUE (Default), then no output from CIRCUITSCAPE will be printed to the console. Only set to FALSE when trying to troubleshoot/debug code.
#' @return Vector of CIRCUITSCAPE resistance distances (lower half of "XXX_resistances.out"). Alternatively, a raster object of the cumulative current map can be returned when \code{CurrentMap=TRUE} and \code{output="raster"}.
#' @usage Run_CS(CS.inputs, GA.inputs, r, CurrentMap = FALSE, EXPORT.dir = GA.inputs$Write.dir, output = "matrix", hidden = TRUE)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Run_CS <- function(CS.inputs,GA.inputs,r,CurrentMap=FALSE,EXPORT.dir=GA.inputs$Write.dir, output="matrix", hidden = TRUE){
if(class(r)[1]!='RasterLayer') {
  R<-raster(r)
  File.name <- basename(r)
  File.name<-sub(".asc", "", File.name) 
  names(R)<-File.name
} else {
  R=r
  
}
  
  if (CurrentMap==FALSE){
    File.name <- R@data@names
    MAP="write_cum_cur_map_only = False"
    CURRENT.MAP="write_cur_maps = False"
    
   } else {
    File.name <- R@data@names
    MAP="write_cum_cur_map_only = True"
    CURRENT.MAP="write_cur_maps = 1"
  }

  #!##!#
  
  if(cellStats(R,"max")>1e6)  R<-SCALE(R,1,1e6) # Rescale surface in case resistances are too high
  R <- reclassify(R, c(-Inf,0, 1))
  
  #     plot(multi_surface)
  
  writeRaster(x=R,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  
  # Modify and write Circuitscape.ini file
  #!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#  
 ifelse(CS.inputs$Neighbor.Connect==4,connect<-"True",connect<-"False")
 if(is.null(CS.inputs$pairs_to_include)){
   PAIRS_TO_INCLUDE <- paste0("included_pairs_file = (Browse for a file with pairs to include or exclude)")
   PAIRS <- paste0("use_included_pairs = False")   
 } else {
   PAIRS_TO_INCLUDE <- paste0("included_pairs_file = ", CS.inputs$pairs_to_include)
   PAIRS <- paste0("use_included_pairs = True") 
 }

  write.CS_4.0(BATCH=paste0(EXPORT.dir,File.name,".ini"),
               OUT=paste0(paste0("output_file = ",EXPORT.dir), File.name,".out"),
               HABITAT=paste0("habitat_file = ",paste0(EXPORT.dir,File.name,".asc")),
               LOCATION.FILE=paste0("point_file = ", CS.inputs$CS_Point.File),
               CONNECTION=paste0("connect_four_neighbors_only=",connect),
               MAP=MAP,
               CURRENT.MAP=CURRENT.MAP,
               PAIRS_TO_INCLUDE=PAIRS_TO_INCLUDE,
               PAIRS=PAIRS)
  #!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#
# Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
# hidden = TRUE  
# Run Circuitscape
if(CS.inputs$platform=="pc"){
  CS.ini <- paste0(EXPORT.dir,File.name,".ini")
    CS.Run.output<-system(paste(CS.inputs$CS.program, CS.ini), hidden)
} else {
  CS.ini <- paste0(EXPORT.dir,File.name,".ini")
  CS.Run.output<-system(paste("python", CS.inputs$CS.program, CS.ini), hidden)
}
 
  
  #!##!##!##!##!##!##!##!##!##!##!##!##!###
  # Run mixed effect model on each Circuitscape effective resistance
  
  CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")

if(output=="raster" & CurrentMap==TRUE){
  rast<-raster(paste0(EXPORT.dir, File.name,"_cum_curmap.asc"))
  NAME <- basename(rast@file@name)
  NAME<-sub("^([^.]*).*", "\\1", NAME) 
  names(rast)<-NAME
  (rast)
} else {
  (cs.matrix<-read.matrix(CS.results))  
}
}
#!##!##!##!##!##!##!##!##!##!###

Run_CS2 <- function(CS.inputs,GA.inputs,r,EXPORT.dir=GA.inputs$Write.dir, File.name){
    File.name <- File.name
  
  # Modify and write Circuitscape.ini file
  #!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#  
  ifelse(CS.inputs$Neighbor.Connect==4,connect<-"True",connect<-"False")
  
   if(is.null(CS.inputs$pairs_to_include)){
    PAIRS_TO_INCLUDE <- paste0("included_pairs_file = (Browse for a file with pairs to include or exclude)")
    PAIRS <- paste0("use_included_pairs = False")   
 } else {
    PAIRS_TO_INCLUDE <- paste0("included_pairs_file = ", CS.inputs$pairs_to_include)
    PAIRS <- paste0("use_included_pairs = True") 
 }

  write.CS_4.0(BATCH=paste0(EXPORT.dir,File.name,".ini"),
               OUT=paste0(paste0("output_file = ",EXPORT.dir), File.name,".out"),
               HABITAT=paste0("habitat_file = ",paste0(EXPORT.dir,File.name,".asc")),
               LOCATION.FILE=paste0("point_file = ", CS.inputs$CS_Point.File),
               CONNECTION=paste0("connect_four_neighbors_only=",connect),
               MAP="write_cum_cur_map_only = False",
               CURRENT.MAP="write_cur_maps = False",
               PAIRS_TO_INCLUDE=PAIRS_TO_INCLUDE,
               PAIRS=PAIRS)
  #!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#
  # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
  hidden = TRUE  
  # Run Circuitscape
  if(CS.inputs$platform=="pc"){
    CS.ini <- paste0(EXPORT.dir,File.name,".ini")
    CS.Run.output<-system(paste(CS.inputs$CS.program, CS.ini), hidden)
  } else {
    CS.ini <- paste0(EXPORT.dir,File.name,".ini")
    CS.Run.output<-system(paste("python", CS.inputs$CS.program, CS.ini), hidden)
  }
  
  CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")
  
  (cs.matrix<-read.matrix(CS.results))
}
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!###

#' Get cost distance using gdistance
#' 
#' Execute gdistance
#' 
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param r Accepts two types of inputs. Provide either the path to the resistance surface file (.asc) or specify an R RasterLayer object
#' @return A costDistance matrix object from gdistance
#' @usage Run_gdistance(gdist.inputs, r)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Run_gdistance <- function(gdist.inputs, r){
  if(class(r)[1]!='RasterLayer') {
    r<-raster(r)   
  }  
  tr <- transition(x=r, transitionFunction=gdist.inputs$transitionFunction,directions=gdist.inputs$directions)
  if(gdist.inputs$longlat==TRUE|gdist.inputs$directions>=8){
    tr <- geoCorrection(tr, "c")
  }
  ret<- costDistance(tr, gdist.inputs$samples)
  return(ret)
}

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# COMBINE RESISTANCE SURFACES--No GAUSSIAN #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#' Combine multiple resistance surfaces together
#' 
#' Combine multiple resistance surfaces into new composite surface based on specified parameters
#' 
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Requires a vector with parameters specified in the order of resistance surfaces
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function. 
#' @param out Directory to write combined .asc file. Default = NULL and no files are exported
#' @param File.name Name of output .asc file. Default is the combination of all surfaces combined, separated by "."
#' @param rescale Locical. If TRUE (default), the values of the combined raster surface will be divided by the minimum value to create a resistance surface with a minimum value = 1.
#' @details \code{PARM} is designed to accept the output of \code{MS_optim}. For continuous surfaces, there are three terms: 1) Transformation, 2) shape, and 3) maximum value. Transformation must be provided as a numeric value:\cr
#' \tabular{ll}{
#'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#'    \tab 2 = "Inverse-Reverse Ricker"\cr
#'    \tab 3 = "Monomolecular"\cr
#'    \tab 4 = "Ricker"\cr
#'    \tab 5 = "Reverse Monomolecular"\cr
#'    \tab 6 = "Reverse Ricker"\cr
#'    \tab 7 = "Inverse Monomolecular"\cr
#'    \tab 8 = "Inverse Ricker"\cr
#'    \tab 9 = "Distance"\cr
#'    }
#' 
#' The Distance transformation sets all values equal to one. Because of the flexibility of the Ricker function to take a monomolecular shape (try \code{Plot.trans(PARM=c(10,100), Resistance=c(1,10), transformation="Ricker")} to see this), whenever a shape parameter >6 is selected in combination with a Ricker family transformation, the transformation reverts to a Distance transformation. In general, it seems that using a combination of intermediate Ricker and Monomolecular transformations provides the best, most flexible coverasge of parameter space.
#' @return R raster object that is the sum all transformed and/or reclassified resistance surfaces provided
#' @usage Combine_Surfaces(PARM, CS.inputs, gdist.inputs, GA.inputs, out, File.name, rescale)
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Combine_Surfaces <- function(PARM, CS.inputs=NULL, gdist.inputs=NULL, GA.inputs, out=NULL, File.name=paste(GA.inputs$parm.type$name,collapse="."), rescale = TRUE){
  if(!is.null(CS.inputs)){
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  response<-CS.inputs$response
  CS_Point.File<-CS.inputs$CS_Point.File
  CS.program<-CS.inputs$CS.program
  EXPORT.dir<-out
  }
  
  if(!is.null(gdist.inputs)){
    ID<-gdist.inputs$ID
    ZZ<-gdist.inputs$ZZ
    response<-gdist.inputs$response
    samples<-gdist.inputs$samples
    EXPORT.dir<-out
  }

  
  #!##!#
  r <- GA.inputs$Resistance.stack
  
  for(i in 1:GA.inputs$n.layers){
    if(GA.inputs$surface.type[i]=="cat"){
      parm <- PARM[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
      parm <- parm/min(parm)
      df <- data.frame(id=unique.rast(r[[i]]),parm) # Data frame with original raster values and replacement values
      r[[i]] <-subs(r[[i]],df)
      
      r[[i]]<-r[[i]]#-1 # Set minimum to 0  

      
    } else {
      rast <-SCALE(data=r[[i]],MIN=0,MAX=10)
      parm <- PARM[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
      
      
      # Set equation for continuous surface
      equation <- floor(parm[1]) # Parameter can range from 1-9.99
      
      # Read in resistance surface to be optimized
      SHAPE <-  (parm[2])
      Max.SCALE <- (parm[3])
      
      # Apply specified transformation
      rick.eq<-(equation==2||equation==4||equation==6||equation==8)
      if(rick.eq==TRUE & SHAPE>5){
        equation<-9
      }
      
      # Apply specified transformation
      if(equation==1){
        r[[i]] <- Inv.Rev.Monomolecular(rast,parm)
        EQ <- "Inverse-Reverse Monomolecular"
        
      } else if(equation==5){
        r[[i]] <- Rev.Monomolecular(rast,parm)
        EQ <- "Reverse Monomolecular"        
        
      } else if(equation==3){
        r[[i]] <- Monomolecular(rast,parm)   
        EQ <- "Monomolecular"
        
      } else if (equation==7) {
        r[[i]] <- Inv.Monomolecular(rast,parm)
        EQ <- "Inverse Monomolecular"        
        
      } else if (equation==8) {
        r[[i]] <- Inv.Ricker(rast,parm)
        EQ <- "Inverse Ricker"  
        
      } else if (equation==4) {
        r[[i]] <- Ricker(rast,parm)
        EQ <- "Ricker"
        
      } else if (equation==6) {
        r[[i]] <- Rev.Ricker(rast,parm)
        EQ <- "Reverse Ricker"        
        
      } else if (equation==2) {
        r[[i]] <- Inv.Rev.Ricker(rast,parm)
        EQ <- "Inverse-Reverse Ricker"
        
      } else {
        r[[i]] <- (rast*0) #+ 1 #  Cancel layer...set to zero
      } # End if-else  
    } # Close parameter type if-else  
  } # Close layer loop
  
  
  File.name <- File.name
  
  multi_surface <- sum(r) #+ 1 # Add all surfaces together
  
  if(rescale==TRUE) multi_surface <- multi_surface/cellStats(multi_surface,"min") # Rescale to min of 1
  
  if(cellStats(multi_surface,"max")>1e6)  multi_surface<-SCALE(multi_surface,1,1e6) # Rescale surface in case resistance are too high
  #       plot(multi_surface)
  
  if(is.null(out)){
    (multi_surface)
  } else {
  writeRaster(x=multi_surface,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  (multi_surface)
  }
}
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!###

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# TRANSFORM RESISTANCE SURFACES #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#' Apply transformation to continuous resistance surface
#' 
#' Apply on the eight resistance transformations to a continuous resistance surface
#' 
#' @param transformation Transformation equation to apply. Can be provided as the name of the transformation or its numeric equivalent (see details)
#' @param shape Value of the shape parameter
#' @param max Value of the maximum value parameter
#' @param r Resistance surface to be transformed. Can be supplied as full path to .asc file or as a raster object
#' @param out Directory to write transformed .asc file. Default is NULL, and will not export .asc file
#' @usage Resistance.tran(transformation, shape, max, r, out)
#' @return R raster object
#' @details Valid arguements for \code{transformation} are:\cr
#' \tabular{ll}{
#'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#'    \tab 2 = "Inverse-Reverse Ricker"\cr
#'    \tab 3 = "Monomolecular"\cr
#'    \tab 4 = "Ricker"\cr
#'    \tab 5 = "Reverse Monomolecular"\cr
#'    \tab 6 = "Reverse Ricker"\cr
#'    \tab 7 = "Inverse Monomolecular"\cr
#'    \tab 8 = "Inverse Ricker"\cr
#'    \tab 9 = "Distance"\cr
#'    }
#'    
#' The Distance transformation sets all values equal to one. Because of the flexibility of the Ricker function to take a monomolecular shape (try \code{Plot.trans(PARM=c(10,100), Resistance=c(1,10), transformation="Ricker")} to see this), whenever a shape parameter >6 is selected in combination with a Ricker family transformation, the transformation reverts to a Distance transformation. In general, it seems that using a combination of intermediate Ricker and Monomolecular transformations provides the best, most flexible coverasge of parameter space.
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>

Resistance.tran <- function(transformation, shape, max, r, out=NULL){
  if(class(r)[1]!='RasterLayer') {
    R<-raster(r)
    NAME <- basename(r)
    NAME<-sub("^([^.]*).*", "\\1", NAME) 
    names(R)<-NAME
  } else {
    R<-r
    NAME <- r@data@names
  }
  if(is.numeric(transformation)){
    parm<-c(transformation,shape, max)    
  } else {
     parm<-c(get.EQ(transformation),shape, max)
  }
  
  EXPORT.dir<-out
  
  #!##!#
  
      r <-SCALE(data=R,MIN=0,MAX=10)

  # Set equation for continuous surface
      equation <- floor(parm[1]) # Parameter can range from 1-9.99
      
#       # Read in resistance surface to be optimized
#       SHAPE <-  (parm[2])
#       Max.SCALE <- (parm[3])
  
  
  #!##!##!##!##!##!#
      
# Apply specified transformation
    if(equation==1){
      r <- Inv.Rev.Monomolecular(r,parm)
      EQ <- "Inverse-Reverse Monomolecular"
      
    } else if(equation==5){
      r <- Rev.Monomolecular(r,parm)      
      EQ <- "Reverse Monomolecular"        
      
    } else if(equation==3){
      r <- Monomolecular(r,parm)      
      EQ <- "Monomolecular"
      
    } else if (equation==7) {
      r <- Inv.Monomolecular(r,parm)      
      EQ <- "Inverse Monomolecular"  	    
      
    } else if (equation==8) {
      r <- Inv.Ricker(r,parm)
      EQ <- "Inverse Ricker"  
      
    } else if (equation==4) {
      r <- Ricker(r,parm)
      EQ <- "Ricker"
      
    } else if (equation==6) {
      r <- Rev.Ricker(r,parm)
      EQ <- "Reverse Ricker"        
      
    } else if (equation==2) {
      r <- Inv.Rev.Ricker(r,parm)
      EQ <- "Inverse-Reverse Ricker"
   
    } else {
      r <- (r*0)+1 #  Distance
      EQ <- "Distance"    
    } # End if-else  
  File.name <- NAME
  
  if(cellStats(r,"max")>1e6)  r<-SCALE(r,1,1e6) # Rescale surface in case resistance are too high

  if(!is.null(out)){
    writeRaster(x=r,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
    return(r)
  } else {
    return(r)
  }
}
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!###

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# MULTISURFACE OPTIMIZATION FUNCTION WITH GA--NO GAUSSIAN #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#' Optimize multiple resistance surfaces simultaneously
#' 
#' Create composite resistance surface by simultaneously optimizing multiple categoricla and continuous surfaces. This optimization function is designed to be called from GA
#' 
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Should be a vector with parameters specified in the order of resistance surfaces. These values are selected during optimization if called within GA function.
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max')
#' @param quiet Logical, if TRUE, AICc and iteration time will not be printed to the screen at the completion of each iteration. Default = FALSE
#' @return AIC value from mixed effect model
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Resistance.Opt_multi <- function(PARM,CS.inputs=NULL, gdist.inputs=NULL, GA.inputs, Min.Max, quiet=FALSE){
  t1<-proc.time()[3]
 
  EXPORT.dir<-GA.inputs$Write.dir
  #!##!#
#   r <- GA.inputs$Resistance.stack
  File.name="resist_surface"
if(!is.null(CS.inputs)){  
Combine_Surfaces(PARM=PARM,CS.inputs=CS.inputs,GA.inputs=GA.inputs,out=GA.inputs$Write.dir,File.name=File.name,rescale = FALSE)  
  
  CS.resist <- Run_CS2(CS.inputs,GA.inputs,r=multi_surface,EXPORT.dir=GA.inputs$Write.dir,File.name=File.name)
  
  # Run mixed effect model on each Circuitscape effective resistance
  AIC.stat <- suppressWarnings(AIC(MLPE.lmm2(resistance=CS.resist,
                              response=CS.inputs$response,
                              ID=CS.inputs$ID,
                              ZZ=CS.inputs$ZZ,
                              REML=FALSE)))
  ROW <- nrow(CS.inputs$ID)
}

if(!is.null(gdist.inputs)){
  r <- Combine_Surfaces(PARM=PARM,gdist.inputs=gdist.inputs,GA.inputs=GA.inputs,out=NULL,File.name=File.name,rescale = FALSE)
  cd <- Run_gdistance(gdist.inputs,r)
   
  AIC.stat <- suppressWarnings(AIC(MLPE.lmm2(resistance=cd,
                              response=gdist.inputs$response,
                              ID=gdist.inputs$ID,
                              ZZ=gdist.inputs$ZZ,
                              REML=FALSE)))
  ROW <- nrow(gdist.inputs$ID)
}

  k<-max(GA.inputs$parm.index)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(ROW-k-1))  
 
rt<-proc.time()[3]-t1
  if(quiet==FALSE){
   cat(paste0("\t", "Iteration took ", round(rt,digits=2), " seconds to complete"),"\n")
    cat(paste0("\t", "AICc = ",round(AICc,4)),"\n","\n")
  }
  
  
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
}

#!#!#!#!#!#!#!#!#!#! MULTI-SURFACE SUB-LANDSCAPE !#!#!#!#!#!#!
Resistance.Opt_multi_sub <- function(PARM,CS.inputs=NULL, gdist.inputs=NULL, GA.inputs, Min.Max, quiet=FALSE){
  t1<-proc.time()[3]
  
  EXPORT.dir<-GA.inputs$Write.dir
  #!##!#
  #   r <- GA.inputs$Resistance.stack
  File.name="resist_surface"
  if(!is.null(CS.inputs)){  
    Combine_Surfaces(PARM=PARM,CS.inputs=CS.inputs,GA.inputs=GA.inputs,out=GA.inputs$Write.dir,File.name=File.name,rescale = FALSE)  
    
    CS.resist <- Run_CS2(CS.inputs,GA.inputs,r=multi_surface,EXPORT.dir=GA.inputs$Write.dir,File.name=File.name)
    
    # Run mixed effect model on each Circuitscape effective resistance
    AIC.stat <- suppressWarnings(AIC(MLPE.lmm.sub(resistance=CS.resist,
                                               response=CS.inputs$response,
                                               ID=CS.inputs$ID,
                                               ZZ=CS.inputs$ZZ,
                                               sub=CS.inputs$sub,
                                               REML=FALSE)))
    ROW <- nrow(CS.inputs$ID)
  }
  
  if(!is.null(gdist.inputs)){
    r <- Combine_Surfaces(PARM=PARM,gdist.inputs=gdist.inputs,GA.inputs=GA.inputs,out=NULL,File.name=File.name,rescale = FALSE)
    cd <- Run_gdistance(gdist.inputs,r)
    
    AIC.stat <- suppressWarnings(AIC(MLPE.lmm2(resistance=cd,
                                               response=gdist.inputs$response,
                                               ID=gdist.inputs$ID,
                                               ZZ=gdist.inputs$ZZ,
                                               REML=FALSE)))
    ROW <- nrow(gdist.inputs$ID)
  }
  
  k<-max(GA.inputs$parm.index)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(ROW-k-1))  
  
  rt<-proc.time()[3]-t1
  if(quiet==FALSE){
    cat(paste0("\t", "Iteration took ", round(rt,digits=2), " seconds to complete"),"\n")
    cat(paste0("\t", "AICc = ",round(AICc,4)),"\n","\n")
  }
  
  
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
}

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!## 
#!##!##!##!# ITERATIVE SINGLE OPTIMIZATION FUNCTION WITH GA---No Gaussian #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!## 
#' Optimize resistance surfaces individually
#' 
#' Optimize all resistance surfaces that are located in the same directory individually. This optimization function is designed to be called from GA
#' 
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. A vector with parameters specified in the order of resistance surfaces.These values are selected during optimization if called within GA function.
#' @param Resistance Resistance surface to be optimized. This should be an R raster object. If not specified, the function will attempt to find the a resistance surface from \code{GA.inputs}
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max'). Default in 'max'
#' @param iter A counter for the number of surfaces that will be optimized
#' @param quiet Logical, if TRUE AICc and iteration duration will not be printed to the screen at the completion of each iteration.
#' @return AIC value from mixed effect model
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Resistance.Opt_single <- function(PARM,Resistance,CS.inputs=NULL, gdist.inputs=NULL, GA.inputs, Min.Max='max',iter, quiet=FALSE){
  t1<-proc.time()[3]

  
  EXPORT.dir<-GA.inputs$Write.dir
  #!##!#
  r <- Resistance
   
  if(GA.inputs$surface.type[iter]=="cat"){
    PARM<-PARM/min(PARM)
    parm <-PARM
    df <- data.frame(id=unique.rast(r),PARM) # Data frame with original raster values and replacement values
    r <-subs(r,df)    
       
  } else {
    r<-SCALE(r,0,10)
    
    # Set equation for continuous surface
    equation <- floor(PARM[1]) # Parameter can range from 1-9.99
    
    # Read in resistance surface to be optimized
    SHAPE <- (PARM[2])
    Max.SCALE <- (PARM[3])
    
    # Apply specified transformation
    rick.eq<-(equation==2||equation==4||equation==6||equation==8)
    if(rick.eq==TRUE & SHAPE>5){
      equation<-9
    }
    
    if(equation==1){
      r <- Inv.Rev.Monomolecular(r,parm=PARM)
      EQ <- "Inverse-Reverse Monomolecular"
      
    } else if(equation==5){
      r <- Rev.Monomolecular(r,parm=PARM)      
      EQ <- "Reverse Monomolecular"        
      
    } else if(equation==3){
      r <- Monomolecular(r,parm=PARM)      
      EQ <- "Monomolecular"
      
    } else if (equation==7) {
      r <- Inv.Monomolecular(r,parm=PARM)      
      EQ <- "Inverse Monomolecular"        
      
    } else if (equation==8) {
      r <- Inv.Ricker(r,parm=PARM)
      EQ <- "Inverse Ricker"  
      
    } else if (equation==4) {
      r <- Ricker(r,parm=PARM)
      EQ <- "Ricker"
      
    } else if (equation==6) {
      r <- Rev.Ricker(r,parm=PARM)
      EQ <- "Reverse Ricker"        
      
    } else if (equation==2) {
      r <- Inv.Rev.Ricker(r,parm=PARM)
      EQ <- "Inverse-Reverse Ricker"
      
    } else {
      r <- (r*0)+1 #  Distance
      EQ <- "Distance"    
    } # End if-else     
  } # Close parameter type if-else 
  
  File.name <- "resist_surface"
  if(cellStats(r,"max")>1e6)  r<-SCALE(r,1,1e6) # Rescale surface in case resistance are too high
  r <- reclassify(r, c(-Inf,1e-06, 1e-06,1e6,Inf,1e6))
  

  if(!is.null(CS.inputs)){
    writeRaster(x=r,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
    CS.resist <- Run_CS2(CS.inputs,GA.inputs,r=r,EXPORT.dir=GA.inputs$Write.dir,File.name=File.name)
    
    # Run mixed effect model on each Circuitscape effective resistance
    AIC.stat <- suppressWarnings(AIC(MLPE.lmm2(resistance=CS.resist,
                                response=CS.inputs$response,
                                ID=CS.inputs$ID,
                                ZZ=CS.inputs$ZZ,
                                REML=FALSE)))
    ROW <- nrow(CS.inputs$ID)
    
  }  
  
  if(!is.null(gdist.inputs)){
    cd <- Run_gdistance(gdist.inputs,r)
    
    AIC.stat <- suppressWarnings(AIC(MLPE.lmm2(resistance=cd,
                                response=gdist.inputs$response,
                                ID=gdist.inputs$ID,
                                ZZ=gdist.inputs$ZZ,
                                REML=FALSE)))
    ROW <- nrow(gdist.inputs$ID)
  } 

  k<-length(PARM)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(ROW-k-1))
  
  rt<-proc.time()[3]-t1  
  if(quiet==FALSE){
    cat(paste0("\t", "Iteration took ", round(rt,digits=2), " seconds to complete"),"\n")
#     cat(paste0("\t", EQ,"; ",round(SHAPE,digits=2),"; ", round(Max.SCALE,digits=2)),"\n")
    cat(paste0("\t", "AICc = ",round(AICc,4)),"\n")
    cat(paste0("\t", EQ, " | Shape = ",PARM[2]," | Max = ",PARM[3]),"\n","\n")

    
  }
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
}


#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!## 
#!##!##!##!# SUBPLOT Single Surface OPTIMIZATION FUNCTION #!##!##!##!# 
#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#

Resistance.Opt_single_sub <- function(PARM,Resistance,CS.inputs=NULL, gdist.inputs=NULL, GA.inputs, Min.Max='max',iter, quiet=FALSE){
  t1<-proc.time()[3]
  
  
  EXPORT.dir<-GA.inputs$Write.dir
  #!##!#
  r <- Resistance
  
  if(GA.inputs$surface.type[iter]=="cat"){
    PARM<-PARM/min(PARM)
    parm <-PARM
    df <- data.frame(id=unique.rast(r),PARM) # Data frame with original raster values and replacement values
    r <-subs(r,df)    
    
  } else {
    r<-SCALE(r,0,10)
    
    # Set equation for continuous surface
    equation <- floor(PARM[1]) # Parameter can range from 1-9.99
    
    # Read in resistance surface to be optimized
    SHAPE <- (PARM[2])
    Max.SCALE <- (PARM[3])
    
    # Apply specified transformation
    rick.eq<-(equation==2||equation==4||equation==6||equation==8)
    if(rick.eq==TRUE & SHAPE>5){
      equation<-9
    }
    
    if(equation==1){
      r <- Inv.Rev.Monomolecular(r,parm=PARM)
      EQ <- "Inverse-Reverse Monomolecular"
      
    } else if(equation==5){
      r <- Rev.Monomolecular(r,parm=PARM)      
      EQ <- "Reverse Monomolecular"        
      
    } else if(equation==3){
      r <- Monomolecular(r,parm=PARM)      
      EQ <- "Monomolecular"
      
    } else if (equation==7) {
      r <- Inv.Monomolecular(r,parm=PARM)      
      EQ <- "Inverse Monomolecular"        
      
    } else if (equation==8) {
      r <- Inv.Ricker(r,parm=PARM)
      EQ <- "Inverse Ricker"  
      
    } else if (equation==4) {
      r <- Ricker(r,parm=PARM)
      EQ <- "Ricker"
      
    } else if (equation==6) {
      r <- Rev.Ricker(r,parm=PARM)
      EQ <- "Reverse Ricker"        
      
    } else if (equation==2) {
      r <- Inv.Rev.Ricker(r,parm=PARM)
      EQ <- "Inverse-Reverse Ricker"
      
    } else {
      r <- (r*0)+1 #  Distance
      EQ <- "Distance"    
    } # End if-else     
  } # Close parameter type if-else 
  
  File.name <- "resist_surface"
  if(cellStats(r,"max")>1e6)  r<-SCALE(r,1,1e6) # Rescale surface in case resistance are too high
  r <- reclassify(r, c(-Inf,1e-06, 1e-06,1e6,Inf,1e6))
  
  
  if(!is.null(CS.inputs)){
    writeRaster(x=r,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
    CS.resist <- Run_CS2(CS.inputs,GA.inputs,r=r,EXPORT.dir=GA.inputs$Write.dir,File.name=File.name)
    
    # Run mixed effect model on each Circuitscape effective resistance
    AIC.stat <- suppressWarnings(AIC(MLPE.lmm.sub(resistance=CS.resist,
                                                  response=CS.inputs$response,
                                                  ID=CS.inputs$ID,
                                                  ZZ=CS.inputs$ZZ,
                                                  REML=FALSE,
                                                  sub=CS.inputs$sub)
                                     ))
    ROW <- nrow(CS.inputs$ID)
    
  }  
  
  if(!is.null(gdist.inputs)){
    cd <- Run_gdistance(gdist.inputs,r)
    
    AIC.stat <- suppressWarnings(AIC(MLPE.lmm.sub(resistance=cd,
                                               response=gdist.inputs$response,
                                               ID=gdist.inputs$ID,
                                               ZZ=gdist.inputs$ZZ,
                                               REML=FALSE,
                                               sub=CS.inputs$sub)))
    ROW <- nrow(gdist.inputs$ID)
  } 
  
  k<-length(PARM)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(ROW-k-1))
  
  rt<-proc.time()[3]-t1  
  if(quiet==FALSE){
    cat(paste0("\t", "Iteration took ", round(rt,digits=2), " seconds to complete"),"\n")
    #     cat(paste0("\t", EQ,"; ",round(SHAPE,digits=2),"; ", round(Max.SCALE,digits=2)),"\n")
    cat(paste0("\t", "AICc = ",round(AICc,4)),"\n")
    cat(paste0("\t", EQ, " | Shape = ",PARM[2]," | Max = ",PARM[3]),"\n","\n")
    
    
  }
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
}

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!# 
#!##!##!##!# PLOT response CURVES #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!# 
#' Plot continuous surface transformation
#' 
#' Plots a transformed continuous resistance surface against the original resistance values
#' 
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. A vector of two parameters is required. The first term is the value of shape parameter (c), and the second term is the value of maximum scale parameter (b)
#' @param Resistance Accepts three types of inputs. Provide either the path to the raw, untransformed resistance surface file or specify an R raster object. Alternatively, supply a vector with the minimum and manximum values (e.g., c(1,10))
#' @param transformation Transformation equation to apply. Can be provided as the name of the transformation or its numeric equivalent (see details)
#' @param print.dir Specify the directory where a .tiff of the transformation will be written (Default = NULL)
#' @param Name Name of resistance surface being transformed (optional). This will be added to the output file name.
#' @return plot of transformed resistance values against original resistance values
#' @details This function will create a ggplot object and plot, so it requires \pkg{ggplot2} to be installed.\cr 
#' Equation names can be:
#' \tabular{ll}{
#'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#'    \tab 2 = "Inverse-Reverse Ricker"\cr
#'    \tab 3 = "Monomolecular"\cr
#'    \tab 4 = "Ricker"\cr
#'    \tab 5 = "Reverse Monomolecular"\cr
#'    \tab 6 = 'Reverse Ricker"\cr
#'    \tab 7 = "Inverse Monomolecular"\cr
#'    \tab 8 = "Inverse Ricker"\cr
#'    \tab 9 = "Distance"\cr
#'    }
#'    
#' The "Distance" equation sets all cell values equal to 1.
#' Because of the flexibility of the Ricker function to take a monomolecular shape (try \code{Plot.trans(PARM=c(10,100), Resistance=c(1,10), transformation="Ricker")} to see this), whenever a shape parameter >6 is selected in combination with a Ricker family transformation, the transformation reverts to a Distance transformation. In general, it seems that using a combination of intermediate Ricker and Monomolecular transformations provides the best, most flexible coverasge of parameter space. This constraint has not been implemented in the \code{Plot.tans} function.
#' @usage Plot.trans(PARM, Resistance, transformation, print.dir, Name)
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>

Plot.trans <- function(PARM,Resistance,transformation, print.dir=NULL, Name="layer"){
    if(length(Resistance)==2) {
        r <- Resistance
        Mn=min(r)
        Mx=max(r) 
        NAME <- Name
    } else if(class(Resistance)[1]!='RasterLayer') {
        r<-raster(Resistance)          
        NAME <- basename(Resistance)
        NAME<-sub("^([^.]*).*", "\\1", NAME) 
        names(r)<-NAME
        Mn=cellStats(r,stat='min')
        Mx=cellStats(r,stat='max')   
          
  } else {
        r <- Resistance
        Mn=cellStats(r,stat='min')
        Mx=cellStats(r,stat='max')
        NAME <- basename(r@file@name)
        NAME<-sub("^([^.]*).*", "\\1", NAME) 
        names(r)<-NAME
  }  
  
  if(Name!="layer"){
    NAME<-Name
  }
  
  # Make vector of data 
  original <- seq(from=Mn,to=Mx,length.out=200)
  dat.t <- SCALE.vector(data=original,0,10)
  
  SHAPE <- PARM[[1]]
  Max.SCALE <- PARM[[2]]
  
  if(is.numeric(transformation)){
   equation<-get.EQ(transformation)
    
  } else {
  equation<-transformation
  }
  # Set equation/name combination
  if(equation=="Distance") {
    Trans.vec <- (dat.t*0)+1 
    TITLE <- "Distance"
    
  } else if(equation=="Inverse-Reverse Monomolecular"){
    SIGN<- -1 # Inverse
    Trans.vec <-  SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN # Monomolecular
    Trans.vec <- rev(SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec))))
    TITLE <- "Inverse-Reverse Monomolecular"
    
  } else if(equation=="Inverse Monomolecular"){
    SIGN<- -1 # Inverse
    Trans.vec <-  SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN # Monomolecular
    Trans.vec <- (SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec))))
    TITLE <- "Inverse Monomolecular"
    
  } else if(equation=="Monomolecular") {
    SIGN<- 1
    Trans.vec <- SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN # Monomolecular
    TITLE <- "Monomolecular"
    
  } else if(equation=="Reverse Monomolecular") {
    SIGN<- 1
    Trans.vec <-  SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN # Monomolecular
    Trans.vec <- rev(Trans.vec) # Reverse
    TITLE <- "Reverse Monomolecular"  
    
  } else if (equation=="Inverse Ricker") {
    SIGN <- -1 # Inverse
    Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
    Trans.vec <- SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec)))
    TITLE <- "Inverse Ricker" 
    
  } else if (equation=="Reverse Ricker") { 
    SIGN <- 1
    Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
    Trans.vec <- rev(Trans.vec) # Reverse
    TITLE <- "Reverse Ricker" 
  
  } else if (equation=="Inverse-Reverse Ricker") { 
    SIGN <- -1 # Inverse
    Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
    Trans.vec <- rev(SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec)))) # Reverse
    TITLE <- "Inverse-Reverse Ricker" 
    
  }  else  {
    SIGN <- 1
    Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
    TITLE <- "Ricker"
  } 
  transformed<-Trans.vec
  plot.data<-data.frame(original,transformed)
  x.break<-pretty(original*1.07)
  y.break<-pretty(transformed*1.07)
  
 ( p<- ggplot(plot.data,aes(x=original,y=transformed)) +
    ggtitle(equation) +
    theme_bw() +
    geom_line(size=1.5) +
    xlab(expression(bold("Original data values"))) +
    ylab(expression(bold("Transformed data values"))) +
    theme(plot.title = element_text(lineheight=2, face="bold",size=20),
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),  
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
     scale_x_continuous(limits=c(min(original),max(original)),breaks=x.break) +
     scale_y_continuous(limits=c(min(transformed),max(transformed)),breaks=y.break) 
 )
  
  if(!is.null(print.dir)){
    tiff(filename=paste0(print.dir,equation,"_Transformation_",NAME ,".tif"),width=160,height=150,units="mm",res=300,compression="lzw")
    print(p)
    dev.off()    
    return(p)
  }
  print(p)
  return(p)
  
}


#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
#!##!##!##!# OPTIMIZATION FUNCTION USING GA STARTS: CONTINUOUS #!##!##!##!# 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##  
Resistance.Optimization_cont.nlm<-function(PARM,Resistance,equation, get.best,CS.inputs=NULL, gdist.inputs=NULL,Min.Max,quiet=FALSE, write.dir) {
  CS.program<-CS.inputs$CS.program
  
      EXPORT.dir<-write.dir  
  
  r <- Resistance
  r<-SCALE(r,0,10)
  Name <- Resistance@data@names
  t1<-Sys.time()
  
  equation<-floor(equation)
  EQ <-get.EQ(equation)
  
  #   if(equation<7){
  SHAPE<-ifelse(exp(PARM[1])>10000,10000,exp(PARM[1])) # Upper boundaries on parameters
  Max.SCALE<-ifelse(exp(PARM[2])>1e6,1e6,exp(PARM[2]))
  cat("\n", Name, as.character(EQ),paste0("| Shape = ",SHAPE,"; "),paste0("Maximum scale = ",Max.SCALE,"\n"))  
   
  # Apply specified transformation
  rick.eq<-(equation==2||equation==4||equation==6||equation==8)
  if(rick.eq==TRUE & SHAPE>6){
    equation<-9
  }
  
  if(equation==1){
    SIGN=-1 # Inverse
    R <- SIGN*Max.SCALE*(1-exp(-1*r/SHAPE))+SIGN # Monomolecular
    R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
    R.vec <- rev(R) # Reverse
    rast.R <- setValues(R,values=R.vec)
    r <- rast.R
    EQ <- "Inverse-Reverse Monomolecular"
    
  } else if(equation==5){
    SIGN=1
    R <- SIGN*Max.SCALE*(1-exp(-1*r/SHAPE))+SIGN # Monomolecular
    R.vec <- rev(R) # Reverse
    rast.R <- setValues(R,values=R.vec)
    r <- rast.R
    EQ <- "Reverse Monomolecular"        
    
  } else if(equation==3){
    SIGN=1
    r <- SIGN*Max.SCALE*(1-exp(-1*r/SHAPE))+SIGN # Monomolecular    
    EQ <- "Monomolecular"
    
  } else if (equation==7) {
    SIGN=-1 #Inverse
    R <- SIGN*Max.SCALE*(1-exp(-1*r/SHAPE))+SIGN # Monomolecular
    r <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
    EQ <- "Inverse Monomolecular"        
    
  } else if (equation==8) {
    SIGN=-1 #Inverse
    R <- SIGN*(Max.SCALE*r*exp(-1*r/SHAPE))+SIGN # Ricker
    r <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
    EQ <- "Inverse Ricker"  
    
  } else if (equation==4) {
    SIGN=1
    r <- SIGN*(Max.SCALE*r*exp(-1*r/SHAPE))+SIGN #  Ricker
    EQ <- "Ricker"
    
  } else if (equation==6) {
    SIGN=1
    R <- SIGN*(Max.SCALE*r*exp(-1*r/SHAPE))+SIGN #  Ricker
    R.vec <- rev(R)
    rast.R <- setValues(R,values=R.vec)
    r <- rast.R
    EQ <- "Reverse Ricker"        
    
  } else if (equation==2) {
    SIGN=-1 # Inverse
    R <- SIGN*(Max.SCALE*r*exp(-1*r/SHAPE))+SIGN # Ricker
    R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
    R.vec <- rev(R) # Reverse
    rast.R <- setValues(R,values=R.vec)
    r <- rast.R
    EQ <- "Inverse-Reverse Ricker"
    
  } else {
    r <- (r*0)+1 #  Distance
    EQ <- "Distance"    
  } # End if-else  
  
  #   File.name <- paste0("exp",TRAN,"_MAX", MAX)
  File.name <- if (get.best==FALSE){
    File.name <- "optim_iter"
  } else {
    File.name <- Name
  }
  
  if(cellStats(r,"max")>1e6)  r<-SCALE(r,1,1e6) # Rescale surface in case resistance are too high
  r <- reclassify(r, c(-Inf,1e-06, 1e-06,1e6,Inf,1e6))
  
  if(!is.null(CS.inputs)){
    writeRaster(x=r,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
    CS.resist <- Run_CS2(CS.inputs,GA.inputs,r=r,EXPORT.dir=GA.inputs$Write.dir,File.name=File.name)
    
    # Run mixed effect model on each Circuitscape effective resistance
    AIC.stat <- suppressWarnings(AIC(MLPE.lmm2(resistance=CS.resist,
                                response=CS.inputs$response,
                                ID=CS.inputs$ID,
                                ZZ=CS.inputs$ZZ,
                                REML=FALSE)))
    ROW <- nrow(CS.inputs$ID)
    
  }  
  
  if(!is.null(gdist.inputs)){
    cd <- Run_gdistance(gdist.inputs,r)
    
    AIC.stat <- suppressWarnings(AIC(MLPE.lmm2(resistance=cd,
                                response=gdist.inputs$response,
                                ID=gdist.inputs$ID,
                                ZZ=gdist.inputs$ZZ,
                                REML=FALSE)))
    ROW <- nrow(gdist.inputs$ID)
  }   
 
  k<-length(PARM)+2
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(ROW-k-1))
    
  t2 <-Sys.time()
  if(quiet==FALSE){    
  cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
  cat(paste0("\t", "AICc = ",round(AICc,3)),"\n")
  cat(paste0("\t", EQ, " | Shape = ",PARM[2]," | Max = ",PARM[3]),"\n","\n")
  
  }
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized    
  
}
#!##!##!##!##!##!##!##!##!##!##!##!##!##!###

# Run Mixed effects models, recovery parameter estimates
# ' Obtain coefficients from maximum likelihood population effects mixed effects model (MLPE)
# ' 
# ' Runs MLPE as detailed by Clarke et al. (2002). This function is designed to generate a table of the fitted mixed effects model coefficients
# ' 
# ' @param resistance Directory containing optimized .asc resistance files (optimized resistance surfaces from CIRCUITSCAPE), or costDistance object from optimizing with gdistance
# ' @param genetic.dist Lower half of pairwise genetic distance matrix
# ' @param out.dir If specified, a .csv table will printed to the specified directory (Default = NULL)
# ' @return A table of MLPE fitted model coefficients
# ' @author Bill Peterman <Bill.Peterman@@gmail.com>
# ' @usage MLPE.lmm_coef(resistance, genetic.dist, out.dir)
# ' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

MLPE.lmm_coef <- function(resistance, genetic.dist,out.dir=NULL, method, ID=NULL, ZZ=NULL, sub=NULL){ 
  if(method=="cs"){
    if(is.null(sub)){
      response=genetic.dist
      resist.mat<-list.files(resistance,pattern="*_resistances.out",full.names=TRUE)
      resist.names<-gsub(pattern="_resistances.out","",x=list.files(resistance,pattern="*_resistances.out"))
      COEF.Table<-array()
         for(i in 1:length(resist.mat)){
               m<-length(read.table(resist.mat[i])[-1,-1])
               mm<-read.table(resist.mat[i])[-1,-1]
               mm <- lower(mm)
               mm <- mm[which(mm!=-1)]
                 if(is.null(ID)){
                     ID<-To.From.ID(POPS=m)                   
                       }
                 if(is.null(ZZ)){
                     ZZ<-ZZ.mat(ID=ID)                         
                       }
                       
               resistance<-scale(mm,center=TRUE,scale=TRUE)
               dat<-data.frame(ID,resistance=resistance,response=response)
               colnames(dat)<-c("pop1","pop2","resistance","response")
                       
               # Assign value to layer
               LAYER<-assign(resist.names[i],value=dat$cs.matrix)
                       
               # Fit model
               mod <- lFormula(response ~ resistance + (1|pop1), data=dat,REML=TRUE)
               mod$reTrms$Zt <- ZZ
               dfun <- do.call(mkLmerDevfun,mod)
               opt <- optimizeLmer(dfun)
               Mod.Summary <- summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
               COEF<-Mod.Summary$coefficients
               row.names(COEF)<-c("Intercept",resist.names[i])
               COEF.Table<-rbind(COEF.Table, COEF)
                     }
    } else { # Get sub-landscape coefficients
      response=genetic.dist
      resist.mat<-list.files(resistance,pattern="*_resistances.out",full.names=TRUE)
      resist.names<-gsub(pattern="_resistances.out","",x=list.files(resistance,pattern="*_resistances.out"))
      COEF.Table<-array()
      for(i in 1:length(resist.mat)){
        m<-length(read.table(resist.mat[i])[-1,-1])
        mm<-read.table(resist.mat[i])[-1,-1]
        mm <- lower(mm)
        mm <- mm[which(mm!=-1)]
        if(is.null(ID)){
          ID<-To.From.ID(POPS=m)
          
        }
        
        if(is.null(ZZ)){
          ZZ<-ZZ.mat(ID=ID)
          
        }
        
        resistance<-scale(mm,center=TRUE,scale=TRUE)
        dat<-data.frame(ID,resistance=resistance,response=response, sub=sub)
        colnames(dat)<-c("pop1","pop2","resistance","response","sub")
        
        # Assign value to layer
        LAYER<-assign(resist.names[i],value=dat$cs.matrix)
        
        # Fit model
        mod <- MLPE.lmm.sub(resistance = resistance,
                            response = response,
                            REML = TRUE,
                            ID = ID,
                            ZZ = ZZ,
                            sub = sub)
        
        Mod.Summary <- summary(mod)
        COEF<-Mod.Summary$coefficients
        row.names(COEF)<-c("Intercept",resist.names[i])
        COEF.Table<-rbind(COEF.Table, COEF)
      }   
    }
  
  } else {   # Use Gdist    
    response=genetic.dist
    resist.mat<-list.files(resistance,pattern="*.rda",full.names=TRUE)
    resist.names<-gsub(pattern=".rda","",x=list.files(resistance,pattern=".rda"))
    COEF.Table<-array()
  for(i in 1:length(resist.mat)){
    load(resist.mat[i])
    mm <-lower(as.matrix(cd))
    m<-attr(cd,"Size")
    ID<-To.From.ID(POPS=m)
    ZZ<-ZZ.mat(ID=ID)
    cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
    cs.unscale<-mm
    dat<-cbind(ID,cs.matrix,response)  
    
    # Assign value to layer
    LAYER<-assign(resist.names[i],value=dat$cs.matrix)
    
    # Fit model
    mod <- lFormula(response ~ LAYER + (1|pop1), data=dat,REML=TRUE)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    Mod.Summary <- summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
    COEF<-Mod.Summary$coefficients
    row.names(COEF)<-c("Intercept",resist.names[i])
    COEF.Table<-rbind(COEF.Table, COEF)
  }
}
    
  
  if(is.null(out.dir)){
    COEF.Table<-(COEF.Table[-1,])
  } else {
    COEF.Table<-COEF.Table[-1,]
    write.table(COEF.Table,file=paste0(out.dir,"MLPE_coeff_Table.csv"),sep=",",row.names=T,col.names=NA)
    return(COEF.Table)
  }
}

#!#!#!#!#!#!#! SUB-LANDSCAPE MLPE COEFFICIENTS

# MLPE.lmm_coef.sub <- function(resistance, genetic.dist,out.dir=NULL, method, ID=NULL, ZZ=NULL, sub=NULL){ 
#   if(method=="cs"){
#     if(is.null(sub)){
#       response=genetic.dist
#       resist.mat<-list.files(resistance,pattern="*_resistances.out",full.names=TRUE)
#       resist.names<-gsub(pattern="_resistances.out","",x=list.files(resistance,pattern="*_resistances.out"))
#       COEF.Table<-array()
#       for(i in 1:length(resist.mat)){
#         m<-length(read.table(resist.mat[i])[-1,-1])
#         mm<-read.table(resist.mat[i])[-1,-1]
#         mm <- lower(mm)
#         mm <- mm[which(mm!=-1)]
#         if(is.null(ID)){
#           ID<-To.From.ID(POPS=m)                   
#         }
#         if(is.null(ZZ)){
#           ZZ<-ZZ.mat(ID=ID)                         
#         }
#         
#         resistance<-scale(mm,center=TRUE,scale=TRUE)
#         dat<-data.frame(ID,resistance=resistance,response=response)
#         colnames(dat)<-c("pop1","pop2","resistance","response")
#         
#         # Assign value to layer
#         LAYER<-assign(resist.names[i],value=dat$cs.matrix)
#         
#         # Fit model
#         mod <- lFormula(response ~ resistance + (1|pop1), data=dat,REML=TRUE)
#         mod$reTrms$Zt <- ZZ
#         dfun <- do.call(mkLmerDevfun,mod)
#         opt <- optimizeLmer(dfun)
#         Mod.Summary <- summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
#         COEF<-Mod.Summary$coefficients
#         row.names(COEF)<-c("Intercept",resist.names[i])
#         COEF.Table<-rbind(COEF.Table, COEF)
#       }
#     } else { # Get sub-landscape coefficients
#       response=genetic.dist
#       resist.mat<-list.files(resistance,pattern="*_resistances.out",full.names=TRUE)
#       resist.names<-gsub(pattern="_resistances.out","",x=list.files(resistance,pattern="*_resistances.out"))
#       COEF.Table<-array()
#       for(i in 1:length(resist.mat)){
#         m<-length(read.table(resist.mat[i])[-1,-1])
#         mm<-read.table(resist.mat[i])[-1,-1]
#         mm <- lower(mm)
#         mm <- mm[which(mm!=-1)]
#         if(is.null(ID)){
#           ID<-To.From.ID(POPS=m)
#           
#         }
#         
#         if(is.null(ZZ)){
#           ZZ<-ZZ.mat(ID=ID)
#           
#         }
#         
#         resistance<-scale(mm,center=TRUE,scale=TRUE)
#         dat<-data.frame(ID,resistance=resistance,response=response, sub=sub)
#         colnames(dat)<-c("pop1","pop2","resistance","response","sub")
#         
#         # Assign value to layer
#         LAYER<-assign(resist.names[i],value=dat$cs.matrix)
#         
#         # Fit model
#         mod <- MLPE.lmm.sub(resistance = resistance,
#                             response = response,
#                             REML = TRUE,
#                             ID = ID,
#                             ZZ = ZZ,
#                             sub = sub)
#         
#         Mod.Summary <- summary(mod)
#         COEF<-Mod.Summary$coefficients
#         row.names(COEF)<-c("Intercept",resist.names[i])
#         COEF.Table<-rbind(COEF.Table, COEF)
#       }   
#     }
#     
#   } else {   # Use Gdist    
#     response=genetic.dist
#     resist.mat<-list.files(resistance,pattern="*.rda",full.names=TRUE)
#     resist.names<-gsub(pattern=".rda","",x=list.files(resistance,pattern=".rda"))
#     COEF.Table<-array()
#     for(i in 1:length(resist.mat)){
#       load(resist.mat[i])
#       mm <-lower(as.matrix(cd))
#       m<-attr(cd,"Size")
#       ID<-To.From.ID(POPS=m)
#       ZZ<-ZZ.mat(ID=ID)
#       cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
#       cs.unscale<-mm
#       dat<-cbind(ID,cs.matrix,response)  
#       
#       # Assign value to layer
#       LAYER<-assign(resist.names[i],value=dat$cs.matrix)
#       
#       # Fit model
#       mod <- lFormula(response ~ LAYER + (1|pop1), data=dat,REML=TRUE)
#       mod$reTrms$Zt <- ZZ
#       dfun <- do.call(mkLmerDevfun,mod)
#       opt <- optimizeLmer(dfun)
#       Mod.Summary <- summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
#       COEF<-Mod.Summary$coefficients
#       row.names(COEF)<-c("Intercept",resist.names[i])
#       COEF.Table<-rbind(COEF.Table, COEF)
#     }
#   }
#   
#   
#   if(is.null(out.dir)){
#     COEF.Table<-(COEF.Table[-1,])
#   } else {
#     COEF.Table<-COEF.Table[-1,]
#     write.table(COEF.Table,file=paste0(out.dir,"MLPE_coeff_Table.csv"),sep=",",row.names=T,col.names=NA)
#     return(COEF.Table)
#   }
# }



# Run Mixed effects models, recovery parameter estimates
#' Run maximum likelihood population effects mixed effects model (MLPE)
#' 
#' Runs MLPE as detailed by Clarke et al. (2002). This function will run the model and return glmer object
#' 
#' @param resistance Path to pairwise resistance distance matrix (resistances.out) from CS results. Alternatively, provide the pairwise resistances created from optimizing with `gdistance` (result of Run_gdistance).
#' @param pairwise.genetic Lower half of pairwise genetic distance matrix
#' @param REML Logical. If TRUE, mixed effects model will be fit using restricted maximum likelihood. Default = FALSE
#' @param ID The to_from ID list for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @param ZZ The sparse matrix object for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @return A lmer object from the fitted model
#' @details An AIC value will only be returned if \code{REML = FALSE}

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage MLPE.lmm(resistance, pairwise.genetic, REML, ID, ZZ)
#' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

MLPE.lmm <- function(resistance, pairwise.genetic, REML=FALSE, ID=NULL, ZZ=NULL){ 
  response=pairwise.genetic

  if(class(resistance)[[1]]=='dist'){
    mm<-lower(as.matrix(resistance))
    m<-attr(resistance,"Size")
    mm <- mm[which(mm!=-1)]
    
    if(is.null(ID)){
    ID<-To.From.ID(POPS=m)
    }
    if(is.null(ZZ)){
    ZZ<-ZZ.mat(ID=ID)
    }
    cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
    
  } else {
    mm<-(read.table(resistance)[-1,-1])
    m<-nrow(mm)
    mm <- lower(mm)
    mm <- mm[which(mm!=-1)]
    
    if(is.null(ID)){
      ID<-To.From.ID(POPS=m)
    }
    if(is.null(ZZ)){
      ZZ<-ZZ.mat(ID=ID)
    }
    cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
  }
  
    dat<-data.frame(ID,resistance=cs.matrix,response=response)
    colnames(dat)<-c("pop1","pop2","resistance","response")
  
    # Assign value to layer
#     LAYER<-assign("Resist",value=dat$cs.matrix)
    
    # Fit model
    mod <- lFormula(response ~ resistance + (1|pop1), data=dat,REML=REML)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    MOD <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))   
  return(MOD)
}


MLPE.lmm2 <- function(resistance, response, REML=FALSE, ID, ZZ){ 
  if(class(resistance)[1]!='dist'){
    resistance <- resistance[which(resistance!=-1)]    
    dat<-data.frame(ID,resistance=resistance,response=response)
    colnames(dat)<-c("pop1","pop2","resistance","response")
    } else {
    resistance <-lower(as.matrix(resistance))
    resistance <- resistance[which(resistance!=-1)]
    dat<-data.frame(ID,resistance=resistance,response=response)
    colnames(dat)<-c("pop1","pop2","resistance","response")
    
  }  
 
  # Assign value to layer
#   LAYER<-assign("Resist",value=dat$resistance)
  
  # Fit model
  mod <- lFormula(response ~ resistance + (1|pop1), data=dat,REML=REML)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  MOD <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))   
  return(MOD)
}

#!#!#!#!#!#!#!#!#!#!#!

MLPE.lmm.sub <- function(resistance, response, REML=FALSE, ID, ZZ, sub){ 
  if(class(resistance)[1]!='dist'){
    resistance <- resistance[which(resistance!=-1)]    
    dat<-data.frame(ID,resistance=resistance,response=response, sub=sub)
    colnames(dat)<-c("pop1","pop2","resistance","response", "sub")
  } else {
    resistance <-lower(as.matrix(resistance))
    resistance <- resistance[which(resistance!=-1)]
    dat<-data.frame(ID,resistance=resistance,response=response, sub=sub)
    colnames(dat)<-c("pop1","pop2","resistance","response","sub")    
  }
  # Fit model
  mod <- lFormula(response ~ resistance + (1|sub/pop1), data=dat,REML=REML)
#   subplots <- length(unique(sub[,1]))
#   sub_dim <- mod$reTrms$Zt@Dim[1]
#   subplot_levels <- mod$reTrms$Zt[c((sub_dim-subplots+1):sub_dim),]
#   ZZ.sub <- rBind(ZZ, subplot_levels, deparse.level = 1)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  MOD <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))   
  return(MOD)
}
#!##!##!##!##!##!##!##!##!##!##!##
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
#' @param sublandscape Default is NULL. If using the "pairs_to_include" advanced option in CIRCUITSCAPE, then this  must be included. Should be in the form of a single-column data frame of the same length as the number of pairwise comparisons, with values indicating the sublandscape within which pairs reside. 
#' @return A multipanel panel .tif including histogram of residuals and qqplot of fitted mixed effects model 

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage Diagnostic.Plots(resistance.mat, genetic.dist, XLAB,YLAB, plot.dir, type, name, ID, ZZ, sublandscape)

Diagnostic.Plots<-function(resistance.mat, genetic.dist, XLAB="Estimated resistance",YLAB ="Genetic distance",plot.dir, type="categorical", name=NULL, ID=NULL, ZZ=NULL, sublandscape=NULL){
  if(is.null(sublandscape)){
  if(length(resistance.mat)==1){
    response=genetic.dist
    if(is.null(name)){
    NAME<-gsub(pattern="*_resistances.out","",x=(basename(resistance.mat)))
    }
    mm<-read.table(resistance.mat)[-1,-1]
    m<-length(mm)
    mm <-lower(mm)
    mm <- mm[which(mm!=-1)]
    
    if(is.null(ID)){
      ID<-To.From.ID(POPS=m)
    }
    if(is.null(ZZ)){
      ZZ<-ZZ.mat(ID=ID)
    }
    
    cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
    cs.unscale<-mm
    dat<-data.frame(ID,cs.matrix=cs.matrix,response=response)
    colnames(dat) <- c("pop1","pop2","cs.matrix","response")
    
    # Assign value to layer
    LAYER<-assign("LAYER",value=dat$cs.matrix)
  
    # Fit model
    mod <- lFormula(response ~ LAYER + (1|pop1), data=dat,REML=TRUE)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    Mod <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  }
  
  
  
  if(length(resistance.mat)>1){
    response=genetic.dist
    if(is.null(name)){
      stop("Output file 'name' must be specified!!!")
    }
    NAME<-name
    mm<-lower(as.matrix(resistance.mat))
    m<-attr(resistance.mat,"Size")
    mm <- mm[which(mm!=-1)]
    
    if(is.null(ID)){
      ID<-To.From.ID(POPS=m)
    }
    if(is.null(ZZ)){
      ZZ<-ZZ.mat(ID=ID)
    }
    cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
    cs.unscale<-mm
    dat<-data.frame(ID,cs.matrix=cs.matrix,response=response)
    colnames(dat) <- c("pop1","pop2","cs.matrix","response")
    
    # Assign value to layer
    LAYER<-assign("LAYER",value=dat$cs.matrix)
    
    # Fit model
    mod <- lFormula(response ~ LAYER + (1|pop1), data=dat,REML=TRUE)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    Mod <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  }
  } # End sublandscapes not used
  
  #!#!#!#!#!#!
  if(!is.null(sublandscape)){
    if(length(resistance.mat)==1){
      response=genetic.dist
      if(is.null(name)){
        NAME<-gsub(pattern="*_resistances.out","",x=(basename(resistance.mat)))
      }
      mm<-read.table(resistance.mat)[-1,-1]
      m<-length(mm)
      mm <-lower(mm)
      mm <- mm[which(mm!=-1)]
      
      if(is.null(ID)){
        ID<-To.From.ID(POPS=m)
      }
      if(is.null(ZZ)){
        ZZ<-ZZ.mat(ID=ID)
      }
      
      cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
      cs.unscale<-mm
#       dat<-data.frame(ID,cs.matrix=cs.matrix,response=response,sub=sublandscape)
#       colnames(dat) <- c("pop1","pop2","cs.matrix","response","sub")
#       
#       # Assign value to layer
#       LAYER<-assign("LAYER",value=dat$cs.matrix)
      
      Mod <- MLPE.lmm.sub(resistance = cs.matrix,
                          response = response,
                          REML = TRUE,
                          ID = ID,
                          ZZ = ZZ,
                          sub = sublandscape)
      
#       # Fit model
#       mod <- lFormula(response ~ LAYER + (1|pop1), data=dat,REML=TRUE)
#       mod$reTrms$Zt <- ZZ
#       dfun <- do.call(mkLmerDevfun,mod)
#       opt <- optimizeLmer(dfun)
#       Mod <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
    }
    
    
    
    if(length(resistance.mat)>1){
      response=genetic.dist
      if(is.null(name)){
        stop("Output file 'name' must be specified!!!")
      }
      NAME<-name
      mm<-lower(as.matrix(resistance.mat))
      m<-attr(resistance.mat,"Size")
      mm <- mm[which(mm!=-1)]
      
      if(is.null(ID)){
        ID<-To.From.ID(POPS=m)
      }
      if(is.null(ZZ)){
        ZZ<-ZZ.mat(ID=ID)
      }
      cs.matrix<-scale(mm,center=TRUE,scale=TRUE)
      cs.unscale<-mm
#       dat<-data.frame(ID,cs.matrix=cs.matrix,response=response)
#       colnames(dat) <- c("pop1","pop2","cs.matrix","response")
#       
#       # Assign value to layer
#       LAYER<-assign("LAYER",value=dat$cs.matrix)
#       
#       # Fit model
#       mod <- lFormula(response ~ LAYER + (1|pop1), data=dat,REML=TRUE)
#       mod$reTrms$Zt <- ZZ
#       dfun <- do.call(mkLmerDevfun,mod)
#       opt <- optimizeLmer(dfun)
#       Mod <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
       
      Mod <- MLPE.lmm.sub(resistance = cs.matrix,
                          response = response,
                          REML = TRUE,
                          ID = ID,
                          ZZ = ZZ,
                          sub = sublandscape)
    }
  } # End sublandscapes used
  #!##!##
  # Make diagnostic plots
  if(type=="categorical"){
    tiff(filename = paste0(plot.dir,NAME,"_DiagnosticPlots.tif"), 
         width = 279, height = 215, units = "mm", 
         compression = c("lzw"),
         bg = "white", res = 300)
    par(mfrow=c(2,1),
        oma = c(0,4,0,0) + 0.1,
        mar = c(4,4,1,1) + 0.1)  
    hist(residuals(Mod),xlab="Residuals",main="")
    qqnorm(resid(Mod),main="")
    qqline(resid(Mod))
    dev.off()
    par(mfrow=c(1,1))  
  } else {
    tiff(filename = paste0(plot.dir,NAME,"_DiagnosticPlots.tif"), 
         width = 279, height = 215, units = "mm", 
         compression = c("lzw"),
         bg = "white", res = 300)
    par(mfrow=c(2,2),
        oma = c(0,4,0,0) + 0.1,
        mar = c(4,4,1,1) + 0.1)
    plot(response~cs.unscale,xlab=XLAB,ylab=YLAB)
    abline(lm(response~cs.unscale))
    plot(residuals(Mod)~cs.unscale,xlab=XLAB,ylab="Residuals")
    abline(lm(residuals(Mod)~cs.unscale))
    hist(residuals(Mod),xlab="Residuals",main="")
    qqnorm(resid(Mod),main="")
    qqline(resid(Mod))
    dev.off()
    par(mfrow=c(1,1))   
  }  
}
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#
# Function to bundle input parameters

#' Prepare and bundle input CIRCUITSCAPE model parameters
#' 
#' This function will prepare objects needed for running optimization functions
#' 
#' @param n.POPS The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param CS_Point.File The path to the Circuitscape formatted point file. See Circuitscape documentation for help.
#' @param CS.program The path to the CIRCUITSCAPE executable file (cs_run.exe) on a Windows PC. See details below. 
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
#' @param pairs_to_include Default is NULL. If you wish to use the advanced CIRCUITSCAPE setting mode to include or exclude certain pairs of sample locations, provide the path to the properly formatted "pairs_to_include.txt" file here. Currently only "include" method is supported.
#' @param sublandscape Default is NULL. If using the "pairs_to_include" advanced option in CIRCUITSCAPE, then this  must be included. This should be a vector of values specifying the number of individuals included in each sublandscape. 
# @param platform What computing platform are you using ("pc", "other"). This code has only been tested on Windows PC!!!
#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage CS.prep(n.POPS, response, CS_Point.File, CS.program, Neighbor.Connect, pairs_to_include, sublandscape)
#' @details \code{CS.program} Example of path to CIRCUITSCAPE executible on Windows: 
#' 
#' '"C:/Program Files/Circuitscape/cs_run.exe"'
#'
#' ***NOTE: Double quotation used***
#' This is the current default for \code{CS.program}, but the directory may need to be changed depending upon your installation of CIRCUITSCAPE

CS.prep <- function(n.POPS, response=NULL,CS_Point.File,CS.program='"C:/Program Files/Circuitscape/cs_run.exe"',Neighbor.Connect=8, pairs_to_include=NULL, sublandscape=NULL){
  CS.exe_Test <- gsub("\"", "", CS.program)
  # Error messages
  if(!file.exists(CS_Point.File)) { stop( "The specified CS_Point.File does not exist" ) }
  if(!file.exists(gsub("\"", "", CS.program))) { stop( "The specified path to 'cs_run.exe' is incorrect" ) }
  
  if(grepl(".asc", x = CS_Point.File)){
    CS_grid<-raster<-raster(CS_Point.File)
    CS_Point.txt <- rasterToPoints(x = CS_grid)
    site<-CS_Point.txt[ , 3]
    CS_Point.txt<-data.frame(site,CS_Point.txt[,c(1,2)])
    CS_Point.txt <- CS_Point.txt[order(site),]
    CS_Point.File <- sub(".asc",".txt", x = CS_Point.File)
    write.table(CS_Point.txt,file = CS_Point.File,col.names = F,row.names = F)
  }
  if(!is.null(response)) {TEST.response <- (is.vector(response) || ncol(response)==1)
                             if(TEST.response==FALSE) {stop("The object 'response' is not in the form of a single column vector")}
  }
  platform="pc"
  
  if(!is.null(pairs_to_include)){
    if(!file.exists(pairs_to_include)) { stop( "The specified pairs_to_include file does not exist") }
        toMatch <- c("min", "max")
        if(grep(read.table(file = pairs_to_include,header = F,sep = "\t")[1,1],pattern = paste(toMatch,collapse = "|"))){
          # Function to make list of observations to include
          Min.MAX <- read.table(file = pairs_to_include,header = F,sep = "\t")[c(1,2),c(1,2)]
          MIN <- Min.MAX[which(Min.MAX[,1]=="min"),2]
          MAX <- Min.MAX[which(Min.MAX[,1]=="max"),2]
          PTI <- read.table(file = pairs_to_include,header = F,sep = "\t")[-c(1:3),]
          site <- PTI[,1]
          PTI <- as.matrix(PTI[,-1])
          
          p_t_i <- list()
          count <- 0
          for(i in 1:(ncol(PTI)-1)){
            for(j in (i+1):ncol(PTI)){
                if(PTI[j,i]>=MIN && PTI[j,i]<=MAX){
                count <- count +1                  
                p_t_i[[count]] <- data.frame(i,j)                 
              } # close if statement
            } # close j loop           
          } # close i loop
          ID <- plyr::ldply(p_t_i,.fun = identity)
          colnames(ID) <- c("pop1","pop2")
#           ID<-arrange(tmp2,as.numeric(pop1),as.numeric(pop2))
          n1 <- table(ID$pop1)[[1]]
          p1<-ID[n1,1]; p2<-ID[n1,2]
          ID[n1,1]<-p2; ID[n1,2]<-p1
          ID$pop1 <- factor(ID$pop1)
          ID$pop2 <- factor(ID$pop2)
          sub <-  sub.pair(sublandscape)
          if(nrow(ID)!=nrow(sub)) warning("The 'sub-landscape' vector is not equal to the number of pairwise comparisons. Please confirm that the number of individuals included in each sublandscape is correct.")
          suppressWarnings(ZZ<-ZZ.mat_sub(ID,sub))

        } # close function                
  } # close pairs to include statement
  
  # Make to-from population list
  if(!exists(x = "ID")){
  ID<-To.From.ID(n.POPS)
  suppressWarnings(ZZ<-ZZ.mat(ID))  
  }
  list(ID=ID,ZZ=ZZ,response=response,CS_Point.File=CS_Point.File,CS.program=CS.program,Neighbor.Connect=Neighbor.Connect,n.POPS=n.POPS,platform=platform,pairs_to_include=pairs_to_include,sub=sub)
}

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#
# Data processing for analyzing multiple layers simultaneously

#' Create R object with genetic algorithm optimization settings
#' 
#' This function prepares and compiles objects and commands for optimization with the GA package
#' 
#' @param ASCII.dir Directory containing all raster objects to be optimized. If optimizing using least cost paths, a RasterStack or RasterLayer object can be specified.
#' @param Results.dir If a RasterStack is provided in place of a directory containing .asc files for ASCII.dir, then a directory to export optimization results must be specified. It is critical that there are NO SPACES in the directory, as this will cause the function to fail.
#' @param min.cat The minimum value to be assessed during optimization of of categorical resistance surfaces (Default = 1e-04)
#' @param max.cat The maximum value to be assessed during optimization of of categorical resistance surfaces (Default = 2500)
#' @param max.cont The maximum value to be assessed during optimization of of continuous resistance surfaces (Default = 2500)
#' @param cont.shape A vector of hypothesized relationships that each continuous resistance surface will have in relation to the genetic distance reposnse (Default = NULL; see details)
#' @param pop.mult Value will be multiplied with number of parameters in surface to determine 'popSize' in GA. By default this is set to 15.
#' @param percent.elite Percent used to determine the number of best fitness individuals to survive at each generation ('elitism' in GA). By default the top 5\% individuals will survive at each iteration.
#' @param type Default is "real-valued"
#' @param population Default is "gareal_Population" from GA
#' @param selection Default is "gareal_lsSelection" from GA
#' @param mutation Default is "gareal_raMutation" from GA
#' @param pcrossover Probability of crossover. Default = 0.85
#' @param pmutation Probability of mutation. Default = 0.125
#' @param crossover Default = "gareal_blxCrossover". This crossover method greatly improved optimization during preliminary testing
#' @param maxiter Maximum number of iterations to run before the GA search is halted (Default = 1000)
#' @param pop.size Number of individuals to create each generation
#' @param parallel A logical argument specifying if parallel computing should be used (TRUE) or not (FALSE, default) for evaluating the fitness function. You can also specifiy the number of cores to use. Parallel processing currently only works when optimizing using least cost paths. It will fail if used with CIRCUITSCAPE, so this is currently not an option.
#' @param run Number of consecutive generations without any improvement in AICc before the GA is stopped (Default = 25)
#' @param keepBest A logical argument specifying if best solutions at each iteration should be saved (Default = TRUE)
#' @param Min.Max Define whether the optimization function should be minimized ('min') or maximized ('max' = Default). Optimization with \code{ga} maximizes the objective criteria
#' @param seed Integer random number seed to replicate \code{ga} optimization
#' @param quiet Logical. If TRUE, AICc and step run time will not be printed to the screen after each step. Only \code{ga} summary information will be printed following each iteration. Default = FALSE
#' @return An R object that is a required input into optimization functions
#' 
#' @details Only files that you wish to optimize, either in isolation or simultaneously, should be included in the specified \code{ASCII.dir}. If you wish to optimize different combinations of surfaces, different directories contaiing these surfaces must be created.
#' 
#' \code{cont.shape} can take values of "Increase", "Decrease", or "Peaked". If you believe a resistance surface is related to your reposnse in a particular way, specifying this here may decrease the time to optimization. \code{cont.shape} is used to generate an initial set of parameter values to test during optimization. If specified, a greater proportion of the starting values will include your believed relatiosnship. If unspecified (the Default), a completely random set of starting values will be generated.
#' 
#' It is recommended to first run GA optimization with the default settings

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage GA.prep(ASCII.dir,
#' Results.dir=NULL,
#' min.cat=1e-04,
#' max.cat=2500,
#' max.cont=2500,
#' cont.shape=NULL,
#' Min.Max="max",
#' pop.mult = 15,
#' percent.elite = 0.05,
#' type= "real-valued",
#' pcrossover=0.85,
#' pmutation=0.125,
#' maxiter=1000,
#' run=25,
#' keepBest=TRUE,
#' population = gaControl(type)$population,
#' selection = gaControl(type)$selection,
#' crossover="gareal_blxCrossover",
#' mutation = gaControl(type)$mutation,
#' pop.size = NULL,
#' parallel = FALSE,
#' seed = NULL,
#' quiet = FALSE)

GA.prep<-function(ASCII.dir,
                  Results.dir = NULL,                  
                  min.cat = 0.0001,
                  max.cat = 2500, 
                  max.cont = 2500,
                  cont.shape = NULL,
                  Min.Max ='max',
                  pop.mult = 15,
                  percent.elite = 0.05,
                  type = "real-valued",
                  pcrossover = 0.85,
                  pmutation = 0.125,
                  maxiter = 1000,
                  run = 25,
                  keepBest = TRUE,
                  population = gaControl(type)$population,
                  selection = gaControl(type)$selection,
                  crossover="gareal_blxCrossover",
                  mutation = gaControl(type)$mutation,
                  pop.size = NULL,
                  parallel = FALSE,
                  seed = NULL,
                  quiet = FALSE) { 
  
  if(!is.null(Results.dir)) {TEST.dir <- !file_test("-d",Results.dir)
                             if(TEST.dir==TRUE) {stop("The specified 'Results.dir' does not exist")}}
   
  if((class(ASCII.dir)[1]=='RasterStack' | class(ASCII.dir)[1]=='RasterLayer') & is.null(Results.dir)){
    warning(paste0("'Results.dir' was not specified. Results will be exported to ", getwd()))
    Results.dir<-getwd()
  }
  
  if(class(ASCII.dir)[1]!='RasterStack' & is.null(Results.dir)){
    Results.dir<-ASCII.dir
  }
    
  if(class(ASCII.dir)[1]=='RasterStack' | class(ASCII.dir)[1]=='RasterLayer'){
     r<-ASCII.dir
     names <- names(r)
     n.layers <-length(names)     
  } else {
     ASCII.list <-list.files(ASCII.dir,pattern="*.asc", full.names=TRUE) # Get all .asc files from directory
     if(length(ASCII.list)==0) {stop("There are no .asc files in the specified 'ASCII.dir")}
     r <- stack(lapply(ASCII.list,raster))
     names <- gsub(pattern="*.asc","",x=(list.files(ASCII.dir,pattern="*.asc")))
     n.layers <-length(ASCII.list) 
  }
   
  if("Results"%in%dir(Results.dir)==FALSE) dir.create(file.path(Results.dir, "Results")) 
  #   dir.create(file.path(ASCII.dir, "Results"),showWarnings = FALSE)
  Results.DIR<-paste0(Results.dir, "Results/")
  if("tmp"%in%dir(Results.dir)==FALSE) dir.create(file.path(Results.dir, "tmp")) 
  #   dir.create(file.path(Results.dir, "tmp"),showWarnings = FALSE)
  Write.dir <-paste0(Results.dir,"tmp/") 
  if("Plots"%in%dir(Results.DIR)==FALSE) dir.create(file.path(Results.DIR, "Plots")) 
  #   dir.create(file.path(Results.dir, "tmp"),showWarnings = FALSE)
  Plots.dir <-paste0(Results.DIR,"Plots/") 
      
  # Determine total number of parameters and types of surfaces included
  parm.type<-data.frame()
  min.list <- list()
  max.list <-list()
  SUGGESTS <- list()
  for(i in 1:n.layers){
    n.levels<-length(unique.rast(r[[i]]))
    if (n.levels <=15){
      Level.val <- unique.rast(r[[i]])      
      parm.type[i,1]<-"cat"
      parm.type[i,2]<-n.levels 
      parm.type[i,3]<-names[i]
      min.list[[i]]<-c(1,rep(min.cat,(n.levels-1))) 
      max.list[[i]]<-c(1,rep(max.cat,(n.levels-1)))
      
    } else {
      parm.type[i,1]<-"cont"
      parm.type[i,2]<-3
      parm.type[i,3]<-names[i]
      min.list[[i]]<-c(1,.001,.001) # eq, shape/gaus.opt, max, gaus.sd
      max.list[[i]]<-c(9.99,15,max.cont)     
    }
  }
  
  colnames(parm.type)<-c("type","n.parm","name")   
  parm.index <- c(0,cumsum(parm.type$n.parm)) 
  ga.min <- unlist(min.list)
  ga.max <- unlist(max.list)
  surface.type <- parm.type$type
  
  if(is.null(pop.size)){
  if (length(ga.min)<10){
    pop.size <- min(c(15*length(ga.min),100))
  } else {
    pop.size <- 10*length(ga.min)
  }
  }
  
  for(i in 1:length(surface.type)){
    if (surface.type[i]=="cat"){
      SUGGESTS[[i]] <- sv.cat(levels=parm.type[i,2],pop.size=pop.size,min.cat,max.cat)
      
    } else if (exists("cont.shape") && length(cont.shape>0)){
      SUGGESTS[[i]] <- sv.cont.nG(cont.shape[1],pop.size=pop.size,max.cont)
      cont.shape<-cont.shape[-1]
    } else {
      SUGGESTS[[i]] <- sv.cont.nG("NA",pop.size=pop.size,max.cont)
    }
  }
  SUGGESTS <-matrix(unlist(SUGGESTS), nrow=nrow(SUGGESTS[[1]]), byrow=F)
  
  list(parm.index=parm.index,ga.min=ga.min,ga.max=ga.max,surface.type=surface.type,parm.type=parm.type,Resistance.stack=r,n.layers=n.layers,layer.names=names,pop.size=pop.size, min.list=min.list,max.list=max.list, SUGGESTS=SUGGESTS,ASCII.dir=ASCII.dir, Results.dir=Results.DIR, Write.dir=Write.dir,Plots.dir=Plots.dir,type= type, pcrossover=pcrossover, pmutation=pmutation, crossover=crossover, maxiter=maxiter, run=run, keepBest=keepBest, population=population,selection=selection,mutation=mutation, parallel=parallel,pop.mult = pop.mult, percent.elite = percent.elite,Min.Max=Min.Max, seed=seed, quiet = quiet)  
  
}
#!##!##!##!##!##!##!##!##!##!##!##!##
#' Make a vector of the lower half of a square distance matrix
#' 
#' This function will prepare and compile objects and commands for optimization with the GA package
#' 
#' @param matrix Square distance matrix with no row names.
#' @return A vector of the lower half of the matrix
#' 
#' @details This is a convenience function to obtain the lower half of a matrix, which is required as input for several other functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>

lower<-function(matrix){
  if(is.vector(matrix)==TRUE || dim(matrix)[1]!=dim(matrix)[2]) {warning("Must provide square distance matrix with no column or row names")}
  lm<-matrix[lower.tri(matrix)]
  return(lm)
}
#!##!##!##!##!##!##!##!##!#

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!###
#!##!##!##!# OTHER NECESSARY FUNCTIONS  #!##!##!##!##!##!##!#
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##

# FUNCTIONS
OPTIM.DIRECTION <- function(x){
  OPTIM<-ifelse(x=='max',-1,1)
  return(OPTIM)
}


Cont.Param <- function(PARM){
  df<-data.frame(PARM[1],PARM[2])
  colnames(df)<-c("shape_opt","max");row.names(df)<-NULL
  return(df) 
}

# read.matrix<-function(cs.matrix){  m<-read.table(cs.matrix)[-1,-1]
#                                    m[lower.tri(m)]}

read.matrix<-function(cs.matrix){  lower(read.table(cs.matrix)[-1,-1])
                                   }

read.matrix2<-function(cs.matrix){  m<-read.table(cs.matrix)[-1,-1]
}
# Make to-from population list
To.From.ID<-function(POPS){
  tmp <- matrix(nrow=POPS,ncol=POPS)
  dimnames(tmp) <- list( 1:POPS, 1:POPS)  
  tmp2 <- as.data.frame( which( row(tmp) < col(tmp), arr.ind=TRUE))  
  tmp2[[2]] <-dimnames(tmp)[[2]][tmp2$col]
  tmp2[[1]] <-dimnames(tmp)[[2]][tmp2$row]
  colnames(tmp2)<-c("pop1","pop2")
  as.numeric(tmp2$pop1);as.numeric(tmp2$pop2)
  ID<-arrange(tmp2,as.numeric(pop1),as.numeric(pop2))
  #   ID<-tmp2[with(tmp2, order(pop1, pop2)), ]
  p1<-ID[POPS-1,1]; p2<-ID[POPS-1,2]
  ID[POPS-1,1]<-p2; ID[POPS-1,2]<-p1
  ID$pop1 <- factor(ID$pop1)
  ID$pop2 <- factor(ID$pop2)
  return(ID)
}

# Create ZZ matrix for mixed effects model
ZZ.mat <- function(ID) {
  Zl <- lapply(c("pop1","pop2"), function(nm) Matrix::fac2sparse(ID[[nm]],"d", drop=FALSE))
  ZZ <- Reduce("+", Zl[-1], Zl[[1]])
  return(ZZ)
}

ZZ.mat_sub <- function(ID,sub) {
#   sub[,1] <- factor(sub[,1])
  ID <- cbind(ID,sub)
  names(ID)<- c("pop1",'pop2',"sub")
  Zl <- lapply(c("pop1","pop2"), function(nm) Matrix::fac2sparse(ID[[nm]],"d", drop=FALSE))
  Zl2 <- lapply(c("sub"), function(nm) Matrix::fac2sparse(ID[[nm]],"d", drop=FALSE))
  ZZ <- Reduce("+", Zl[-1], Zl[[1]])
  ZZ.sub <- Matrix::rBind(ZZ,Zl2[[1]],deparse.level = 1)
  return(ZZ.sub)
}

# Rescale function
SCALE.vector <-function(data,MIN,MAX){Mn=min(data)
                                      Mx=max(data)
                                      (MAX-MIN)/(Mx-Mn)*(data-Mx)+MAX}

# Define scaling function
# This will rescale from 1 to specified MAX
SCALE <-function(data,MIN,MAX){Mn=cellStats(data,stat='min')
                               Mx=cellStats(data,stat='max')
                               (MAX-MIN)/(Mx-Mn)*(data-Mx)+MAX}

# # Function to write .ini file for Circuitscape 
# write.CS_3.5.8 <- function(BATCH,OUT,HABITAT,LOCATION.FILE,CONNECTION,MAP="write_cum_cur_map_only=False"){
# sink(BATCH)
# cat("[Options for advanced mode]
# ground_file_is_resistances=True
# source_file=(Browseforacurrentsourcefile)
# remove_src_or_gnd=keepall
# ground_file=(Browseforagroundpointfile)
# use_unit_currents=False
# use_direct_grounds=False
# 
# [Calculation options]
# low_memory_mode=False
# solver=cg+amg
# print_timings=True
# 
# [Options for pairwise and one-to-all and all-to-one modes]
# included_pairs_file=None
# point_file_contains_polygons=False
# use_included_pairs=False")
# cat("\n")
# cat(LOCATION.FILE)
# cat("\n")
# 
# cat("
# [Output options]")
# cat("\n")
# cat(MAP)
# cat("\n")
# cat("log_transform_maps=False
# set_focal_node_currents_to_zero=False
# write_max_cur_maps=False
# write_volt_maps=False
# set_null_currents_to_nodata=True
# set_null_voltages_to_nodata=True
# compress_grids=False
# write_cur_maps=False")
# cat("\n")
# cat(OUT)
# cat("\n")
# cat("\n")
# cat("[Shortcircuit regions(aka polygons)]
# use_polygons=False
# polygon_file=(Browse for a short-circuit region file)
# 
# [Connection scheme for raster habitat data]")
# cat("\n")
# cat(CONNECTION)
# cat("\n")
# cat("connect_using_avg_resistances=True
# 
# [Habitat raster or graph]
# habitat_map_is_resistances=True")
# cat("\n")
# cat(HABITAT)
# cat("\n")
# cat("\n")
# cat("[Options for one-to-all and all-to-one modes]
# use_variable_source_strengths=False
# variable_source_file=None
# 
# [Version]
# version = 3.5.8
# 
# [Maskfile]
# use_mask=False
# mask_file=None
# 
# [Circuitscape mode]
# data_type=raster
# scenario=pairwise")
# sink()
# }


write.CS_4.0 <- function(BATCH,OUT,HABITAT,LOCATION.FILE,CONNECTION,CURRENT.MAP="write_cur_maps = False", MAP="write_cum_cur_map_only = False",PARALLELIZE="parallelize = False",CORES="max_parallel = 0", PAIRS_TO_INCLUDE, PAIRS){
sink(BATCH)
cat("[Options for advanced mode]
ground_file_is_resistances = True
remove_src_or_gnd = rmvsrc
ground_file = (Browse for a raster mask file)
use_unit_currents = False
source_file = (Browse for a raster mask file)
use_direct_grounds = False

[Mask file]
mask_file = (Browse for a raster mask file)
use_mask = False

[Calculation options]
low_memory_mode = False")
cat("\n")
cat(PARALLELIZE)
cat("\n")
cat(CORES)
cat("\n")
cat("
solver = cg+amg
print_timings = False
preemptive_memory_release = False
print_rusages = False

[Short circuit regions (aka polygons)]
polygon_file = (Browse for a short-circuit region file)
use_polygons = False

[Options for one-to-all and all-to-one modes]
use_variable_source_strengths = False
variable_source_file = (Browse for a short-circuit region file)

[Output options]
set_null_currents_to_nodata = True
set_focal_node_currents_to_zero = False
set_null_voltages_to_nodata = True
compress_grids = False
write_volt_maps = False")
cat("\n")
cat(CURRENT.MAP)
cat("\n")
cat(OUT)
cat("\n")
cat(MAP)
cat("\n")
cat("
log_transform_maps = False
write_max_cur_maps = False

[Version]
version = 4.0.5

[Options for reclassification of habitat data]
reclass_file = (Browse for file with reclassification data)
use_reclass_table = False

[Logging Options]
log_level = INFO
log_file = None
profiler_log_file = None
screenprint_log = False

[Options for pairwise and one-to-all and all-to-one modes]")
cat("\n")
cat(PAIRS_TO_INCLUDE)
cat("\n")
cat(PAIRS)
cat("\n")
cat(LOCATION.FILE)
cat("\n")
cat("\n")

cat("
[Connection scheme for raster habitat data]
connect_using_avg_resistances = True")
cat("\n")
cat(CONNECTION)
cat("\n")
cat("\n")
cat("
[Habitat raster or graph]
habitat_map_is_resistances = True")
cat("\n")
cat(HABITAT)
cat("\n")
cat("\n")
cat("
[Circuitscape mode]
data_type = raster
scenario = pairwise")
sink()
}

#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##
# Sample values for suggests
sv.cat <-function(levels,pop.size,min,max){
  cat.starts<-matrix(nrow=pop.size,ncol=levels)
  for(r in 1:pop.size){
    L<-list()
    for (i in 1:levels){
      if(runif(1)<.5){
        z<-runif(1)
      } else {
        z<-runif(1,min,max)
      }
      L[[i]]<-z
    } 
    #   uz<-unlist(L)
    cat.starts[r,]<-(unlist(L))    
  }
  cat.starts[,1]<-1
  return(cat.starts)
}
#!##!##!##!##!##!##!##!##!##
# No Gaussian distribution
sv.cont.nG <-function(direction,pop.size,max){
  inc<-c(1,3)
  dec<-c(7,5)
  peak<-c(2,4,6,8)
  L<-list()
  cont.starts<-matrix(nrow=pop.size,ncol=3)  
  for(r in 1:pop.size){
    if(runif(1)<.5 && direction=="Increase"){
      #       z1<-c(sample(inc,1)
      z<-Increase.starts.nG(sample(inc,1))
    } else if(runif(1)<.5 && direction=="Decrease") {
      z<-c(sample(dec,1),runif(1,.01,10),runif(1,1,max))
    } else if (runif(1)<.5 && direction=="Peaked") {
      z<-c(sample(peak,1),runif(1,.01,10),runif(1,1,max))
    } else {
      z<-c(runif(1,1,9.99),runif(1,.01,10),runif(1,1,max))
    }
    cont.starts[r,]<-z
  } 
  cont.starts
}

#!##!##!##!##!##!##!##!##!##
Increase.starts<-function(x){
  if(x==1){
    z<-c(x,runif(1,.01,10),runif(1,1,10),1)
  } else {
    z<-c(x,runif(1,.01,10),runif(1,1,100),1)
  }
}

Increase.starts.nG<-function(x){
  if(x==1){
    z<-c(x,runif(1,.01,10),runif(1,1,10))
  } else {
    z<-c(x,runif(1,.01,10),runif(1,1,100))
  }
}
#!##!##!##!##!##!##!##!##!#
unique.rast<-raster::unique

get.EQ <-function(equation){   # Apply specified transformation
  if(is.numeric(equation)){
  equation=floor(equation)
  if(equation==1){
    EQ <- "Inverse-Reverse Monomolecular"
    
  } else if(equation==5){
    EQ <- "Reverse Monomolecular"        
    
  } else if(equation==3){
    EQ <- "Monomolecular"
    
  } else if (equation==7) {
    EQ <- "Inverse Monomolecular"  	    
    
  } else if (equation==8) {
    EQ <- "Inverse Ricker"  
    
  } else if (equation==4) {
    EQ <- "Ricker"
    
  } else if (equation==6) {
    EQ <- "Reverse Ricker"
  
  } else if (equation==2) {
    EQ <- "Inverse-Reverse Ricker"
  
  } else {
    EQ <- "Distance"  	
  }
 
  (EQ)
  } else {
  if(equation=="Inverse-Reverse Monomolecular"){
    EQ <- 1
    
  } else if(equation=="Reverse Monomolecular"){
    EQ <- 5        
    
  } else if(equation=="Monomolecular"){
    EQ <- 3
    
  } else if (equation=="Inverse Monomolecular") {
    EQ <- 7
    
  } else if (equation=="Inverse Ricker") {
    EQ <- 8
    
  } else if (equation=="Ricker") {
    EQ <- 4
    
  } else if (equation=="Reverse Ricker") {
    EQ <- 6
  
  } else if (equation=="Inverse-Reverse Ricker") {
    EQ <- 2
  
  } else {
    EQ <- 9  	
  }
 
  (EQ)  
  }
}

Result.txt <- function(GA.results, GA.inputs, method, Run.Time){
  summary.file<-paste0(GA.inputs$Results.dir,"Multisurface_Optim_Summary.txt")
  AICc<--GA.results@fitnessValue*-1
  AICc<-round(AICc,digits=4)
  ELITE<-floor(GA.inputs$percent.elite*GA.inputs$pop.size)
#   mlpe.results<-MLPE.lmm_coef(GA.inputs$Results.dir,genetic.dist=CS.inputs$response)
  
sink(summary.file)
cat(paste0("Summary from multisurface optimization run conducted on ",Sys.Date()),"\n")
cat("\n")
cat(paste0("Optimized using:",method),"\n")
cat("\n")
cat(paste0("Surfaces included in optimization:"),"\n")
cat(GA.inputs$parm.type$name,"\n")
cat("\n")
cat("Genetic Algorithm optimization settings:")
cat("\n")
cat(paste0("Type of genetic algorithm used: ", GA.results@type),"\n")
cat(paste0("popSize at each iteration: ", GA.results@popSize),"\n")
cat(paste0("Maximum number of iterations: ", GA.results@maxiter),"\n")
cat(paste0("Number individuals retained each generation (elistism): ", ELITE),"\n")
cat(paste0("Crossover probability: ", GA.results@pcrossover),"\n")
cat(paste0("Mutation probability: ", GA.results@pmutation),"\n")
cat("\n")
cat(paste0("The Genetic Algorithm completed after ",GA.results@iter," iterations"),"\n")
cat("\n")
cat(paste0("Minimum AICc: ",-1*AICc),"\n")
cat("\n")
cat(paste0("Optimized values for each surface:"),"\n")
cat(GA.results@solution,"\n")
cat("\n")
cat(paste0("Optimization took ",Run.Time," seconds to complete"),"\n")
sink()
}
 
#!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!##!#

get.name<-function(x){
   nm <-deparse(substitute(x))
   return(nm)
   }

# Function to determine identity of pairwise sublandscape elements
sub.pair <- function(x){
  sub <- data.frame()
 for(i in seq_along(x)){
    m <- matrix(0,x[i],x[i])
    N <- length(lower(m))   
    N.out <- as.data.frame(rep(i,times = N))
    sub <- rbind(sub,N.out)    
 }
 colnames(sub) <- "sub"
 return(sub)
}
#!##!##!##!##!##!##!##!##!#
Monomolecular <- function(r, parm){
  parm[3]*(1-exp(-1*r/parm[2]))+1 # Monomolecular
}

Inv.Monomolecular <- function(r, parm){
  if(class(r)=="RasterLayer") {
    R <- parm[3]*(exp(-1*r/parm[2])) 
    (R <- (R-cellStats(R,stat = "min"))+1)  
  } else {
    R <- parm[3]*(exp(-1*r/parm[2])) 
    (R <- (R-min(R))+1)
     }
}

Inv.Rev.Monomolecular <- function(r, parm){
  if(class(r)=="RasterLayer") {
    rev.rast <- SCALE((-1*r),0,10)
    Inv.Monomolecular(rev.rast,parm)
  } else {
    rev.rast <- SCALE.vector((-1*r),0,10)
    Inv.Monomolecular(rev.rast,parm)
  }
}

Rev.Monomolecular <- function(r, parm){
  if(class(r)=="RasterLayer") {
    rev.rast <- SCALE((-1*r),0,10)
    Monomolecular(rev.rast,parm)
  } else {
    rev.rast <- SCALE.vector((-1*r),0,10)
    Monomolecular(rev.rast,parm)
  }
}


Ricker <- function(r,parm){
  parm[3]*r*exp(-1*r/parm[2])+1 # Ricker
  }

Inv.Ricker <- function(r,parm){
  if(class(r)=="RasterLayer") {
    R <- (-1*parm[3])*r*exp(-1*r/parm[2])-1 # Ricker
    R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
  } else {    
    R <- (-1*parm[3])*r*exp(-1*r/parm[2])-1 # Ricker
    R <- SCALE.vector(R,MIN=abs(max(R)),MAX=abs(min(R))) # Rescale
  }
}

Inv.Rev.Ricker <- function(r,parm){
  if(class(r)=="RasterLayer") {
    rev.rast <- SCALE((-1*r),0,10)
    Inv.Ricker(rev.rast,parm)
  } else {
    rev.rast <- SCALE.vector((-1*r),0,10)
    Inv.Ricker(rev.rast,parm)
  }
}

Rev.Ricker <- function(r,parm){
  if(class(r)=="RasterLayer") {
    rev.rast <- SCALE((-1*r),0,10)
    Ricker(rev.rast,parm)
  } else {
    rev.rast <- SCALE.vector((-1*r),0,10)
    Ricker(rev.rast,parm)
  }
}

#!##!##!##!###
# TEST
# plot(r)
# Mono <- Monomolecular(r = r,parm = parm); plot(Mono)
# Rev.Mono <- Rev.Monomolecular(r,parm); plot(Rev.Mono)
# Inv.Rev.Mono <- Inv.Rev.Monomolecular(r,parm); plot(Inv.Rev.Mono)
# Inv.Mono <- Inv.Monomolecular(r,parm); plot(Inv.Mono)
# 
# Rick <- Ricker(r = r,parm = parm); plot(Rick)
# Rev.Rick <- Rev.Ricker(r,parm); plot(Rev.Rick)
# Inv.Rev.Rick <- Inv.Rev.Ricker(r,parm); plot(Inv.Rev.Rick)
# Inv.Rick <- Inv.Ricker(r,parm); plot(Inv.Rick)

