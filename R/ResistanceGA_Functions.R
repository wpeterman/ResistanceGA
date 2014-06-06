########################################################################################  
############ CONDUCT GRID SEARCH OF PARAMETER SPACE TO VIEW response SURFACE ########### 
########################################################################################  
#' Conduct grid search of response surface 
#' 
#' Visualize the AICc response surface
#' 
#' @param shape A vector of values for the shape parameter
#' @param max A vector of values for the maximum value parameter
#' @param transformation Transformation to apply. Can be either numeric or character of transformation name
#' @param Resistance An R Raster object, or path to a .asc file
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @usage Grid.Search(shape, max, transformation, Resistance, CS.inputs)
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @return This function will return values that can be plotted to visualize the response surface
#' @details This function will perform a full factorial grid search of the values provided in the shape and max.scale vectors. Depending on the number of values provided for each, and the time it takes to run each iteration, this process may take a while to complete. \cr Suitable values for transformation:\cr
#' \tabular{ll}{
#'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#'    \tab 2 = "Inverse-Reverse Ricker"\cr
#'    \tab 3 = "Monomolecular"\cr
#'    \tab 4 = "Ricker"\cr"
#'    \tab 5 = "Reverse Monomolecular"\cr
#'    \tab 6 = "Reverse Ricker"\cr
#'    \tab 7 = "Inverse Monomolecular"\cr
#'    \tab 8 = "Inverse Ricker"\cr
#'    \tab 9 = "Distance"\cr
#'    }



Grid.Search <- function(shape, max, transformation, Resistance, CS.inputs) {
  if(class(Resistance)[1]!='RasterLayer') {  
    r <- raster(Resistance)
    r <- SCALE(r,0,10)
  } else {    
    r <-SCALE(Resistance,0,10)
  }
   
  GRID <- expand.grid(shape,max)
  RESULTS <- matrix(nrow=nrow(GRID),ncol=3); colnames(RESULTS)<-c("shape","max","AICc")
  EQ<-get.EQ(transformation)

for(i in 1:nrow(GRID)){
  AICc<-Resistance.Optimization_cont.nlm(PARM=log(c(t(GRID[i,]))),Resistance=r,equation=EQ, get.best=FALSE,CS.inputs,Min.Max='min')
  
  results<-as.matrix(cbind(GRID[i,],AICc))
  
  RESULTS[i,]<-results  
}
RESULTS <- data.frame(RESULTS)
Results.mat <- interp(RESULTS$shape,RESULTS$max,RESULTS$AICc,duplicate='strip')

AICc<-RESULTS
colnames(AICc)<-c("shape","max","AICc")
Results.mat<-list(Plot.data=Results.mat,AICc=AICc)

return(Results.mat)
}
########################################################################################  
############ Single command function to execute single surface optimization ############ 
########################################################################################  
#' Single surface optimization
#' 
#' Optimize all surfaces contained in a directory using a genetic algorithm executed with the \code{\link[GA]{ga}} function in the Genetic Algorithms package \pkg{GA}
#' 
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param nlm Logical, if TRUE, the final step of optimization will use nlm to fine-tune parameter estimates. This may lead to overfitting in some cases. Default = FALSE.
#' @return This function optimizes multiple resistance surfaces. Following optimization, several summary objects are created.\cr
#' \enumerate{
#' \item Diagnostic plots of model fit are output to the "Results/Plots" directory that is automatically generated within the folder containing the optimized ASCII files.
#' \item A .csv file with the Maximum Likelihood Population Effects mixed effects model coefficient estimates (MLPE_coeff_Table.csv)
#' \item Three summary .csv files are generated: CategoricalResults.csv, ContinuousResults.csv, & All_Results_AICc.csv. These tables contain AICc values and optimization summaries for each surface.
#' }
#' All results tables are also summarized in a named list ($ContinuousResults, $CategoricalResults, $AICc, $MLPE)
#' @usage SS_optim(CS.inputs, GA.inputs, nlm)
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @export
SS_optim <- function(CS.inputs,GA.inputs, nlm=FALSE){
  RESULTS.cat <- list() # List to store categorical results within
  RESULTS.cont <-list() # List to store continuous results within
  cnt1<-0
  cnt2<-0
  # Optimize each surface in turn
  for (i in 1:GA.inputs$n.layers){
    r<-GA.inputs$Resistance.stack[[i]]
    names(r)<-GA.inputs$layer.names[i]
    
    # Processing of categorical surfaces    
    if (GA.inputs$surface.type[i]=='cat'){
      cnt1 <- cnt1+1    
      r<-GA.inputs$Resistance.stack[[i]]
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
      
      
      df <- data.frame(id=unique.rast(r),single.GA@solution) 
      r <-subs(r,df)
      names(r)<-GA.inputs$layer.names[i]
      
      Run_CS(CS.inputs,GA.inputs,r,EXPORT.dir=GA.inputs$Results.dir)
      
      Diagnostic.Plots(cs.resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),genetic.dist=CS.inputs$response,plot.dir=GA.inputs$Plots.dir)
  
   
      RS <- data.frame(GA.inputs$layer.names[i], -single.GA@fitnessValue,single.GA@solution)
      k=GA.inputs$parm.type$n.parm[i]
      Features <- matrix()
      for(z in 1:(k)){
        feature <- paste0("Feature",z)
        Features[z]<-feature
      }
      
      colnames(RS)<-c("Surface","AICc", Features)
      
      RESULTS.cat[[cnt1]]<-RS
      # Processing of continuous surfaces    
    } else {
      cnt2 <- cnt2+1    
      r<-SCALE(GA.inputs$Resistance.stack[[i]],0,10)
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
      Optim.nlm <-nlm(Resistance.Optimization_cont.nlm, log(start.vals), Resistance=r, equation=single.GA@solution[1],get.best=FALSE,CS.inputs=CS.inputs,Min.Max='min')
      
      OPTIM <- Resistance.Optimization_cont.nlm(PARM=(Optim.nlm$estimate),Resistance=r, equation=single.GA@solution[1],get.best=TRUE,CS.inputs=CS.inputs,Min.Max='min')
      
      Diagnostic.Plots(cs.resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),genetic.dist=CS.inputs$response,plot.dir=GA.inputs$Plots.dir)
      
      Plot.trans(PARM=exp(Optim.nlm$estimate), Resistance=GA.inputs$Resistance.stack[[i]], transformation=EQ, print.dir=GA.inputs$Plots.dir,Name=GA.inputs$layer.names[i])
      
      RS<-data.frame(GA.inputs$layer.names[i],Optim.nlm$minimum,EQ,Cont.Param(exp(Optim.nlm$estimate)))
      colnames(RS) <- c("Surface","AICc","Equation","shape","max")
      RESULTS.cont[[cnt2]] <- RS
     
      } else {
        EQ <-get.EQ(single.GA@solution[1])
        r.tran <- Resistance.tran(transformation=single.GA@solution[1],shape=single.GA@solution[2],max=single.GA@solution[3],r=r) 
        names(r.tran)<-GA.inputs$layer.names[i]
      
        Run_CS(CS.inputs,GA.inputs,r.tran,EXPORT.dir=GA.inputs$Results.dir)
      
        Diagnostic.Plots(cs.resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),genetic.dist=CS.inputs$response,plot.dir=GA.inputs$Plots.dir)
        
        Plot.trans(PARM=single.GA@solution[-1], 
                   Resistance=GA.inputs$Resistance.stack[[i]], 
                   transformation=EQ, 
                   print.dir=GA.inputs$Plots.dir)
          
        RS <- data.frame(GA.inputs$layer.names[i], -single.GA@fitnessValue,get.EQ(single.GA@solution[1]),single.GA@solution[2],single.GA@solution[3])
      colnames(RS) <- c("Surface","AICc","Equation","shape","max")
      RESULTS.cont[[cnt2]] <- RS
      }     
    } # Close if-else    
  } # Close ascii loop
  
  ####################################################
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
  ##################################
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
  write.table(Results.cat,paste0(GA.inputs$Results.dir,"CategoricalResults.csv"),sep=",",col.names=T,row.names=F)
  }
  
  if(ncol(Results.cont)>0){    
    colnames(Results.cont)<-c("Surface","AICc","Equation","shape","max")
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
 
  cat("\n")
  cat("\n")
  write.table(Results.All,paste0(GA.inputs$Results.dir,"All_Results_AICc.csv"),sep=",",col.names=T,row.names=F)
  
  # Get parameter estimates
  MLPE.results<-MLPE.lmm_coef(resist.dir=GA.inputs$Results.dir,genetic.dist=CS.inputs$response,out.dir=GA.inputs$Results.dir)
  
  # Full Results
  if(nrow(Results.cat)>0 & nrow(Results.cont)>0){
    RESULTS<-list(ContinuousResults=Results.cont, CategoricalResults=Results.cat,AICc=Results.All,MLPE=MLPE.results)
  } else if(nrow(Results.cat)<1 & nrow(Results.cont)>0){
    RESULTS<-list(ContinuousResults=Results.cont, CategoricalResults=NULL,AICc=Results.All,MLPE=MLPE.results)
  } else {
    RESULTS<-list(ContinuousResults=NULL, CategoricalResults=NULL,AICc=Results.All,MLPE=MLPE.results)
  }
  return(RESULTS)
  ###############################################################################################################
}


########################################################################################  
############ Single command function to execute mulit surface optimization ############ 
########################################################################################  
#' Simultaneous optimization of multiple resistance surfaces
#' 
#' Optimize multiple resistance surfaces simultaneously using genetic algorithms
#' 
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @return This function optimizes multiple resistance surfaces, returning a Genetic Algorithm (GA) object with summary information. Diagnostic plots of model fit are output to the "Results/Plots" folder that is automatically generated within the folder containing the optimized ASCII files. A text summary of the optimization settings and results is printed to the results folder.
#' @usage MS_optim(CS.inputs, GA.inputs)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
MS_optim<-function(CS.inputs,GA.inputs){
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
                   keepBest=GA.inputs$keepBest,
                   seed = GA.inputs$seed,
                   suggestions=GA.inputs$SUGGESTS,
                   quiet = GA.inputs$quiet) 

  #####  RUN BRENT OPTIMIZATION ####
#   Run second optimization to determine if maximum resistance values should be adjusted
#   GA.opt = Multi.Surface_optim@solution
  Parm.multiplier <- optim(par=1,
                           fn = Max.optim_Brent,
                           method = "Brent",
                           lower = 0,
                           upper = 25,
                           GA.inputs = GA.inputs,
                           CS.inputs = CS.inputs,
                           GA.opt = multi.GA_nG@solution)

PARM<-Parm.multiplier$par
Opt.parm <- GA.opt <- multi.GA_nG@solution
for(i in 1:GA.inputs$n.layers){
    if(GA.inputs$surface.type[i]=="cat"){
      ga.p <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
      parm <- ((ga.p-1)*PARM[1])+1
      Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
      
    } else {
      parm <- GA.opt[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]
      mx<-parm[3]*PARM[1]
      parm[3]<-mx
      Opt.parm[(GA.inputs$parm.index[i]+1):(GA.inputs$parm.index[i+1])]<-parm
    }
}
# ####
multi.GA_nG@solution <- Opt.parm
multi.GA_nG@fitnessValue <- Parm.multiplier$value
  
  RAST<-Combine_Surfaces(PARM=multi.GA_nG@solution,CS.inputs=CS.inputs,GA.inputs=GA.inputs)
  NAME<-paste(GA.inputs$parm.type$name,collapse=".")
  names(RAST)<-NAME
  Run_CS(CS.inputs,GA.inputs,r=RAST,CurrentMap=FALSE,EXPORT.dir=GA.inputs$Results.dir)
  
  Diagnostic.Plots(cs.resistance.mat=paste0(GA.inputs$Results.dir,NAME,"_resistances.out"),genetic.dist=CS.inputs$response,plot.dir=GA.inputs$Plots.dir)
  
  # Get parameter estimates
  MLPE.results<-MLPE.lmm_coef(resist.dir=GA.inputs$Results.dir,genetic.dist=CS.inputs$response,out.dir=GA.inputs$Results.dir)  
  
  Result.txt(GA.results=multi.GA_nG,GA.inputs=GA.inputs, CS.inputs=CS.inputs) 
  return(multi.GA_nG)
}

###############################################################################  
Max.optim_Brent <- function(PARM,CS.inputs,GA.inputs, Min.Max='min', quiet=FALSE, GA.opt){
  t1<-Sys.time()
  
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  response<-CS.inputs$response
  Opt.parm <-vector(length=length(PARM))

  CS_Point.File<-CS.inputs$CS_Point.File
  CS.exe<-CS.inputs$CS.exe
  EXPORT.dir<-GA.inputs$Write.dir
  GA.params<-GA.inputs
  ######
  r <- GA.params$Resistance.stack
  
  for(i in 1:GA.params$n.layers){
    if(GA.params$surface.type[i]=="cat"){
      ga.p <- GA.opt[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]
      parm <- ((ga.p-1)*PARM[1])+1
      Opt.parm[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]<-parm
      df <- data.frame(id=unique.rast(r[[i]]),parm) # Data frame with original raster values and replacement values
      r[[i]] <-subs(r[[i]],df)
      
      r[[i]]<-r[[i]]-(cellStats(x=r[[i]],stat="min"))
      
    } else {
      rast <-SCALE(data=r[[i]],MIN=0,MAX=10)
      parm <- GA.opt[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]
      mx<-parm[3]*PARM[1]
      parm[3]<-mx
      Opt.parm[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]<-parm
      
      # Set equation for continuous surface
      equation <- floor(parm[1]) # Parameter can range from 1-9.99
      
      # Read in resistance surface to be optimized
      SHAPE <-  (parm[2])
      Max.SCALE <- (parm[3])
      
      rick.eq<-(equation==2||equation==4||equation==6||equation==8)
      if(rick.eq==TRUE & SHAPE>6){
        equation<-9
      }
      
      # Apply specified transformation
      if(equation==1){
        SIGN=-1 # Inverse
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE))+SIGN # Monomolecular
        R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Inverse-Reverse Monomolecular"
        
      } else if(equation==5){
        SIGN=1
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE))+SIGN # Monomolecular
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Reverse Monomolecular"        
        
      } else if(equation==3){
        SIGN=1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE))+SIGN # Monomolecular    
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Monomolecular"
        
      } else if (equation==7) {
        SIGN=-1 #Inverse
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE))+SIGN # Monomolecular
        r[[i]] <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Inverse Monomolecular"        
        
      } else if (equation==8) {
        SIGN=-1 #Inverse
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE))+SIGN # Ricker
        r[[i]] <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Inverse Ricker"  
        
      } else if (equation==4) {
        SIGN=1
        r[[i]] <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE))+SIGN #  Ricker
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Ricker"
        
      } else if (equation==6) {
        SIGN=1
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE))+SIGN #  Ricker
        R.vec <- rev(R)
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Reverse Ricker"        
        
      } else if (equation==2) {
        SIGN=-1 # Inverse
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE))+SIGN # Ricker
        R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Inverse-Reverse Ricker"
        
      } else {
        r[[i]] <- (rast*0) #  Cancel layer...set to zero
      } # End if-else  
    } # Close parameter type if-else  
  } # Close layer loop
  
  
  File.name <- "multi_surface"
  
  multi_surface <- sum(r)+1 # Add all surfaces together (+1 for distance)
  if(cellStats(multi_surface,"max")>5e5)  multi_surface<-SCALE(multi_surface,1,5e5) # Rescale surface in case resistance are too high
  
  writeRaster(x=multi_surface,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  
  # Modify and write Circuitscape.ini file
  #############################################################################################  
  BATCH<-paste0(EXPORT.dir,File.name,".ini")        
  OUT<-paste0(paste0("output_file = ",EXPORT.dir), File.name,".out")
  HABITAT<-paste0("habitat_file = ",paste0(EXPORT.dir,File.name,".asc"))
  LOCATION.FILE <- paste0("point_file = ", CS.inputs$CS_Point.File)
  ifelse(CS.inputs$Neighbor.Connect==4,connect<-"True",connect<-"False")
  CONNECTION=paste0("connect_four_neighbors_only=",connect)
  
  #   if(CS.version=='3.5.8'){
  #     write.CS_3.5.8(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,VERSION=VERSION)
  #   } else {
  write.CS_4.0(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION)    
  #   }
  
  ##########################################################################################
  # Run Circuitscape
  CS.exe<-CS.exe
  
  # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
  hidden = TRUE
  
  CS.ini <- paste0(EXPORT.dir,File.name,".ini")
  CS.Run.output<-system(paste(CS.exe, CS.ini), hidden)
  
  #########################################
  # Run mixed effect model on each Circuitscape effective resistance
  
  CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")
  
  # Get AIC statistic for transformed-scaled resistance surface
  cs.matrix<-read.matrix(CS.results)
  cs.matrix<-scale(cs.matrix,center=TRUE,scale=TRUE)
  # cs.matrix2<-round(read.matrix(CS.results),4)
  
  data<-cbind(ID,cs.matrix,response)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=data$cs.matrix)
  
  # Fit model
  mod <- lFormula(response ~ LAYER + (1|pop1), data=data,REML=FALSE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  #    summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  
  k<-max(GA.params$parm.index)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
  
 
  t2 <-Sys.time()
  if(quiet==FALSE){
  cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
  cat(paste0("\t", "AICc = ",round(AICc,4)),"\n","\n")
  }
  
  
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
}

###############################################################################  
############ Create continuous surface response figures ############ 
############################################################################### 
Response.Figs<- function(Optim.input){
  Top.params <- read.csv(file=paste0(Optim.input$Results.dir,"/TopModel_Optimization_Parameters.csv"),header=T)
  dir.create(file.path(Optim.input$Results.dir, "Plots"))
  
  PDF.dir <- paste0(Optim.input$Results.dir,"Plots/")
  for (i in 1:nrow(Top.params)){
    ASCII.file <- list.files(Optim.input$ASCII.dir,pattern=paste0(Top.params[i,1],".asc"),full.names=TRUE)
    if(Top.params[i,5]<1e-5 | Top.params[i,6]<1e-5) {
      cat("\n","\n",paste0("Plotting of ", Top.params[i,1]," could not be completed due to extremely small parameter estimates."),"\n")
      next # Use 'next' command to avoid testing invalid combinations    
    } else {
      PLOT.response(PARM=Top.params[i,c(5,6)],Resistance=ASCII.file,equation=Top.params[i,2],AIC=Top.params[i,7], OutputFolder=PDF.dir)
    }
  }
}

###############################################################################  
############ RUN CIRCUITSCAPE ############ 
###############################################################################  
#' Run CIRCUITSCAPE in R
#' 
#' Execute CS from R
#' 
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param r Accepts two types of inputs. Provide either the path to the raw, untransformed resistance surface file or specify an R raster object
#' @param CurrentMap Logical. If TRUE, the cumulative resistance map will be generated during the CS run (Default = FALSE)
#' @param EXPORT.dir Directory where CS results should be written (Default = GA.inputs$Write.dir, which is a temporary directory for reading/writing CS results)
#' @return Vector of CIRCUITSCAPE resistance distances (lower half of "XXX_resistances.out")
#' @usage Run_CS(CS.inputs, GA.inputs, r, CurrentMap, EXPORT.dir)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Run_CS <- function(CS.inputs,GA.inputs,r,CurrentMap=FALSE,EXPORT.dir=GA.inputs$Write.dir){
if(class(r)[1]!='RasterLayer') {
  R<-raster(r)
  NAME <- basename(r)
  NAME<-sub("^([^.]*).*", "\\1", NAME) 
  names(R)<-NAME
}
  
  if (CurrentMap==FALSE){
    File.name <- basename(r@data@names)
    MAP="write_cum_cur_map_only = False"
    CURRENT.MAP="write_cur_maps = False"
    
  } else {
    File.name <- basename(r@data@names)
    MAP="write_cum_cur_map_only = True"
    CURRENT.MAP="write_cur_maps = 1"
  }
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  response<-CS.inputs$response
  CS_Point.File<-CS.inputs$CS_Point.File
  CS.exe<-CS.inputs$CS.exe
  
  ######
  multi_surface=r
  
  if(cellStats(multi_surface,"max")>5e4)  multi_surface<-SCALE(multi_surface,1,5e4) # Rescale surface in case resistances are too high
  multi_surface <- reclassify(multi_surface, c(-Inf,0, 1))
  
  #     plot(multi_surface)
  
  writeRaster(x=multi_surface,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  
  # Modify and write Circuitscape.ini file
  #############################################################################################  
  BATCH<-paste0(EXPORT.dir,File.name,".ini")        
  OUT<-paste0(paste0("output_file = ",EXPORT.dir), File.name,".out")
  HABITAT<-paste0("habitat_file = ",paste0(EXPORT.dir,File.name,".asc"))
  LOCATION.FILE <- paste0("point_file = ", CS.inputs$CS_Point.File)
  ifelse(CS.inputs$Neighbor.Connect==4,connect<-"True",connect<-"False")
  CONNECTION=paste0("connect_four_neighbors_only=",connect)
  
  #     if(CS.version=='3.5.8'){
  #       write.CS_3.5.8(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION,CURRENT.MAP=CURRENT.MAP,MAP=MAP)
  #     } else {
  write.CS_4.0(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION,MAP=MAP,CURRENT.MAP=CURRENT.MAP)    
  #     }
  
  ##########################################################################################
  # Run Circuitscape
  CS.exe<-CS.exe
  
  # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
  hidden = TRUE
  
  CS.ini <- paste0(EXPORT.dir,File.name,".ini")
  CS.Run.output<-system(paste(CS.exe, CS.ini), hidden)
  # CS.Run.output<-system(paste(CS.exe, CS.ini), hidden,minimized=FALSE) 
  
  #########################################
  # Run mixed effect model on each Circuitscape effective resistance
  
  CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")
  
  (cs.matrix<-read.matrix(CS.results))
  #   cs.matrix<-scale(cs.matrix,center=TRUE,scale=TRUE)
}

###############################################################################  
############ COMBINE RESISTANCE SURFACES--No GAUSSIAN ############ 
###############################################################################  
#' Combine multiple resistance surfaces together
#' 
#' Combine multiple resistance surfaces into new composite surface based on specified parameters
#' 
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Requires a vector with parameters specified in the order of resistance surfaces
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param out Directory to write combined .asc file. Default is GA.inputs$Results.dir
#' @details \code{PARM} is designed to accept the output of \code{MS_optim}. For continuous surfaces, there are three terms: 1) Transformation, 2) shape, and 3) maximum value. Transformation must be provided as a numeric value:\cr
#' \tabular{ll}{
#'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#'    \tab 2 = "Inverse-Reverse Ricker"\cr
#'    \tab 3 = "Monomolecular"\cr
#'    \tab 4 = "Ricker"\cr"
#'    \tab 5 = "Reverse Monomolecular"\cr
#'    \tab 6 = "Reverse Ricker"\cr
#'    \tab 7 = "Inverse Monomolecular"\cr
#'    \tab 8 = "Inverse Ricker"\cr
#'    \tab 9 = "Distance"\cr
#'    }
#' 
#' @return R raster object that is the sum all transformed and/or reclassified resistance surfaces provided
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Combine_Surfaces <- function(PARM,CS.inputs,GA.inputs, out=GA.inputs$Results.dir){
  t1<-Sys.time()
  GA.params<-GA.inputs
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  response<-CS.inputs$response
  #   CS.version<-CS.inputs$CS.version
  CS_Point.File<-CS.inputs$CS_Point.File
  CS.exe<-CS.inputs$CS.exe
  EXPORT.dir<-out
  
  ######
  r <- GA.params$Resistance.stack
  
  for(i in 1:GA.params$n.layers){
    if(GA.params$surface.type[i]=="cat"){
      parm <- PARM[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]
      df <- data.frame(id=unique.rast(r[[i]]),parm) # Data frame with original raster values and replacement values
      r[[i]] <-subs(r[[i]],df)
      
      r[[i]]<-r[[i]]-(cellStats(x=r[[i]],stat="min"))   

      
    } else {
      rast <-SCALE(data=r[[i]],MIN=0,MAX=10)
      parm <- PARM[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]
      
      
      # Set equation for continuous surface
      equation <- floor(parm[1]) # Parameter can range from 1-9.99
      
      # Read in resistance surface to be optimized
      SHAPE <-  (parm[2])
      Max.SCALE <- (parm[3])
      
      # Apply specified transformation
      if(equation==1){
        SIGN=-1 # Inverse
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE)) # Monomolecular
        R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- rast.R
        EQ <- "Inverse-Reverse Monomolecular"
        
      } else if(equation==5){
        SIGN=1
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE)) # Monomolecular
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- rast.R
        EQ <- "Reverse Monomolecular"        
        
      } else if(equation==3){
        SIGN=1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE)) # Monomolecular    
        EQ <- "Monomolecular"
        
      } else if (equation==7) {
        SIGN=-1 #Inverse
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE)) # Monomolecular
        r[[i]] <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
        EQ <- "Inverse Monomolecular"        
        
      } else if (equation==8) {
        SIGN=-1 #Inverse
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE)) # Ricker
        r[[i]] <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
        EQ <- "Inverse Ricker"  
        
      } else if (equation==8) {
        SIGN=1
        r[[i]] <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE)) #  Ricker
        EQ <- "Ricker"
        
      } else if (equation==6) {
        SIGN=1
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE)) #  Ricker
        R.vec <- rev(R)
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- rast.R
        EQ <- "Reverse Ricker"        
        
      } else if (equation==2) {
        SIGN=-1 # Inverse
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE)) # Ricker
        R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- rast.R
        EQ <- "Inverse-Reverse Ricker"
        
      } else {
        r[[i]] <- (rast*0) #  Cancel layer...set to zero
      } # End if-else  
    } # Close parameter type if-else  
  } # Close layer loop
  
  
  File.name <- paste(GA.inputs$parm.type$name,collapse=".")
  
  multi_surface <- sum(r)+1 # Add all surfaces together (+1 for distance)
  if(cellStats(multi_surface,"max")>5e4)  multi_surface<-SCALE(multi_surface,1,5e4) # Rescale surface in case resistance are too high
  #       plot(multi_surface)
  
  if(is.null(out)){
    (multi_surface)
  } else {
  writeRaster(x=multi_surface,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  (multi_surface)
  }
}
###################################################################################

###############################################################################  
############ TRANSFORM RESISTANCE SURFACES ############ 
###############################################################################  
#' Apply transformation to continuous resistance surface
#' 
#' Apply on the six resistance transformations to a continuous resistance surface
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
#'    \tab 4 = "Ricker"\cr"
#'    \tab 5 = "Reverse Monomolecular"\cr
#'    \tab 6 = "Reverse Ricker"\cr
#'    \tab 7 = "Inverse Monomolecular"\cr
#'    \tab 8 = "Inverse Ricker"\cr
#'    \tab 9 = "Distance"\cr
#'    }
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
  
  ######
  
      r <-SCALE(data=R,MIN=0,MAX=10)

  # Set equation for continuous surface
      equation <- floor(parm[1]) # Parameter can range from 1-9.99
      
      # Read in resistance surface to be optimized
      SHAPE <-  (parm[2])
      Max.SCALE <- (parm[3])
      
# Apply specified transformation
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
  File.name <- NAME
  
  if(cellStats(r,"max")>5e5)  r<-SCALE(r,1,5e5) # Rescale surface in case resistance are too high

  if(!is.null(out)){
    writeRaster(x=r,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
    return(r)
  } else {
    return(r)
  }
}
###################################################################################

###############################################################################  
############ MULTISURFACE OPTIMIZATION FUNCTION WITH GA--NO GAUSSIAN ############ 
###############################################################################  
#' Optimize multiple resistance surfaces simultaneously
#' 
#' Create composite resistance surface by simultaneously optimizing multiple categoricla and continuous surfaces. This optimization function is designed to be called from GA
#' 
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Should be a vector with parameters specified in the order of resistance surfaces. These values are selected during optimization if called within GA function.
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max')
#' @param quiet Logical, if TRUE, AICc and iteration time will not be printed to the screen at the completion of each iteration. Default = FALSE
#' @return AIC value from mixed effect model
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Resistance.Opt_multi <- function(PARM,CS.inputs,GA.inputs, Min.Max, quiet=FALSE){
  t1<-Sys.time()
  
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  response<-CS.inputs$response
  #   CS.version<-CS.inputs$CS.version
  CS_Point.File<-CS.inputs$CS_Point.File
  CS.exe<-CS.inputs$CS.exe
  EXPORT.dir<-GA.inputs$Write.dir
  GA.params<-GA.inputs
  ######
  r <- GA.params$Resistance.stack
  
  for(i in 1:GA.params$n.layers){
    if(GA.params$surface.type[i]=="cat"){
      parm <- PARM[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]
      df <- data.frame(id=unique.rast(r[[i]]),parm) # Data frame with original raster values and replacement values
      r[[i]] <-subs(r[[i]],df)
      
      r[[i]]<-r[[i]]-(cellStats(x=r[[i]],stat="min"))
      
      
      #   	  cat(GA.params$layer.names[i],"\n")
      
    } else {
      rast <-SCALE(data=r[[i]],MIN=0,MAX=10)
      parm <- PARM[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]
      
      
      # Set equation for continuous surface
      equation <- floor(parm[1]) # Parameter can range from 1-9.99
      
      # Read in resistance surface to be optimized
      SHAPE <-  (parm[2])
      Max.SCALE <- (parm[3])
      
      rick.eq<-(equation==2||equation==4||equation==6||equation==8)
      if(rick.eq==TRUE & SHAPE>6){
        equation<-9
      }
      
      # Apply specified transformation
      if(equation==1){
        SIGN=-1 # Inverse
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE)) # Monomolecular
        R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Inverse-Reverse Monomolecular"
        
      } else if(equation==5){
        SIGN=1
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE)) # Monomolecular
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Reverse Monomolecular"        
        
      } else if(equation==3){
        SIGN=1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE)) # Monomolecular    
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Monomolecular"
        
      } else if (equation==7) {
        SIGN=-1 #Inverse
        R <- SIGN*Max.SCALE*(1-exp(-1*rast/SHAPE)) # Monomolecular
        r[[i]] <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min')))# Rescale
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Inverse Monomolecular"        
        
      } else if (equation==8) {
        SIGN=-1 #Inverse
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE)) # Ricker
        r[[i]] <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Inverse Ricker"  
        
      } else if (equation==4) {
        SIGN=1
        r[[i]] <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE)) #  Ricker
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Ricker"
        
      } else if (equation==6) {
        SIGN=1
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE)) #  Ricker
        R.vec <- rev(R)
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Reverse Ricker"        
        
      } else if (equation==2) {
        SIGN=-1 # Inverse
        R <- SIGN*(Max.SCALE*rast*exp(-1*rast/SHAPE)) # Ricker
        R <- SCALE(R,MIN=abs(cellStats(R,stat='max')),MAX=abs(cellStats(R,stat='min'))) # Rescale
        R.vec <- rev(R) # Reverse
        rast.R <- setValues(R,values=R.vec)
        r[[i]] <- reclassify(rast.R, c(-Inf,1e-05, 1e-05,1e6,Inf,1e6))
        EQ <- "Inverse-Reverse Ricker"
        
      } else {
        r[[i]] <- (rast*0) #  Cancel layer...set to zero
      } # End if-else  
    } # Close parameter type if-else  
  } # Close layer loop
  
  
  File.name <- "multi_surface"
  
  multi_surface <- sum(r)+1 # Add all surfaces together (+1 for distance)
  if(cellStats(multi_surface,"max")>5e5)  multi_surface<-SCALE(multi_surface,1,5e5) # Rescale surface in case resistance are too high
  
  writeRaster(x=multi_surface,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  
  # Modify and write Circuitscape.ini file
  #############################################################################################  
  BATCH<-paste0(EXPORT.dir,File.name,".ini")        
  OUT<-paste0(paste0("output_file = ",EXPORT.dir), File.name,".out")
  HABITAT<-paste0("habitat_file = ",paste0(EXPORT.dir,File.name,".asc"))
  LOCATION.FILE <- paste0("point_file = ", CS.inputs$CS_Point.File)
  ifelse(CS.inputs$Neighbor.Connect==4,connect<-"True",connect<-"False")
  CONNECTION=paste0("connect_four_neighbors_only=",connect)
  
  #   if(CS.version=='3.5.8'){
  #     write.CS_3.5.8(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,VERSION=VERSION)
  #   } else {
  write.CS_4.0(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION)    
  #   }
  
  ##########################################################################################
  # Run Circuitscape
  CS.exe<-CS.exe
  
  # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
  hidden = TRUE
  
  CS.ini <- paste0(EXPORT.dir,File.name,".ini")
  CS.Run.output<-system(paste(CS.exe, CS.ini), hidden)
  
  #########################################
  # Run mixed effect model on each Circuitscape effective resistance
  
  CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")
  
  # Get AIC statistic for transformed-scaled resistance surface
  cs.matrix<-read.matrix(CS.results)
  cs.matrix<-scale(cs.matrix,center=TRUE,scale=TRUE)
  # cs.matrix2<-round(read.matrix(CS.results),4)
  
  data<-cbind(ID,cs.matrix,response)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=data$cs.matrix)
  
  # Fit model
  mod <- lFormula(response ~ LAYER + (1|pop1), data=data,REML=FALSE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  #    summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  
  k<-max(GA.params$parm.index)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
  
 
  t2 <-Sys.time()
  if(quiet==FALSE){
  cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
  cat(paste0("\t", "AICc = ",round(AICc,4)),"\n","\n")
  }
  
  
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
}

################################################################### 
############ ITERATIVE SINGLE OPTIMIZATION FUNCTION WITH GA---No Gaussian ############ 
################################################################### 
#' Optimize resistance surfaces individually
#' 
#' Optimize all resistance surfaces that are located in the same directory individually. This optimization function is designed to be called from GA
#' 
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Should be a vector with parameters specified in the order of resistance surfaces.These values are selected during optimization if called within GA function.
#' @param Resistance Resistance surface to be optimized. This should be an R raster object. If not specified, the function will attempt to find the a resistance surface from \code{GA.inputs}
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max'). Default in 'max'
#' @param iter A counter for the number of surfaces that will be optimized
#' @param quiet Logical, if TRUE AICc and iteration duration will not be printed to the screen at the completion of each iteration.
#' @return AIC value from mixed effect model
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Resistance.Opt_single <- function(PARM,Resistance,CS.inputs,GA.inputs, Min.Max='max',iter, quiet=FALSE){
  t1<-Sys.time()
  
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  response<-CS.inputs$response
  #   CS.version<-CS.inputs$CS.version
  CS_Point.File<-CS.inputs$CS_Point.File
  CS.exe<-CS.inputs$CS.exe
  EXPORT.dir<-GA.inputs$Write.dir
  GA.params<-GA.inputs
  ######
  r <- Resistance
  
  if(GA.params$surface.type[iter]=="cat"){
    #       parm <- PARM[(GA.params$parm.index[iter]+1):(GA.params$parm.index[iter+1])]
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
  } # Close parameter type if-else  
  
  File.name <- "resist_surface"
  
  if(cellStats(r,"max")>1e6)  r<-SCALE(r,1,1e6) # Rescale surface in case resistance are too high
  r <- reclassify(r, c(-Inf,1e-06, 1e-06,1e6,Inf,1e6))
  
  writeRaster(x=r,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  
  # Modify and write Circuitscape.ini file
  #############################################################################################  
  BATCH<-paste0(EXPORT.dir,File.name,".ini")        
  OUT<-paste0(paste0("output_file = ",EXPORT.dir), File.name,".out")
  HABITAT<-paste0("habitat_file = ",paste0(EXPORT.dir,File.name,".asc"))
  LOCATION.FILE <- paste0("point_file = ", CS.inputs$CS_Point.File)
  ifelse(CS.inputs$Neighbor.Connect==4,connect<-"True",connect<-"False")
  CONNECTION=paste0("connect_four_neighbors_only=",connect)
  
  # if(CS.version=='3.5.8'){
  #   write.CS_3.5.8(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,VERSION=VERSION)
  # } else {
  write.CS_4.0(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION)    
  # }
  ##########################################################################################
  # Run Circuitscape
  CS.exe<-CS.exe
  
  # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
  hidden = TRUE
  
  # CS.Batch<- list.files(path=EXPORT.dir, pattern = "\\.ini$") # Make list of all files with '.asc' extension
  CS.ini <- paste0(EXPORT.dir,File.name,".ini")
  CS.Run.output<-system(paste(CS.exe, CS.ini), hidden) 
  
  #########################################
  # Run mixed effect model on each Circuitscape effective resistance
  
  CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")
  
  # Get AIC statistic for transformed-scaled resistance surface
  cs.matrix<-scale(read.matrix(CS.results),center=TRUE,scale=TRUE)
  
  data<-cbind(ID,cs.matrix,response)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=data$cs.matrix)
  
  # Fit model
  mod <- lFormula(response ~ LAYER + (1|pop1), data=data,REML=FALSE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  
  k<-length(PARM)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))

  
  t2 <-Sys.time()
  if(quiet==FALSE){
    cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
#     cat(paste0("\t", EQ,"; ",round(SHAPE,digits=2),"; ", round(Max.SCALE,digits=2)),"\n")
    cat(paste0("\t", "AICc = ",round(AICc,4)),"\n","\n")
    
  }
  
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
}

################################################### 
############ PLOT response CURVES ############ 
################################################### 
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
    Trans.vec <- rev(SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec)))) # Reverse
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


############################################################################  
############ OPTIMIZATION FUNCTION USING GA STARTS: CONTINUOUS ############ 
############################################################################  
Resistance.Optimization_cont.nlm<-function(PARM,Resistance,equation, get.best,CS.inputs,Min.Max,quiet=FALSE) {
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  response<-CS.inputs$response
  #   CS.version<-CS.inputs$CS.version
  CS_Point.File<-CS.inputs$CS_Point.File
  CS.exe<-CS.inputs$CS.exe
  
  EXPORT.dir<-  if(get.best==FALSE){
    EXPORT.dir<-GA.inputs$Write.dir} else{
      EXPORT.dir<-GA.inputs$Results.dir
    }
  
  
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
  
  writeRaster(x=r,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  
  #Modify and write Circuitscape.ini file
  #############################################################################################  
  BATCH<-paste0(EXPORT.dir,File.name,".ini")        
  OUT<-paste0(paste0("output_file = ",EXPORT.dir), File.name,".out")
  HABITAT<-paste0("habitat_file = ",paste0(EXPORT.dir,File.name,".asc"))
  LOCATION.FILE <- paste0("point_file = ", CS.inputs$CS_Point.File)
  ifelse(CS.inputs$Neighbor.Connect==4,connect<-"True",connect<-"False")
  CONNECTION=paste0("connect_four_neighbors_only=",connect)
  
  # if(CS.version=='3.5.8'){
  #   write.CS_3.5.8(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,VERSION=VERSION)
  # } else {
  write.CS_4.0(BATCH=BATCH,OUT=OUT,HABITAT=HABITAT,LOCATION.FILE=LOCATION.FILE,CONNECTION=CONNECTION)    
  # }
  ##########################################################################################
  # Run Circuitscape
  CS.exe<-CS.exe
  
  # Keep status of each run hidden? Set to either 'TRUE' or 'FALSE'; If 'FALSE' updates will be visible on screen
  hidden = TRUE
  
  # CS.Batch<- list.files(path=EXPORT.dir, pattern = "\\.ini$") # Make list of all files with '.asc' extension
  CS.ini <- paste0(EXPORT.dir,File.name,".ini")
  CS.Run.output<-system(paste(CS.exe, CS.ini), hidden) 
  
  ##########################################################################################
  # Run mixed effect model on each Circuitscape effective resistance
  
  CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")
  
  # Get AIC statistic for transformed-scaled resistance surface
  cs.matrix<-scale(read.matrix(CS.results),center=TRUE,scale=TRUE)
  
  data<-cbind(ID,cs.matrix,response)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=data$cs.matrix)
  
  # Fit model
  mod <- lFormula(response ~ LAYER + (1|pop1), data=data,REML=FALSE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  
  k<-length(PARM)+2
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
    
  t2 <-Sys.time()
  if(quiet==FALSE){    
  cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
  cat(paste0("\t", "AICc = ",round(AICc,3)),"\n")
  }
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized    
  
}
############################################

# Run Mixed effects models, recovery parameter estimates
#' Obtain coefficients from maximum likelihood population effects mixed effects model (MLPE)
#' 
#' Runs MLPE as detailed by Clarke et al. (2002). This function is designed to generate a table of the fitted mixed effects model coefficients
#' 
#' @param resist.dir Directory containing optimized .asc resistance files (optimized resistance surfaces)
#' @param genetic.dist Lower half of pairwise genetic distance matrix
#' @param out.dir If specified, a .csv table will printed to the specified directory (Default = NULL)
#' @return A table of MLPE fitted model coefficients

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage MLPE.lmm_coef(resist.dir, genetic.dist, out.dir)
#' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

MLPE.lmm_coef <- function(resist.dir, genetic.dist,out.dir=NULL){ 
  response=genetic.dist
  resist.mat<-list.files(resist.dir,pattern="*_resistances.out",full.names=TRUE)
  resist.names<-gsub(pattern="_resistances.out","",x=list.files(resist.dir,pattern="*_resistances.out"))
  COEF.Table<-array()
  for(i in 1:length(resist.mat)){
    m<-length(read.table(resist.mat[i])[-1,-1])
    mm<-read.table(resist.mat[i])[-1,-1]
    ID<-To.From.ID(POPS=m)
    ZZ<-ZZ.mat(ID=ID)
    cs.matrix<-scale(mm[lower.tri(mm)],center=TRUE,scale=TRUE)
    
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
  if(is.null(out.dir)){
    COEF.Table<-(COEF.Table[-1,])
  } else {
    COEF.Table<-COEF.Table[-1,]
    write.table(COEF.Table,file=paste0(out.dir,"MLPE_coeff_Table.csv"),sep=",",row.names=T,col.names=NA)
    return(COEF.Table)
  }
}



# Run Mixed effects models, recovery parameter estimates
#' Run maximum likelihood population effects mixed effects model (MLPE)
#' 
#' Runs MLPE as detailed by Clarke et al. (2002). This function will run the model and return glmer object
#' 
#' @param cs.resistance Pairwise resistance distance matrix (resistances.out) from CS results
#' @param pairwise.genetic Lower half of pairwise genetic distance matrix
#' @param REML Logical. If TRUE, mixed effects model will be fit using restricted maximum likelihood. Default = FALSE
#' @return A glmer object from the fitted model
#' @details An AIC value will only be returned if \code{REML = FALSE}

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage MLPE.lmm(cs.resistance, pairwise.genetic, REML)
#' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

MLPE.lmm <- function(cs.resistance, pairwise.genetic, REML=FALSE){ 
  response=pairwise.genetic

    mm<-(read.table(cs.resistance)[-1,-1])
    m<-nrow(mm)
    ID<-To.From.ID(POPS=m)
    ZZ<-ZZ.mat(ID=ID)
    cs.matrix<-scale(lower(mm),center=TRUE,scale=TRUE)
    
    dat<-cbind(ID,cs.matrix,response)
    
    # Assign value to layer
    LAYER<-assign("Resist",value=dat$cs.matrix)
    
    # Fit model
    mod <- lFormula(response ~ LAYER + (1|pop1), data=dat,REML=REML)
    mod$reTrms$Zt <- ZZ
    dfun <- do.call(mkLmerDevfun,mod)
    opt <- optimizeLmer(dfun)
    MOD <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))   
  return(MOD)
}

##################################
#' Create diagnostic plots 
#' 
#' This function will generate mixed effect model diagnostic plots following optimization
#' 
#' @param cs.resistance.mat Path to CIRCUITSCAPE "_resistances.out" file
#' @param genetic.dist Vector of pairwise genetic distances (lower half of pairwise matrix). Can be input as CS.inputs$response
#' @param XLAB Label for x-axis (Defaults to "Estimated resistance")
#' @param YLAB Label for y-axis (Defaults to "Genetic distance")
#' @param plot.dir Directory to output PDF of diagnostic plots
#' @return A four panel PDF including residual scatterplot, historgram of residuals, qqplot, and 

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage Diagnostic.Plots(cs.resistance.mat, genetic.dist, XLAB,YLAB, plot.dir)

Diagnostic.Plots<-function(cs.resistance.mat, genetic.dist, XLAB="Estimated resistance",YLAB ="Genetic distance",plot.dir){
  response=genetic.dist
  NAME<-gsub(pattern="*_resistances.out","",x=(basename(cs.resistance.mat)))
  mm<-read.table(cs.resistance.mat)[-1,-1]
  m<-length(mm)
  ID<-To.From.ID(POPS=m)
  ZZ<-ZZ.mat(ID=ID)
  cs.matrix<-scale(lower(mm),center=TRUE,scale=TRUE)
  cs.unscale<-lower(mm)
  dat<-cbind(ID,cs.matrix,response)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=dat$cs.matrix)
  
  # Fit model
  mod <- lFormula(response ~ LAYER + (1|pop1), data=dat,REML=TRUE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  Mod <- (mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  #######
  # Make diagnostic plots
  #   par(mfrow=c(2,2))
  tiff(filename = paste0(plot.dir,NAME,"_DiagnosticPlots.tif"), 
       width = 279, height = 215, units = "mm", 
       compression = c("lzw"),
       bg = "white", res = 300)
  par(mfrow=c(2,2),
      oma = c(0,4,0,0) + 0.1,
      mar = c(4,4,1,1) + 0.1)
  plot(genetic.dist~cs.unscale,xlab=XLAB,ylab=YLAB)
  abline(lm(genetic.dist~cs.unscale))
  plot(residuals(Mod)~cs.unscale,xlab=XLAB,ylab="Residuals")
  abline(lm(residuals(Mod)~cs.unscale))
  hist(residuals(Mod),xlab="Residuals",main="")
  qqnorm(resid(Mod),main="")
  qqline(resid(Mod))
  dev.off()
  par(mfrow=c(1,1))  
}
##########################################################################################
# Function to bundle input parameters

#' Prepare and bundle input CIRCUITSCAPE model parameters
#' 
#' This function will prepare objects needed for running optimization functions
#' 
#' @param n.POPS The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param CS_Point.File The path to the Circuitscape formatted point file. See Circuitscape documentation for help.
#' @param CS.exe The path to the CIRCUITSCAPE executable file (cs_run.exe). See details below. 
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage CS.prep(n.POPS, response, CS_Point.File, CS.exe, Neighbor.Connect)
#' @details \code{CS.exe} Example of path to CIRCUITSCAPE executible: 
#' 
#' '"C:/Program Files/Circuitscape/cs_run.exe"'
#'
#' ***NOTE: Double quotation used***
CS.prep <- function(n.POPS, response=NULL,CS_Point.File,CS.exe,Neighbor.Connect=8){# Make to-from population list
  ID<-To.From.ID(n.POPS)
  ZZ<-ZZ.mat(ID)
  list(ID=ID,ZZ=ZZ,response=response,CS_Point.File=CS_Point.File,CS.exe=CS.exe,Neighbor.Connect=Neighbor.Connect,n.POPS=n.POPS)
}

##########################################################################################
# Data processing for analyzing multiple layers simultaneously

#' Create R object with genetic algorithm optimization settings
#' 
#' This function prepares and compiles objects and commands for optimization with the GA package
#' 
#' @param ASCII.dir Directory containing all raster objects to optimized.
#' @param min.cat The minimum value to be assessed during optimization of of categorical resistance surfaces (Default = 0)
#' @param max.cat The maximum value to be assessed during optimization of of categorical resistance surfaces (Default = 2500)
#' @param max.cont The maximum value to be assessed during optimization of of continuous resistance surfaces (Default = 2500)
#' @param cont.shape A vector of hypothesized relationships that each continuous resistance surface will have in relation to the genetic distance reposnse (Default = NULL; see details)
#' @param pop.mult Value will be multiplied with number of parameters in surface to determine 'popSize' in GA. By default this is set to 15.
#' @param percent.elite Percent used to determine the number of best fitness individuals to survive at each generation ('elitism' in GA). By default the top 5\% individuals will survive at each iteration.
#' @param type Default is "real-valued"
#' @param population Default is gareal_Population from GA
#' @param selection Default is gareal_lsSelection from GA
#' @param mutation Default is gareal_raMutation from GA
#' @param pcrossover Probability of crossover. Default = 0.85
#' @param pmutation Probability of mutation. Default = 0.1
#' @param crossover Default = "gareal_blxCrossover". This crossover method greatly improved optimization during preliminary testing
#' @param maxiter Maximum number of iterations to run before the GA search is halted (Default = 1000)
#' @param run Number of consecutive generations without any improvement in AICc before the GA is stopped (Default = 25)
#' @param keepBest A logical argument specifying if best solutions at each iteration should be saved (Default = TRUE)
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max' = Default). Optimization with \code{ga} maximizes the objective criteria
#' @param seed Integer random number seed to replicate \code{ga} optimization
#' @param quiet Logical. If TRUE, AICc and iteration time will not be printed to the screen after each iteration. Default = FALSE
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
#' Min.Max="max",
#' min.cat=0,
#' max.cat=2500,
#' max.cont=2500,
#' cont.shape=NULL,
#' pop.mult = 15,
#' percent.elite = 0.05,
#' type= "real-valued",
#' pcrossover=0.85,
#' pmutation=0.1,
#' maxiter=1000,
#' run=25,
#' keepBest=TRUE,
#' population = gaControl(type)$population,
#' selection = gaControl(type)$selection,
#' crossover="gareal_blxCrossover",
#' mutation = gaControl(type)$mutation,
#' seed = NULL,
#' quiet = FALSE)

GA.prep<-function(ASCII.dir,
                  Min.Max ='max',
                  min.cat = 0,
                  max.cat = 2500, 
                  max.cont = 2500,
                  cont.shape = NULL,
                  pop.mult = 15,
                  percent.elite = 0.05,
                  type = "real-valued",
                  pcrossover = 0.85,
                  pmutation = 0.1,
                  maxiter = 1000,
                  run = 25,
                  keepBest = TRUE,
                  population = gaControl(type)$population,
                  selection = gaControl(type)$selection,
                  crossover="gareal_blxCrossover",
                  mutation = gaControl(type)$mutation,
                  seed = NULL,
                  quiet = FALSE) {   
  
  ASCII.list <-list.files(ASCII.dir,pattern="*.asc", full.names=TRUE) # Get all .asc files from directory
  
  r <- stack(lapply(ASCII.list,raster))
  
  names <- gsub(pattern="*.asc","",x=(list.files(ASCII.dir,pattern="*.asc")))
  n.layers <-length(ASCII.list) 
  
  if("Results"%in%dir(ASCII.dir)==FALSE) dir.create(file.path(ASCII.dir, "Results")) 
  #   dir.create(file.path(ASCII.dir, "Results"),showWarnings = FALSE)
  Results.dir<-paste0(ASCII.dir, "Results/")
  if("tmp"%in%dir(ASCII.dir)==FALSE) dir.create(file.path(ASCII.dir, "tmp")) 
  #   dir.create(file.path(Results.dir, "tmp"),showWarnings = FALSE)
  Write.dir <-paste0(ASCII.dir,"tmp/") 
  if("Plots"%in%dir(Results.dir)==FALSE) dir.create(file.path(Results.dir, "Plots")) 
  #   dir.create(file.path(Results.dir, "tmp"),showWarnings = FALSE)
  Plots.dir <-paste0(Results.dir,"Plots/") 
  
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
  
  if (length(ga.min)<10){
    pop.size <- min(c(15*length(ga.min),100))
  } else {
    pop.size <- 10*length(ga.min)
  }
  
  for(i in 1:length(surface.type)){
    if (surface.type[i]=="cat"){
      SUGGESTS[[i]] <- sv.cat(levels=parm.type[i,2],pop.size=pop.size)
      
    } else if (exists("cont.shape") && length(cont.shape>0)){
      SUGGESTS[[i]] <- sv.cont.nG(cont.shape[1],pop.size=pop.size)
      cont.shape<-cont.shape[-1]
    } else {
      SUGGESTS[[i]] <- sv.cont.nG("NA",pop.size=pop.size)
    }
  }
  SUGGESTS <-matrix(unlist(SUGGESTS), nrow=nrow(SUGGESTS[[1]]), byrow=F)
  
  list(parm.index=parm.index,ga.min=ga.min,ga.max=ga.max,surface.type=surface.type,parm.type=parm.type,Resistance.stack=r,n.layers=n.layers,layer.names=names,pop.size=pop.size, min.list=min.list,max.list=max.list, SUGGESTS=SUGGESTS,ASCII.dir=ASCII.dir, Results.dir=Results.dir, Write.dir=Write.dir,Plots.dir=Plots.dir,type= type, pcrossover=pcrossover, pmutation=pmutation, crossover=crossover, maxiter=maxiter, run=run, keepBest=keepBest, population=population,selection=selection,mutation=mutation,pop.mult = pop.mult, percent.elite = percent.elite,Min.Max=Min.Max, seed=seed, quiet = quiet)  
  
}
#####################################
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
###########################

##############################################################
############ OTHER NECESSARY FUNCTIONS  #####################
#############################################################

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

read.matrix<-function(cs.matrix){  m<-read.table(cs.matrix)[-1,-1]
                                   m[lower.tri(m)]}

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

# Rescale function
SCALE.vector <-function(data,MIN,MAX){Mn=min(data)
                                      Mx=max(data)
                                      (MAX-MIN)/(Mx-Mn)*(data-Mx)+MAX}

# Define scaling function
# This will rescale from 1 to specified MAX
SCALE <-function(data,MIN,MAX){Mn=cellStats(data,stat='min')
                               Mx=cellStats(data,stat='max')
                               (MAX-MIN)/(Mx-Mn)*(data-Mx)+MAX}
# Function to write .ini file for Circuitscape 
write.CS_3.5.8 <- function(BATCH,OUT,HABITAT,LOCATION.FILE,CONNECTION,MAP="write_cum_cur_map_only=False"){
sink(BATCH)
cat("[Options for advanced mode]
ground_file_is_resistances=True
source_file=(Browseforacurrentsourcefile)
remove_src_or_gnd=keepall
ground_file=(Browseforagroundpointfile)
use_unit_currents=False
use_direct_grounds=False

[Calculation options]
low_memory_mode=False
solver=cg+amg
print_timings=True

[Options for pairwise and one-to-all and all-to-one modes]
included_pairs_file=None
point_file_contains_polygons=False
use_included_pairs=False")
cat("\n")
cat(LOCATION.FILE)
cat("\n")

cat("
[Output options]")
cat("\n")
cat(MAP)
cat("\n")
cat("log_transform_maps=False
set_focal_node_currents_to_zero=False
write_max_cur_maps=False
write_volt_maps=False
set_null_currents_to_nodata=True
set_null_voltages_to_nodata=True
compress_grids=False
write_cur_maps=False")
cat("\n")
cat(OUT)
cat("\n")
cat("\n")
cat("[Shortcircuit regions(aka polygons)]
use_polygons=False
polygon_file=(Browse for a short-circuit region file)

[Connection scheme for raster habitat data]")
cat(CONNECTION)
cat("\n")
cat("connect_using_avg_resistances=True

[Habitat raster or graph]
habitat_map_is_resistances=True")
cat("\n")
cat(HABITAT)
cat("\n")
cat("\n")
cat("[Options for one-to-all and all-to-one modes]
use_variable_source_strengths=False
variable_source_file=None

[Version]
version = 3.5.8

[Maskfile]
use_mask=False
mask_file=None

[Circuitscape mode]
data_type=raster
scenario=pairwise")
sink()
}


write.CS_4.0 <- function(BATCH,OUT,HABITAT,LOCATION.FILE,CONNECTION,CURRENT.MAP="write_cur_maps = False", MAP="write_cum_cur_map_only = False",PARALLELIZE="parallelize = False",CORES="max_parallel = 0"){
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
version = 4.0-beta

[Options for reclassification of habitat data]
reclass_file = (Browse for file with reclassification data)
use_reclass_table = False

[Logging Options]
log_level = INFO
log_file = None
profiler_log_file = None
screenprint_log = False

[Options for pairwise and one-to-all and all-to-one modes]
included_pairs_file = (Browse for a file with pairs to include or exclude)
use_included_pairs = False")
cat("\n")
cat(LOCATION.FILE)
cat("\n")
cat("\n")

cat("
[Connection scheme for raster habitat data]
connect_using_avg_resistances = True")
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

#############################################################
# Sample values for suggests
sv.cat <-function(levels,pop.size){
  cat.starts<-matrix(nrow=pop.size,ncol=levels)
  for(r in 1:pop.size){
    L<-list()
    for (i in 1:levels){
      if(runif(1)<.5){
        z<-runif(1)
      } else {
        z<-runif(1,10,250)
      }
      L[[i]]<-z
    } 
    #   uz<-unlist(L)
    cat.starts[r,]<-(unlist(L))    
  }
  cat.starts[,1]<-1
  return(cat.starts)
}
############################
# No Gaussian distribution
sv.cont.nG <-function(direction,pop.size){
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
      z<-c(sample(dec,1),runif(1,.01,10),runif(1,1,100))
    } else if (runif(1)<.5 && direction=="Peaked") {
      z<-c(sample(peak,1),runif(1,.01,10),runif(1,1,100))
    } else {
      z<-c(runif(1,1,9.99),runif(1,.01,10),runif(1,1,100))
    }
    cont.starts[r,]<-z
  } 
  cont.starts
}

############################
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
###########################
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

Result.txt <- function(GA.results, GA.inputs, CS.inputs){
  summary.file<-paste0(GA.inputs$Results.dir,"Multisurface_Optim_Summary.txt")
  AICc<-GA.results@fitnessValue
  AICc<-round(AICc,digits=4)
  ELITE<-floor(GA.inputs$percent.elite*GA.inputs$pop.size)
#   mlpe.results<-MLPE.lmm_coef(GA.inputs$Results.dir,genetic.dist=CS.inputs$response)
  
sink(summary.file)
cat(paste0("Summary from multisurface optimization run conducted on ",Sys.Date()),"\n")
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
cat(paste0("Minimum AICc: ",AICc),"\n")
cat("\n")
cat(paste0("Optimized values for each surface:"),"\n")
cat(GA.results@solution,"\n")
cat("\n")
sink()
}
# Optimiazation preparation
Optim.input<-function(Response,n.Pops,ASCII.dir,CS_Point.File,CS.exe,Neighbor.Connect=8,Constrained.Max=100,Initial.shape=c(seq(0.2,1,by=0.2),seq(1.25,10.75,by=0.75)),Bootstrap=FALSE,boot.iters=10000,Sample_Proportion=0.75){
  # Install necessary packages
  libs=c("raster", "lme4", "plyr")
  CheckInstallPackage(packages=libs)
  
  #####################
  if(is.vector(Response)==TRUE || dim(Response)[1]!=dim(Response)[2]) {warning("Must provide square distance matrix with no column or row names")}
  Response.vec<-Response[lower.tri(Response)]
  ID <- To.From.ID(n.Pops)
  ZZ<-ZZ.mat(ID)
  ASCII.files<-list.files(ASCII.dir,pattern="*.asc",full.names=TRUE)
  ASCII.names<-gsub(pattern="*.asc","",x=(list.files(ASCII.dir,pattern="*.asc")))
  dir.create(file.path(ASCII.dir, "Results"))
  Results.dir<-paste0(ASCII.dir, "/Results/")
  dir.create(file.path(Results.dir, "tmp"))
  Write.dir <-paste0(Results.dir,"/tmp/")  
  
  list(Response.vec=Response.vec,Response.mat=Response,n.Pops=n.Pops,ID=ID,ZZ=ZZ,ASCII.files=ASCII.files,ASCII.names=ASCII.names,ASCII.dir=ASCII.dir,Write.dir=Write.dir,Results.dir=Results.dir,CS_Point.File=CS_Point.File,CS.exe=CS.exe,Neighbor.Connect=Neighbor.Connect,Constrained.Max=Constrained.Max,Initial.shape=Initial.shape,Bootstrap=Bootstrap,boot.iters=boot.iters,Sample_Proportion=Sample_Proportion)
}

#########################################
# Install necessary packages
CheckInstallPackage <- function(packages, repos="http://cran.r-project.org") {
  installed=as.data.frame(installed.packages())
  for(p in packages) {
    if(is.na(charmatch(p, installed[,1]))) { 
      install.packages(p, repos=repos) 
    }
  }
} 
##########################################################################################

get.name<-function(x){
   nm <-deparse(substitute(x))
   return(nm)
   }
