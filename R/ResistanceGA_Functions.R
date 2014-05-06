########################################################################################  
############ Single command function to execute single surface optimization ############ 
########################################################################################  
#' Single surface optimization
#' 
#' Optimize all surfaces contained in a directory using a genetic algorithm executed with the \code{\link[GA]{ga}} function in the Genetic Algorithms package \pkg{GA}
#' 
#' @param CS.inputs Object created from running \code{\link{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link{GA.prep}} function
#' @return This function optimizes multiple resistance surfaces. Following optimization, several summary objects are created.\cr
#' \enumerate{
#' \item Diagnostic plots of model fit are output to the "Results/Plots" directory that is automatically generated within the folder containing the optimized ASCII files.
#' \item A .csv file with the Maximum Likelihood Population Effects mixed effects model coefficient estimates
#' \item Three summary .csv files are generated: CategoricalResults.csv, ContinuousResults.csv, & All_Results_AICc.csv. These tables contain AICc values and optimization summaries for each surface.
#' }
#' @usage SS_optim(CS.inputs, GA.inputs)

#' @export
SS_optim <- function(CS.inputs,GA.inputs){
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
      
      single.GA <-ga(type= "real-valued",fitness=Resistance.Opt_single,Resistance=r, 
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
                     iter=i) 
      
      df <- data.frame(id=unique.rast(r),single.GA@solution) 
      r <-subs(r,df)
      names(r)<-GA.inputs$layer.names[i]
      
      Run_CS(CS.inputs,GA.inputs,r,EXPORT.dir=GA.inputs$Results.dir)
      
      Diagnostic.Plots(cs.resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),genetic.dist=CS.inputs$RESPONSE,plot.dir=GA.inputs$Plots.dir)
  
   
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
      
      single.GA <-ga(type= "real-valued",fitness=Resistance.Opt_single,Resistance=r, 
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
                     iter=i)  
      
      # Using GA results, optimize with nlm  
      start.vals <- single.GA@solution[-1]
      
      # Informed start values; these are the optimized values from the single parameter optimization
      EQ <-get.EQ(single.GA@solution[1])
      Optim.nlm <-nlm(Resistance.Optimization_cont.nlm, log(start.vals), Resistance=r, equation=single.GA@solution[1],get.best=FALSE,CS.inputs=CS.inputs,Min.Max='min')
      
      OPTIM <- Resistance.Optimization_cont.nlm(PARM=(Optim.nlm$estimate),Resistance=r, equation=single.GA@solution[1],get.best=TRUE,CS.inputs=CS.inputs,Min.Max='min')
      
      Diagnostic.Plots(cs.resistance.mat=paste0(GA.inputs$Results.dir,GA.inputs$layer.names[i],"_resistances.out"),genetic.dist=CS.inputs$RESPONSE,plot.dir=GA.inputs$Plots.dir)
      
      RS<-data.frame(GA.inputs$layer.names[i],Optim.nlm$minimum,EQ,Cont.Param_nG(exp(Optim.nlm$estimate)))
      colnames(RS) <- c("Surface","AICc","Equation","shape","max")
      RESULTS.cont[[cnt2]] <- RS
      
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
  Features <- array()
  for(i in 1:ncol(Results.cat)-2){
    feature <- paste0("Feature",i)
    Features[i]<-feature
  }
  colnames(Results.cat)<-c("Surface","AICc", Features)
  colnames(Results.cont)<-c("Surface","AICc","Equation","shape","max")
  Results.cat
  Results.cont
  
  # Full Results
  Results.All<-rbind(Results.cat[,c(1,2)],Results.cont[,c(1,2)])
  
  write.table(Results.cat,paste0(GA.inputs$Results.dir,"CategoricalResults.csv"),sep=",",col.names=T,row.names=F)
  write.table(Results.cont,paste0(GA.inputs$Results.dir,"ContinuousResults.csv"),sep=",",col.names=T,row.names=F)
  write.table(Results.All,paste0(GA.inputs$Results.dir,"All_Results_AICc.csv"),sep=",",col.names=T,row.names=F)
  
  # Get parameter estimates
  MLPE.results<-MLPE.lmm(resist.dir=GA.inputs$Results.dir,genetic.dist=CS.inputs$RESPONSE,out.dir=GA.inputs$Results.dir)
  return(Results.All)
  ###############################################################################################################
}


########################################################################################  
############ Single command function to execute mulit surface optimization ############ 
########################################################################################  
#' Simultaneous optimization of multiple resistance surfaces
#' 
#' Optimize multiple resistance surfaces simultaneously using genetic algorithms
#' 
#' @param CS.inputs Object created from running \code{\link{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link{GA.prep}} function
#' @return This function optimizes multiple resistance surfaces, returning a Genetic Algorithm (GA) object with summary information. Diagnostic plots of model fit are output to the "Results/Plots" folder that is automatically generated within the folder containing the optimized ASCII files. A text summary of the optimization settings and results is printed to the results folder.
#' @usage MS_optim(CS.inputs, GA.inputs)

#' @export
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
                   suggestions=GA.inputs$SUGGESTS) 
  
  RAST<-Combine_Surfaces(multi.GA_nG@solution,CS.inputs,GA.inputs)
  NAME<-paste(GA.inputs$parm.type$name,collapse=".")
  name(RAST)<-NAME
  Run_CS(CS.inputs,GA.inputs,r=RAST,CurrentMap=FALSE)
  
  NAME<-paste(GA.inputs$parm.type$name,collapse=".")
  Diagnostic.Plots(cs.resistance.mat=paste0(GA.inputs$Results.dir,NAME,"_resistances.out"),genetic.dist=CS.inputs$RESPONSE,plot.dir=GA.inputs$Plots.dir)
  
  Result.txt(GA.results=multi.GA_nG,GA.inputs=GA.inputs) 
  return(multi.GA_nG)
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
#' @param CS.inputs Object created from running \code{\link{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link{GA.prep}} function
#' @param r Raster resistance surface
#' @param CurrentMap Logical. If TRUE, the cumulative resistance map will be generated during the CS run (Default = FALSE)
#' @param EXPORT.dir Directory where CS results should be written (Default = GA.inputs$Write.dir, which is a temporary directory for reading/writing CS results)
#' @return Vector of CIRCUITSCAPE resistance distances (lower half of "_resistances.out")
#' @usage Run_CS(CS.inputs, GA.inputs, r, CurrentMap, EXPORT.dir)

#' @export
  Run_CS <- function(CS.inputs,GA.inputs,r,CurrentMap=FALSE,EXPORT.dir=GA.inputs$Write.dir){
  #   t1<-Sys.time()
  
  if (CurrentMap==FALSE){
    File.name <- r@data@names
    MAP="write_cum_cur_map_only = False"
    CURRENT.MAP="write_cur_maps = False"
    
  } else {
    File.name <- r@data@names
    MAP="write_cum_cur_map_only = True"
    CURRENT.MAP="write_cur_maps = 1"
  }
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  RESPONSE<-CS.inputs$RESPONSE
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
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Should be a vector with parameters specified in the order of resistance surfaces
#' @param CS.inputs Object created from running \code{\link{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link{GA.prep}} function
#' @param out Directory to write combined .asc file. Default is GA.inputs$Results.dir
#' @return R raster object
#' @export

Combine_Surfaces <- function(PARM,CS.inputs,GA.inputs, out=GA.inputs$Results.dir){
  t1<-Sys.time()
  GA.params<-GA.inputs
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  RESPONSE<-CS.inputs$RESPONSE
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
      
      
      #       cat(GA.params$layer.names[i],"\n")
      
    } else {
      r[[i]] <-SCALE(data=r[[i]],MIN=0,MAX=10)
      parm <- PARM[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]
      
      
      # Set equation for continuous surface
      equation <- floor(parm[1]) # Parameter can range from 1-6.99
      
      # Read in resistance surface to be optimized
      SHAPE <-  (parm[2])
      Max.SCALE <- (parm[3])
      
      # Apply specified transformation
      if(equation==1){
        SIGN=-1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(r[[i]]/SHAPE)) # Inverse-Reverse Monomolecular
        r[[i]] <- SCALE(r[[i]],MIN=abs(cellStats(r[[i]],stat='max')),MAX=abs(cellStats(r[[i]],stat='min')))
        EQ <- "Inverse-Reverse Monomolecular"
        
      } else if(equation==2){
        SIGN=1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(r[[i]]/SHAPE)) # Reverse Monomolecular
        r[[i]] <- SCALE(r[[i]],MIN=abs(cellStats(r[[i]],stat='max')),MAX=abs(cellStats(r[[i]],stat='min')))
        EQ <- "Reverse Monomolecular"        
        
      } else if(equation==3){
        SIGN=1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(-1*r[[i]]/SHAPE)) # Monomolecular
        EQ <- "Monomolecular"
        
      } else if (equation==4) {
        SIGN=-1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(-1*r[[i]]/SHAPE)) # Inverse Monomolecular
        r[[i]] <- SCALE(r[[i]],MIN=abs(cellStats(r[[i]],stat='max')),MAX=abs(cellStats(r[[i]],stat='min')))
        EQ <- "Inverse Monomolecular"        
        
      } else if (equation==5) {
        SIGN=-1
        r[[i]] <- SIGN*(Max.SCALE*r[[i]]*exp(-1*r[[i]]/SHAPE)) # Inverse Ricker
        r[[i]] <- SCALE(r[[i]],MIN=abs(cellStats(r[[i]],stat='max')),MAX=abs(cellStats(r[[i]],stat='min')))
        EQ <- "Inverse Ricker"  
        
      } else if (equation==6) {
        SIGN=1
        r[[i]] <- SIGN*(Max.SCALE*r[[i]]*exp(-1*r[[i]]/SHAPE)) #  Ricker
        EQ <- "Ricker"     
        
      } else {
        r[[i]] <- reclassify(r[[i]], c(-Inf,Inf, 0)) # Cancel surface...set to zero  
      } # End if-else
    } # Close parameter type if-else  
  } # Close layer loop
  
  
  File.name <- paste(GA.inputs$parm.type$name,collapse=".")
  
  multi_surface <- sum(r)+1 # Add all surfaces together (+1 for distance)
  if(cellStats(multi_surface,"max")>5e4)  multi_surface<-SCALE(multi_surface,1,5e4) # Rescale surface in case resistance are too high
  #       plot(multi_surface)
  
  writeRaster(x=multi_surface,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  (multi_surface)
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
#' @param CS.inputs Object created from running \code{\link{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max')
#' @return AIC value from mixed effect model
#' @export
Resistance.Opt_multi <- function(PARM,CS.inputs,GA.inputs, Min.Max){
  t1<-Sys.time()
  
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  RESPONSE<-CS.inputs$RESPONSE
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
      r[[i]] <-SCALE(data=r[[i]],MIN=0,MAX=10)
      parm <- PARM[(GA.params$parm.index[i]+1):(GA.params$parm.index[i+1])]
      
      
      # Set equation for continuous surface
      equation <- floor(parm[1]) # Parameter can range from 1-6.99
      
      # Read in resistance surface to be optimized
      SHAPE <- (parm[2])
      Max.SCALE <- (parm[3])
      
      
      # Apply specified transformation
      if(equation==1){
        SIGN=-1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(r[[i]]/SHAPE)) # Inverse-Reverse Monomolecular
        r[[i]] <- SCALE(r[[i]],MIN=abs(cellStats(r[[i]],stat='max')),MAX=abs(cellStats(r[[i]],stat='min')))
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-06, 1e-06,10e6,Inf,10e6))
        
        EQ <- "Inverse-Reverse Monomolecular"
        
      } else if(equation==2){
        SIGN=1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(r[[i]]/SHAPE)) # Reverse Monomolecular
        r[[i]] <- SCALE(r[[i]],MIN=abs(cellStats(r[[i]],stat='max')),MAX=abs(cellStats(r[[i]],stat='min')))
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-06, 1e-06,10e6,Inf,10e6))
        EQ <- "Reverse Monomolecular"  	    
        
      } else if(equation==3){
        SIGN=1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(-1*r[[i]]/SHAPE)) # Monomolecular
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-06, 1e-06,10e6,Inf,10e6))
        EQ <- "Monomolecular"
        
      } else if (equation==4) {
        SIGN=-1
        r[[i]] <- SIGN*Max.SCALE*(1-exp(-1*r[[i]]/SHAPE)) # Inverse Monomolecular
        r[[i]] <- SCALE(r[[i]],MIN=abs(cellStats(r[[i]],stat='max')),MAX=abs(cellStats(r[[i]],stat='min')))
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-06, 1e-06,10e6,Inf,10e6))
        EQ <- "Inverse Monomolecular"  	    
        
      } else if (equation==5) {
        SIGN=-1
        r[[i]] <- SIGN*(Max.SCALE*r[[i]]*exp(-1*r[[i]]/SHAPE)) # Inverse Ricker
        r[[i]] <- SCALE(r[[i]],MIN=abs(cellStats(r[[i]],stat='max')),MAX=abs(cellStats(r[[i]],stat='min')))
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-06, 1e-06,10e6,Inf,10e6))
        EQ <- "Inverse Ricker"  
        
      } else if (equation==6) {
        SIGN=1
        r[[i]] <- SIGN*(Max.SCALE*r[[i]]*exp(-1*r[[i]]/SHAPE)) #  Ricker
        r[[i]] <- reclassify(r[[i]], c(-Inf,1e-06, 1e-06,10e6,Inf,10e6))
        EQ <- "Ricker"
        
      } else {
        r[[i]] <- reclassify(r[[i]], c(-Inf,Inf, 0)) # Cancel surface...set to zero
      } # End if-else
    } # Close parameter type if-else  
  } # Close layer loop
  
  
  File.name <- "multi_surface"
  
  multi_surface <- sum(r)+1 # Add all surfaces together (+1 for distance)
  if(cellStats(multi_surface,"max")>5e4)  multi_surface<-SCALE(multi_surface,1,5e4) # Rescale surface in case resistance are too high
  
  #       plot(multi_surface)
  
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
  # CS.Run.output<-system(paste(CS.exe, CS.ini), hidden,minimized=FALSE) 
  
  #########################################
  # Run mixed effect model on each Circuitscape effective resistance
  
  CS.results<-paste0(EXPORT.dir,File.name,"_resistances.out")
  
  # Get AIC statistic for transformed-scaled resistance surface
  cs.matrix<-read.matrix(CS.results)
  cs.matrix<-scale(cs.matrix,center=TRUE,scale=TRUE)
  # cs.matrix2<-round(read.matrix(CS.results),4)
  
  data<-cbind(ID,cs.matrix,RESPONSE)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=data$cs.matrix)
  
  # Fit model
  mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=data,REML=FALSE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  #    summary(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  
  k<-max(GA.params$parm.index)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
  
  t2 <-Sys.time()
  cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
  cat(paste0("\t", "AIC = ",round(AICc,4)),"\n","\n")
  
  
  
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
#' @param CS.inputs Object created from running \code{\link{CS.prep}} function
#' @param GA.inputs Object created from running \code{\link{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max')
#' @return AIC value from mixed effect model
#' @export
Resistance.Opt_single <- function(PARM,Resistance,CS.inputs,GA.inputs, Min.Max,iter){
  t1<-Sys.time()
  
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  RESPONSE<-CS.inputs$RESPONSE
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
    
    #      r<-r-(cellStats(x=r,stat="min"))
    
    
    #   	  cat(GA.params$layer.names[i],"\n")
    
  } else {
    #      r <-SCALE(data=r,MIN=0,MAX=10)
    #       parm <- PARM[(GA.params$parm.index[iter]+1):(GA.params$parm.index[iter+1])]
    
    
    # Set equation for continuous surface
    equation <- floor(PARM[1]) # Parameter can range from 1-6.99
    
    # Read in resistance surface to be optimized
    SHAPE <- (PARM[2])
    Max.SCALE <- (PARM[3])
    
    # Apply specified transformation
    if(equation==1){
      SIGN=-1
      r <- SIGN*Max.SCALE*(1-exp(r/SHAPE))+SIGN # Inverse-Reverse Monomolecular
      r <- SCALE(r,MIN=abs(cellStats(r,stat='max')),MAX=abs(cellStats(r,stat='min')))
      EQ <- "Inverse-Reverse Monomolecular"
      
    } else if(equation==2){
      SIGN=1
      r <- SIGN*Max.SCALE*(1-exp(r/SHAPE))+SIGN # Reverse Monomolecular
      r <- SCALE(r,MIN=abs(cellStats(r,stat='max')),MAX=abs(cellStats(r,stat='min')))
      EQ <- "Reverse Monomolecular"  	    
      
    } else if(equation==3){
      SIGN=1
      r <- SIGN*Max.SCALE*(1-exp(-1*r/SHAPE))+SIGN # Monomolecular
      EQ <- "Monomolecular"
      
    } else if (equation==4) {
      SIGN=-1
      r <- SIGN*Max.SCALE*(1-exp(-1*r/SHAPE))+SIGN # Inverse Monomolecular
      r <- SCALE(r,MIN=abs(cellStats(r,stat='max')),MAX=abs(cellStats(r,stat='min')))
      EQ <- "Inverse Monomolecular"  	    
      
    } else if (equation==5) {
      SIGN=-1
      r <- SIGN*(Max.SCALE*r*exp(-1*r/SHAPE))+SIGN # Inverse Ricker
      r <- SCALE(r,MIN=abs(cellStats(r,stat='max')),MAX=abs(cellStats(r,stat='min')))
      EQ <- "Inverse Ricker"  
      
    } else if (equation==6) {
      SIGN=1
      r <- SIGN*(Max.SCALE*r*exp(-1*r/SHAPE))+SIGN #  Ricker
      EQ <- "Ricker"
    } else {
      r <- (r*0)+1 #  Distance
      EQ <- "Distance"    
    } # End if-else
  } # Close parameter type if-else  
  
  File.name <- "resist_surface"
  
  if(cellStats(r,"max")>5e4)  r<-SCALE(r,1,5e4) # Rescale surface in case resistance are too high
  r <- reclassify(r, c(-Inf,1e-06, 1e-06,10e6,Inf,10e6))
  
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
  
  data<-cbind(ID,cs.matrix,RESPONSE)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=data$cs.matrix)
  
  # Fit model
  mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=data,REML=FALSE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  
  k<-length(PARM)+1
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
  t2 <-Sys.time()
  cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
  cat(paste0("\t", "AICc = ",round(AICc,4)),"\n","\n")
  
  
  OPTIM.DIRECTION(Min.Max)*(AICc) # Function to be minimized/maximized      
}

################################################### 
############ PLOT RESPONSE CURVES ############ 
################################################### 
#' Plot continuous surface transformation
#' 
#' Plots a transformed continuous resistance surface against the original resistance values
#' 
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. A vector of two parameters is required. The first term isthe value of shape parameter (c), and the second term is the value of maximum scale parameter (b)
#' @param Resistance Path to raw, untransformed resistance surface
#' @param equation Name of the transformation equation to use:\cr
#' "Inverse-Reverse Monomolecular"\cr"Inverse Monomolecular"\cr"Monomolecular"\cr"Reverse Monomolecular"\cr"Inverse Ricker"\cr"Ricker",\cr"Distance"
#' @param print.dir Specify the directory where a .tiff of the transformation will be written (Default = NULL)
#' @return plot of transformed resistance values against original resistance values
#' @details This function will create a ggplot object and plot, so it requires \pkg{ggplot2} to be installed.\cr Equation names can be "Inverse-Reverse Monomolecular", "Inverse Monomolecular", "Monomolecular", "Reverse Monomolecular", "Inverse Ricker", "Ricker", or "Distance". The "Distance" equation sets all cell values equal to 1.
#' @usage PLOT.trans(PARM, Resistance, equation, print.dir)
#' @export
#' @import ggplot2

PLOT.trans <- function(PARM,Resistance,equation, print.dir=NULL){
  Resistance <- raster(Resistance)
  Mn=cellStats(Resistance,stat='min')
  Mx=cellStats(Resistance,stat='max') 
  
  # Make vector of data
  dat.o <- seq(from=Mn,to=Mx,length.out=1000)
  dat.t <- SCALE.vector(data=dat.o,0,10)
  
  SHAPE <- PARM[[1]]
  Max.SCALE <- PARM[[2]]
  
  # Set equation/name combination
  if(equation=="Distance") {
    Trans.vec <- (dat.t*0)+1 
    TITLE <- "Distance"
  } else if(equation=="Inverse-Reverse Monomolecular"){
    SIGN<- -1
    Trans.vec <- SCALE.vector(SIGN*PARM[[2]]*(1-exp(dat.t/PARM[[1]]))+SIGN,MIN=1,MAX=PARM[[2]]+1) # Inverse-Reverse Monomolecular
    TITLE <- "Inverse-Reverse Monomolecular"
  } else if(equation=="Inverse Monomolecular"){
    SIGN<- -1
    Trans.vec <- SCALE.vector(SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN,MIN=1,MAX=PARM[[2]]+1) # Inverse Monomolecular
    TITLE <- "Inverse Monomolecular"
  } else if(equation=="Monomolecular") {
    SIGN<- 1
    Trans.vec <- SCALE.vector(SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN,MIN=1,MAX=PARM[[2]]+1) # Monomolecular
    TITLE <- "Monomolecular"
  } else if(equation=="Reverse Monomolecular") {
    SIGN<- 1
    Trans.vec <- SCALE.vector(SIGN*PARM[[2]]*(1-exp(dat.t/PARM[[1]]))+SIGN,MIN=1,MAX=PARM[[2]]+1) # Reverse Monomolecular
    TITLE <- "Reverse Monomolecular"   
  } else if (equation=="Inverse Ricker") {
    SIGN <- -1
    Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN # Inverse Ricker
    Trans.vec <- SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec)))
    TITLE <- "Inverse Ricker" 
    
  }  else  {
    SIGN <- 1
    Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
    TITLE <- "Ricker"
  } 
  
  plot.data<-data.frame(dat.o,Trans.vec)
  
  p<- ggplot(plot.data,aes(x=dat.o,y=Trans.vec)) +
    ggtitle(paste(equation,"Transformation",sep="")) +
    theme_bw() +
    geom_line(size=1.5) +
    xlab(expression(bold("Original data values"))) +
    ylab(expression(bold("Tansformed data values"))) +
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
    scale_x_continuous(limits=c(min(dat.o),max(dat.o)*1.07)) +
    scale_y_continuous(limits=c(min(Trans.vec),max(Trans.vec)*1.07)) 
  
  if(file_test("-d",print.dir)){
    EXPORT.dir<-print.dir
    tiff(filename=paste0(print.dir,equation,"_Transformation_plot.tif"),width=160,height=150,units="mm",res=300,compression="lzw")
    print(p)
    dev.off()    
    return(p)
  }
  return(p)
  
}


############################################################################  
############ OPTIMIZATION FUNCTION USING GA STARTS: CONTINUOUS ############ 
############################################################################  
Resistance.Optimization_cont.nlm<-function(PARM,Resistance,equation, get.best,CS.inputs,Min.Max) {
  ID<-CS.inputs$ID
  ZZ<-CS.inputs$ZZ
  RESPONSE<-CS.inputs$RESPONSE
  #   CS.version<-CS.inputs$CS.version
  CS_Point.File<-CS.inputs$CS_Point.File
  CS.exe<-CS.inputs$CS.exe
  
  EXPORT.dir<-  if(get.best==FALSE){
    EXPORT.dir<-GA.inputs$Write.dir} else{
      EXPORT.dir<-GA.inputs$Results.dir
    }
  
  
  RESIST <- Resistance
  Name <- Resistance@data@names
  t1<-Sys.time()
  
  equation<-floor(equation)
  EQ <-get.EQ(equation)
  
  #   if(equation<7){
  SHAPE<-ifelse(exp(PARM[1])>10000,10000,exp(PARM[1])) # Upper boundaries on parameters
  Max.SCALE<-ifelse(exp(PARM[2])>10e6,10e6,exp(PARM[2]))
  cat("\n", Name, as.character(EQ),paste0("| Shape = ",SHAPE,"; "),paste0("Maximum scale = ",Max.SCALE,"\n"))  
  #   } else {
  #     Max.SCALE<-ifelse(exp(PARM[1])>10e6,log(10e6),exp(PARM[1]))
  #     OPT<-ifelse(exp(PARM[2])>25,log(25),exp(PARM[2])) # Upper boundaries on parameters
  #     SD<-ifelse(exp(PARM[3])>25,log(25),exp(PARM[3])) # Upper boundaries on parameters
  #     cat("\n", Name, as.character(EQ),paste0("| Max value = ",exp(PARM[1]),"; "),paste0("Gaus. peak = ",exp(PARM[2])),paste0("Gaus. sd = ",exp(PARM[3])),"\n")  
  #   }
  
  # Apply specified transformation
  if(equation==1){
    SIGN=-1
    RESIST <- SIGN*Max.SCALE*(1-exp(RESIST/SHAPE))+SIGN # Inverse-Reverse Monomolecular
    #     EQ <- "Inverse-Reverse Monomolecular"
    
  } else if(equation==2){
    SIGN=1
    R1 <- SIGN*Max.SCALE*(1-exp(RESIST/SHAPE))+SIGN # Reverse Monomolecular
    RESIST <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
    #     EQ <- "Reverse Monomolecular"  	    
    
  } else if(equation==3){
    SIGN=1
    RESIST <- SIGN*Max.SCALE*(1-exp(-1*RESIST/SHAPE))+SIGN # Monomolecular
    #     EQ <- "Monomolecular"
    
  } else if (equation==4) {
    SIGN=-1
    R1 <- SIGN*Max.SCALE*(1-exp(-1*RESIST/SHAPE))+SIGN # Inverse Monomolecular
    RESIST <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
    #     EQ <- "Inverse Monomolecular"  	    
    
  } else if (equation==5) {
    SIGN=-1
    R1 <- SIGN*(Max.SCALE*RESIST*exp(-1*RESIST/SHAPE))+SIGN # Inverse Ricker
    RESIST <- SCALE(R1,MIN=abs(cellStats(R1,stat='max')),MAX=abs(cellStats(R1,stat='min')))
    #     EQ <- "Inverse Ricker"  
    
  } else if (equation==6) {
    SIGN=1
    RESIST <- SIGN*(Max.SCALE*RESIST*exp(-1*RESIST/SHAPE))+SIGN #  Ricker
    #     EQ <- "Ricker"
    
  } else {
    RESIST <- (RESIST*0)+1 #  Distance
    #     EQ <- "Distance"    
  } # End if-else
  
  #   File.name <- paste0("exp",TRAN,"_MAX", MAX)
  File.name <- if (get.best==FALSE){
    File.name <- "optim_iter"
  } else {
    File.name <- Name
  }
  if(cellStats(RESIST,"max")>5e4)  RESIST<-SCALE(RESIST,1,5e4) # Rescale surface in case resistance are too high
  RESIST <- reclassify(RESIST, c(-Inf,1e-06, 1e-06,10e6,Inf,10e6))
  
  writeRaster(x=RESIST,filename=paste0(EXPORT.dir,File.name,".asc"), overwrite=TRUE)
  
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
  
  data<-cbind(ID,cs.matrix,RESPONSE)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=data$cs.matrix)
  
  # Fit model
  mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=data,REML=FALSE)
  mod$reTrms$Zt <- ZZ
  dfun <- do.call(mkLmerDevfun,mod)
  opt <- optimizeLmer(dfun)
  AIC.stat <- AIC(mkMerMod(environment(dfun), opt, mod$reTrms,fr = mod$fr))
  
  k<-length(PARM)+2
  AICc <- (AIC.stat)+(((2*k)*(k+1))/(nrow(CS.inputs$ID)-k-1))
  
  t2 <-Sys.time()
  cat(paste0("\t", "Iteration took ", round(t2-t1,digits=2), " seconds to complete"),"\n")
  cat(paste0("\t", "AIC = ",round(AIC.stat,3)),"\n")
  
  OPTIM.DIRECTION(Min.Max)*(AIC.stat) # Function to be minimized    
  
}
############################################

# Run Mixed effects models, recovery parameter estimates
#' Run maximum likelihood population effects mixed effects model (MLPE)
#' 
#' Runs MLPE as detailed by Clarke et al. (2002). This function is designed to generate a table of the fitted mixed effects model coefficients
#' 
#' @param resist.dir Directory containing .asc resistance files (optimized resistance surfaces)
#' @param genetic.dist Lower half of pairwise genetic distance matrix
#' @param out.dir If specified, a .csv table will printed to the specified directory (Default = NULL)
#' @return A table of MLPE fitted model coefficients

#' @export
#' @usage MLPE.lmm(resist.dir, genetic.dist, out.dir)
#' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

MLPE.lmm <- function(resist.dir, genetic.dist,out.dir=NULL){ 
  RESPONSE=genetic.dist
  resist.mat<-list.files(resist.dir,pattern="*_resistances.out",full.names=TRUE)
  resist.names<-gsub(pattern="_resistances.out","",x=list.files(resist.dir,pattern="*_resistances.out"))
  COEF.Table<-array()
  for(i in 1:length(resist.mat)){
    m<-length(read.table(resist.mat[i])[-1,-1])
    mm<-read.table(resist.mat[i])[-1,-1]
    ID<-To.From.ID(POPS=m)
    ZZ<-ZZ.mat(ID=ID)
    cs.matrix<-scale(mm[lower.tri(mm)],center=TRUE,scale=TRUE)
    
    dat<-cbind(ID,cs.matrix,RESPONSE)
    
    # Assign value to layer
    LAYER<-assign(resist.names[i],value=dat$cs.matrix)
    
    # Fit model
    mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=dat,REML=TRUE)
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
    (COEF.Table)
  }
}

##################################
#' Create diagnostic plots 
#' 
#' This function will generate mixed effect model diagnostic plots following optimization
#' 
#' @param cs.resistance.mat Path to CIRCUITSCAPE "_resistances.out" file
#' @param genetic.dist Vector of pairwise genetic distances (lower half of pairwise matrix). Can be input as CS.inputs$RESPONSE
#' @param XLAB Label for x-axis (Defaults to "Estimated resistance")
#' @param YLAB Label for y-axis (Defaults to "Genetic distance")
#' @param plot.dir Directory to output PDF of diagnostic plots
#' @return A four panel PDF including residual scatterplot, historgram of residuals, qqplot, and 

#' @export
#' @usage Diagnostic.Plots(resist.layer.path, genetic.dist, XLAB,YLAB, plot.dir)

Diagnostic.Plots<-function(cs.resistance.mat, genetic.dist, XLAB="Estimated resistance",YLAB ="Genetic distance",plot.dir){
  RESPONSE=genetic.dist
  NAME<-gsub(pattern="*_resistances.out","",x=(basename(cs.resistance.mat)))
  mm<-read.table(cs.resistance.mat)[-1,-1]
  m<-length(mm)
  ID<-To.From.ID(POPS=m)
  ZZ<-ZZ.mat(ID=ID)
  cs.matrix<-scale(lower(mm),center=TRUE,scale=TRUE)
  cs.unscale<-lower(mm)
  dat<-cbind(ID,cs.matrix,RESPONSE)
  
  # Assign value to layer
  LAYER<-assign("LAYER",value=dat$cs.matrix)
  
  # Fit model
  mod <- lFormula(RESPONSE ~ LAYER + (1|pop1), data=dat,REML=TRUE)
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
#' @param RESPONSE Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param CS_Point.File The path to the Circuitscape formatted point file. See Circuitscape documentation for help.
#' @param CS.exe The path to the CIRCUITSCAPE executable file (cs_run.exe). See details below. 
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
#' @return An R object that is a required input into optimization functions

#' @export
#' @usage CS.prep(n.POPS, RESPONSE, CS_Point.File, CS.exe, Neighbor.Connect)
#' @details \code{CS.exe} Example of path to CIRCUITSCAPE executible: 
#' 
#' '"C:/Program Files/Circuitscape/cs_run.exe"'
#'
#' ***NOTE: Double quotation used***
CS.prep <- function(n.POPS, RESPONSE=NULL,CS_Point.File,CS.exe,Neighbor.Connect=8){# Make to-from population list
  ID<-To.From.ID(n.POPS)
  ZZ<-ZZ.mat(ID)
  list(ID=ID,ZZ=ZZ,RESPONSE=RESPONSE,CS_Point.File=CS_Point.File,CS.exe=CS.exe,Neighbor.Connect=Neighbor.Connect,n.POPS=n.POPS)
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
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
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
#' @param run Number of consecutive generations without any improvement in the best fitness value before the GA is stopped (Default = 25)
#' @param keepBest A logical argument specifying if best solutions at each iteration should be saved (Default = TRUE)
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max' = Default). Optimization with \code{ga} maximizes the objective criteria
#' @return An R object that is a required input into optimization functions
#' 
#' @details Only files that you wish to optimize, either in isolation or simultaneously, should be included in the specified \code{ASCII.dir}. If you wish to optimize different combinations of surfaces, different directories contaiing these surfaces must be created.
#' 
#' \code{cont.shape} can take values of "Increase", "Decrease", or "Peaked". If you believe a resistance surface is related to your reposnse in a particular way, specifying this here may decrease the time to optimization. \code{cont.shape} is used to generate an initial set of parameter values to test during optimization. If specified, a greater proportion of the starting values will include your believed relatiosnship. If unspecified (the Default), a completely random set of starting values will be generated.
#' 
#' It is recommended to first run GA optimization with the default settings

#' @export
#' @autoImports
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
#' mutation = gaControl(type)$mutation)

GA.prep<-function(ASCII.dir,
                  Min.Max='max',
                  min.cat=0,
                  max.cat=2500, 
                  max.cont=2500,
                  cont.shape=NULL,
                  pop.mult = 15,
                  percent.elite = 0.05,
                  type= "real-valued",
                  pcrossover=0.85,
                  pmutation=0.1,
                  maxiter=1000,
                  run=25,
                  keepBest=TRUE,
                  population = gaControl(type)$population,
                  selection = gaControl(type)$selection,
                  crossover="gareal_blxCrossover",
                  mutation = gaControl(type)$mutation) {   
  
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
      max.list[[i]]<-c(7.99,15,max.cont)     
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
  
  list(parm.index=parm.index,ga.min=ga.min,ga.max=ga.max,surface.type=surface.type,parm.type=parm.type,Resistance.stack=r,n.layers=n.layers,layer.names=names,pop.size=pop.size, min.list=min.list,max.list=max.list, SUGGESTS=SUGGESTS,ASCII.dir=ASCII.dir, Results.dir=Results.dir, Write.dir=Write.dir,Plots.dir=Plots.dir,type= type, pcrossover=pcrossover, pmutation=pmutation, crossover=crossover, maxiter=maxiter, run=run, keepBest=keepBest, population=population,selection=selection,mutation=mutation,pop.mult = pop.mult, percent.elite = percent.elite,Min.Max=Min.Max)  
  
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
  equation<-floor(PARM[1])
  
  # Apply specified transformation
  if(equation==1){
    EQ <- "Inverse-Reverse Monomolecular"
    
  } else if(equation==2){
    EQ <- "Reverse Monomolecular"  	    
    
  } else if(equation==3){
    EQ <- "Monomolecular"
    
  } else if (equation==4) {
    EQ <- "Inverse Monomolecular"  	    
    
  } else if (equation==5) {
    EQ <- "Inverse Ricker"  
    
  } else if (equation==6) {
    EQ <- "Ricker"
    
  } else if (equation==7) {
    EQ <- "Inverse Gaussian"  	
    
  } else {
    EQ <- "Gaussian"  
  } # End if-else
  
  if(equation<7){
    df<-data.frame(EQ,PARM[2],PARM[3],NA)
    colnames(df)<-c("Equation","shape_opt","max","gaus.sd");row.names(df)<-NULL
    return(df)
  } else {
    df<-data.frame(EQ,PARM[2],PARM[3],PARM[4])
    colnames(df)<-c("Equation","shape_opt","max","gaus.sd");row.names(df)<-NULL
    return(df)
  }
}

Cont.Param_nG <- function(PARM){
  #   equation<-floor(PARM[1])
  
  #   # Apply specified transformation
  #   if(equation==1){
  #     EQ <- "Inverse-Reverse Monomolecular"
  #     
  #   } else if(equation==2){
  #     EQ <- "Reverse Monomolecular"        
  #     
  #   } else if(equation==3){
  #     EQ <- "Monomolecular"
  #     
  #   } else if (equation==4) {
  #     EQ <- "Inverse Monomolecular"  	    
  #     
  #   } else if (equation==5) {
  #     EQ <- "Inverse Ricker"  
  #     
  #   } else {
  #     EQ <- "Ricker"    
  #   
  #   } # End if-else
  
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
  Zl <- lapply(c("pop1","pop2"), function(nm) Matrix:::fac2sparse(ID[[nm]],"d", drop=FALSE))
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
  dec<-c(2,4)
  peak<-c(5,6)
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
      z<-c(runif(1,1,7.99),runif(1,.01,10),runif(1,1,100))
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
  equation=floor(equation)
  if(equation==1){
    EQ <- "Inverse-Reverse Monomolecular"
    
  } else if(equation==2){
    EQ <- "Reverse Monomolecular"        
    
  } else if(equation==3){
    EQ <- "Monomolecular"
    
  } else if (equation==4) {
    EQ <- "Inverse Monomolecular"  	    
    
  } else if (equation==5) {
    EQ <- "Inverse Ricker"  
    
  } else if (equation==6) {
    EQ <- "Ricker"
    
  } else {
    EQ <- "Distance"  	
  }
  #   } else {
  #     EQ <- "Gaussian"  
  #   } # End if-else
  (EQ)                                       
}

Result.txt <- function(GA.results, GA.inputs){
  SUMMARY<-summary(GA.results)
  
  sink(GA.inputs$Results.dir)
  cat(paste0("Summary from multisurface optimization run conducted on ",Sys.Date()),"\n")
  cat("\n")
  cat("\n",paste0("Surfaces included in optimization:"),"\n")
  cat(GA.inputs$parm.type$name,"\n")
  cat("Genetic Algorithm optimization settings:")
  cat("\n")
  cat(paste0("type: ", SUMMARY$type),"\n")
  cat(paste0("popSize: ", SUMMARY$popSize),"\n")
  cat(paste0("maxiter: ", SUMMARY$maxiter),"\n")
  cat(paste0("Retained elistism: ", SUMMARY$elistism),"\n")
  cat(paste0("pcrossover: ", SUMMARY$pcrossover),"\n")
  cat(paste0("pmutation: ", SUMMARY$pmutation),"\n")
  cat("\n")
  cat(paste0("The Genetic Algorithm completed after ",SUMMARY$iter," iterations"),"\n")
  cat("\n",paste0("Minimum AICc: -",SUMMARY$fitness),"\n")
  cat("\n",paste0("Optimized values for each surface:"),"\n")
  cat(SUMMARY$solution,"\n")
  sink()
}
# Optimiazation preparation
Optim.input<-function(Response,n.Pops,ASCII.dir,CS_Point.File,CS.exe,Neighbor.Connect=8,Constrained.Max=100,Initial.shape=c(seq(0.2,1,by=0.2),seq(1.25,10.75,by=0.75)),Bootstrap=FALSE,boot.iters=10000,Sample_Proportion=0.75){
  # Install necessary packages
  libs=c("raster", "lme4", "plyr")
  CheckInstallPackage(packages=libs)
  
  # Load libraries
  require(raster)
  require(lme4)
  require(plyr)
  require(GA)
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
