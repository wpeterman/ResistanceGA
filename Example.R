#'@examples
#'\dontrun{
#' Code to conduct analysis of GA to optimize multiple surface simultaneously
require(RandomFields)
require(ResistanceGA)

rm(list = ls())
set.seed(123321)
###################################
# Where are ASCII files and batch files going to be written to?
if("ResistanceGA_Example"%in%dir("C:/")==FALSE) dir.create(file.path("C:/", "ResistanceGA_Example")) 

write.dir <- "C:/ResistanceGA_Example/" # Create a directory to write .asc files and results

# Random point and raster parameters
r.dim <- 60 # number of cells on a side
cell.size <- 0.025
min.point <- 0.25*(r.dim*cell.size) # minimum coordinate for generating random points (multiplied by 0.25 to prevent edge effects)
max.point <- (r.dim*cell.size)-min.point # maximum coordinate for generating random points
n <- 50 # Number of random "Sample locations" to generate
x <- seq(from=min.point,max.point,by=cell.size) # set x & y range to draw random samples from
y <- seq(from=min.point,max.point,by=cell.size)
xy <-cbind(x,y)
SAMPLE<-SpatialPoints(xy)
COORD <- spsample(Spatial(bbox=bbox(SAMPLE)), n, type="random")
coord.id <-cbind((1:n),COORD@coords)

write.table(coord.id,file=paste0(write.dir,"samples.txt"),sep="\t",col.names=F,row.names=F)

#####################################################################
#####################################################################
# Create random field to represent continuous and categorical landscape surface
model<-RMexp() +
  RMtrend(mean=5)

grid.vars <- GridTopology(cellcentre.offset=c(cell.size/2, cell.size/2),
                          cellsize=c(cell.size, cell.size),
                          cells.dim=rep(r.dim,2))

rf.sim <- RFsimulate(model, x=grid.vars, n=3) # Simulate 3 surfaces

cont.rf <- raster(rf.sim[1]) # Define the first as a continuous surface
names(cont.rf)<-"cont"
plot(cont.rf)
plot(COORD, pch=16, col="blue", add=TRUE)

cat.rf <- raster(rf.sim[2]) # Make the second surface a 3-category surface
names(cat.rf)<-"cat"
cat.cut <- summary(cat.rf) # Define quartiles

cat.rf[cat.rf<=cat.cut[2]] <- 0
cat.rf[cat.rf>0 & cat.rf<=cat.cut[4]] <- 1
cat.rf[!cat.rf%in%c(0,1)] <- 2
plot(cat.rf)
plot(COORD, pch=16, col="blue", add=TRUE)

# rand.rf <- raster(rf.sim[3]) # Define the third as a continuous surface
# names(rand.rf)<-"rand"
# plot(rand.rf)
# plot(COORD, pch=16, col="blue", add=TRUE)

##############################
# Make landscape feature
feature <- matrix(0,r.dim,r.dim)
feature[25,] <- 1
feature[,30:32] <- 1
feature <- raster(feature)
extent(feature)<-extent(cat.rf)
plot(feature)
names(feature)<-"feature"
plot(feature)
plot(COORD, pch=16, col="blue", add=TRUE)

# Make random noise surface
# model<-RMexp(scale=0.01) + 
#   RMtrend(mean=0)
# 
# grid.vars <- GridTopology(cellcentre.offset=c(cell.size/2, cell.size/2),
#                           cellsize=c(cell.size, cell.size),
#                           cells.dim=rep(r.dim,2))
# 
# NOISE <- RFsimulate(model, x=grid.vars, n=1) # Simulate 1 random noise surface
# NOISE<-raster(NOISE)
# plot(NOISE)

# Write original rasters to file for use with CIRCUITSCAPE
writeRaster(cat.rf,filename=paste0(write.dir,"cat.asc"),overwrite=TRUE)
writeRaster(cont.rf,filename=paste0(write.dir,"cont.asc"),overwrite=TRUE)
# writeRaster(rand.rf,filename=paste0(write.dir,"rand.asc"),overwrite=TRUE)
writeRaster(feature,filename=paste0(write.dir,"feature.asc"),overwrite=TRUE)

#########################################
# Run prep functions
GA.inputs<-GA.prep(ASCII.dir=write.dir,
                   pop.mult=5,
                   min.cat=0,
                   max.cat=500,
                   max.cont=500,
                   run=1) # Only single run selected...THIS WILL NOT OPTIMIZE, done for demostration only

CS.inputs<-CS.prep(n.POPS=n,
                   CS_Point.File=paste0(write.dir,"samples.txt"),
                   CS.exe=paste('"C:/Program Files/Circuitscape/4.0/cs_run.exe"')) # Note that RESPONSE is omittted because it has not been made yet

# Run CIRCUITSCAPE to generate pairwise matrix of effective resistance distance
# Only continuous surface affects resistance in the first example
CS.Resist<- Run_CS(CS.inputs=CS.inputs,GA.inputs=GA.inputs,r=cont.rf)

NOISE <- rnorm(n=length(CS.Resist),mean=0,(0.05*max(CS.Resist))) # Generate random noise matrix
CS.response<-round((CS.Resist+NOISE),digits=4) # Add random noise to resistance output. Use this as the RESPONSE for single surface testing
write.table(CS.response,paste0(write.dir,"CS.response.txt"),col.names=F,row.names=F)


CS.inputs<-CS.prep(n.POPS=n,
                   RESPONSE=CS.response,
                   CS_Point.File=paste0(write.dir,"samples.txt"),
                   CS.exe=paste('"C:/Program Files/Circuitscape/4.0/cs_run.exe"'))

# Single surface optimization
SS_optim(CS.inputs=CS.inputs,
         GA.inputs=GA.inputs)





#############################
# Test and visualize transformation of continuous surface
dat<-seq(0,10,length.out=100) # data ranging from 0-10
SIGN=1
Max.SCALE=5
SHAPE=2
Tran<-SIGN*Max.SCALE*(1-exp(-1*dat/SHAPE)) # Monomolecular transformation
plot(Tran~dat)

# Combine raster surfaces, apply transformation to continuous surface and change values of categorical and feature surfaces
PARM=c(0,50,100,3,2,100,0,100)
# PARM=c(0, # First feature of categorical
#        50, # Second feature of categorical
#        100, # Third feature of categorical
#        3, # Transformation equation for continuous surface
#        2, #  Shape parameter
#        100, # Scale parameter
#        0, # First feature of feature surface
#        100) # Second feature of feature surface

PARM=c(0,50,100,3,2,100,0,100)

Resist<-Combine_Surfaces(PARM=PARM,CS.inputs=CS.inputs,GA.inputs=Sim.GA)
# Resist.Noise<-Resist+NOISE
CS.Resist<- Run_CS(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,r=Resist)
CS.TRUTH<- Run_CS(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,r=Resist)

NOISE <- rnorm(n=length(CS.Resist),mean=0,(0.15*max(CS.Resist)))
plot(Resist)

CS.Resist<-round((CS.Resist+NOISE),digits=4)

write.table(CS.Resist,file=paste0(write.dir,"RESPONSE_Med.csv"),sep=",",row.names=F,col.names=F)

Model.Params<-Function.Params(n.POPS=nrow(coord.id),RESPONSE=CS.Resist,CS.version=cs.vers,CS_Point.File=cs.point,CS.exe=cs,EXPORT.dir=exp.dir)

AIC.Med <-Resistance_multi.GA_noGaus(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,Min.Max="min")
mantel(CS.Resist~CS.TRUTH)
plot(CS.Resist~CS.TRUTH)

# Medium High
# CONTRAST=150
# PARM=c(0,3,CONTRAST,3,2,CONTRAST,0,CONTRAST)
# Resist<-Combine_Surfaces(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA)
# 
# CS.Resist<- Run_CS(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA)
# 
# CS.Resist<-CS.Resist+NOISE
# 
# write.table(CS.Resist,file=paste0(write.dir,"RESPONSE_MedHigh.csv"),sep=",",row.names=F,col.names=F)
# 
# Model.Params<-Function.Params(n.POPS=nrow(coord.id),RESPONSE=CS.Resist,CS.version=cs.vers,CS_Point.File=cs.point,CS.exe=cs,EXPORT.dir=exp.dir)
# 
# AIC.MedHigh <-Resistance_multi.GA_noGaus(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,Min.Max="min")

#  High Noise
CONTRAST=100
PARM=c(0,50,CONTRAST,3,2,CONTRAST,0,CONTRAST,7,1,1)
Resist<-Combine_Surfaces.nG(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA)
# Resist.Noise<-Resist+NOISE
CS.Resist<- Run_CS(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,r=Resist)
CS.TRUTH<- Run_CS(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,r=Resist)

NOISE <- rnorm(n=length(CS.Resist),mean=0,(0.25*max(CS.Resist)))
plot(Resist)

CS.Resist<-round((CS.Resist+NOISE),digits=4)

write.table(CS.Resist,file=paste0(write.dir,"RESPONSE_High.csv"),sep=",",row.names=F,col.names=F)

Model.Params<-Function.Params(n.POPS=nrow(coord.id),RESPONSE=CS.Resist,CS.version=cs.vers,CS_Point.File=cs.point,CS.exe=cs,EXPORT.dir=exp.dir)

AIC.High <-Resistance_multi.GA_noGaus(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,Min.Max="min")
mantel(CS.Resist~CS.TRUTH)
plot(CS.Resist~CS.TRUTH)
#####################################
#  TRUTH
# CONTRAST=100
# PARM=c(0,50,CONTRAST,3,2,CONTRAST,0,CONTRAST,7,1,1)
# Resist<-Combine_Surfaces.nG(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA)
# Resist.Noise<-Resist+NOISE
# CS.Resist<- Run_CS(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,r=Resist)
# CS.Resist<-round(CS.Resist,digits=4)
# 
# write.table(CS.Resist,file=paste0(write.dir,"RESPONSE_TRUTH.csv"),sep=",",row.names=F,col.names=F)
# 
# Model.Params<-Function.Params(n.POPS=nrow(coord.id),RESPONSE=CS.Resist,CS.version=cs.vers,CS_Point.File=cs.point,CS.exe=cs,EXPORT.dir=exp.dir)
# 
# AIC.TRUTH <-Resistance_multi.GA_noGaus(PARM=PARM,Model.Params=Model.Params,GA.params=Sim.GA,Min.Max="min")
##################################
##################################
# LEVEL <- c("Low", "Med.Low", "Med", "Med.High","High")
# AIC.val <- c(AIC.Low,AIC.MedLow,AIC.Med,AIC.MedHigh,AIC.High)
LEVEL <- c("Low", "Med", "High")
AIC.val <- c(AIC.Low,AIC.Med,AIC.High)
(Sim.AIC.Table <- data.frame(LEVEL,AIC.val))
write.table(Sim.AIC.Table,file=paste0(write.dir,"Sim_AIC_Values.txt"),row.names=F,col.names=T,sep="\t")

rm(list=".Random.seed", envir=globalenv()) 
