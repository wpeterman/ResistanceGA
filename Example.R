require(RandomFields)
require(ResistanceGA)
require(ggplot2)


rm(list = ls())
set.seed(12345)
###################################
# Where are ASCII files and batch files going to be written to?
if("ResistanceGA_Example"%in%dir("C:/")==FALSE) dir.create(file.path("C:/", "ResistanceGA_Example")) 

write.dir <- "C:/ResistanceGA_Example/" # Create a directory to write .asc files and results

# Random point and raster parameters
r.dim <- 60 # number of cells on a side
cell.size <- 0.025
min.point <- 0.25*(r.dim*cell.size) # minimum coordinate for generating random points (multiplied by 0.25 to prevent edge effects)
max.point <- (r.dim*cell.size)-min.point # maximum coordinate for generating random points
n <- 49 # Number of "Sample locations" to generate. This example will generate points on a square grid, so choose a number that has an even square root
x <- seq(from=min.point,max.point,by=cell.size) # set x & y range to draw random samples from
y <- seq(from=min.point,max.point,by=cell.size)
# xy <-c(min.point,max.point,min.point,max.point)
# SAMPLE<-rSSI(.05,n=n,win=xy,giveup=1000) 
xy <-cbind(x,y)
SAMPLE<-SpatialPoints(xy)
COORD <- spsample(Spatial(bbox=bbox(SAMPLE)), n, type="regular") # Generate regularly spaced points on a grid
min(lower(pointDistance(COORD,lonlat=F)))
coord.id <-cbind((1:n),COORD@coords)

write.table(coord.id,file=paste0(write.dir,"samples.txt"),sep="\t",col.names=F,row.names=F)

#####################################################################
#####################################################################
# Using random fields create two resistance surfaces
model<-RMexp() +
  RMtrend(mean=10)

grid.vars <- GridTopology(cellcentre.offset=c(cell.size/2, cell.size/2),
                          cellsize=c(cell.size, cell.size),
                          cells.dim=rep(r.dim,2))

rf.sim <- RFsimulate(model, x=grid.vars, n=2) # Simulate 2 surfaces

cont.rf <- raster(rf.sim[1]) # Define the first as a continuous surface
names(cont.rf)<-"cont"
plot(cont.rf)
plot(COORD, pch=16, col="blue", add=TRUE) # Add randomly generated points

cat.rf <- raster(rf.sim[2]) # Make the second surface a 3-class categorical surface
names(cat.rf)<-"cat"
cat.cut <- summary(cat.rf) # Define quartiles, use these to define categories

cat.rf[cat.rf<=cat.cut[2]] <- 0
cat.rf[cat.rf>0 & cat.rf<=cat.cut[4]] <- 1
cat.rf[!cat.rf%in%c(0,1)] <- 2
plot(cat.rf)
plot(COORD, pch=16, col="blue", add=TRUE)

##############################
# Make categorical feature class (like a road)
feature <- matrix(0,r.dim,r.dim)
feature[25,] <- 1
feature[,30:31] <- 1
feature <- raster(feature)
extent(feature)<-extent(cat.rf)
plot(feature)
names(feature)<-"feature"
plot(feature)
plot(COORD, pch=16, col="blue", add=TRUE)

# Write original rasters to file for use with CIRCUITSCAPE
writeRaster(cat.rf,filename=paste0(write.dir,"cat.asc"),overwrite=TRUE)
writeRaster(cont.rf,filename=paste0(write.dir,"cont.asc"),overwrite=TRUE)
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
# Only continuous the surface will affect resistance in the first example
CS.Resist<- Run_CS(CS.inputs=CS.inputs,GA.inputs=GA.inputs,r=cont.rf)

# We add some random noise to the response (i.e. CS resistance output)
NOISE <- rnorm(n=length(CS.Resist),mean=0,(0.05*max(CS.Resist))) # Generate random noise matrix
CS.response<-round((CS.Resist+NOISE),digits=4) # Add random noise to resistance output. Use this as the RESPONSE for single surface testing
plot(CS.response~CS.Resist)
write.table(CS.response,paste0(write.dir,"CS.response.txt"),col.names=F,row.names=F)

# Remake the CS.inputs file, including the newly created CS.Response
CS.inputs<-CS.prep(n.POPS=n,
                   RESPONSE=CS.response,
                   CS_Point.File=paste0(write.dir,"samples.txt"),
                   CS.exe=paste('"C:/Program Files/Circuitscape/4.0/cs_run.exe"'))

# Single surface optimization
SS_optim(CS.inputs=CS.inputs,
         GA.inputs=GA.inputs)





#############################
# Visualize transformation of continuous surface. The first term of the PARM function refers to the shape parameter, and the second term refers to the maximum value parameter.

PLOT.trans(PARM=c(2,10),Resistance="C:/ResistanceGA_Example/cont.asc",equation="Reverse Monomolecular",print.dir="C:/ResistanceGA_Example/Results/Plots/")

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

PARM=c(0,50,150,3,2,300,0,500)

# Combine resistance surfaces
Resist<-Combine_Surfaces(PARM=PARM,CS.inputs=CS.inputs,GA.inputs=GA.inputs) # Make Combine_Surfaces so that it can take both an R raster object or read a .asc file

# Generate new CS response surface by using Run_CS
CS.Resist<- Run_CS(CS.inputs=CS.inputs,GA.inputs=GA.inputs,r=Resist)

# Generate some random noise nad add it to the resistance surface
NOISE <- rnorm(n=length(CS.Resist),mean=0,(0.1*max(CS.Resist)))
plot(Resist)

CS.response<-round((CS.Resist+NOISE),digits=4)
plot(CS.response~CS.Resist)
write.table(CS.response,file=paste0(write.dir,"ombined_response.csv"),sep=",",row.names=F,col.names=F)

# Run prep functions
GA.inputs<-GA.prep(ASCII.dir=write.dir,
                   pop.mult=5,
                   min.cat=0,
                   max.cat=500,
                   max.cont=750,
                   run=1) # Only single run selected...THIS WILL NOT OPTIMIZE, done for demostration only

CS.inputs<-CS.prep(n.POPS=n,
                      RESPONSE=CS.response,
                      CS_Point.File=paste0(write.dir,"samples.txt"),
                      CS.exe=paste('"C:/Program Files/Circuitscape/4.0/cs_run.exe"'))

Multi.Surface_optim <-MS_optim(CS.inputs=CS.inputs,GA.inputs=GA.inputs)

mantel(CS.Resist~CS.TRUTH)
plot(CS.Resist~CS.TRUTH)


rm(list=".Random.seed", envir=globalenv()) 
