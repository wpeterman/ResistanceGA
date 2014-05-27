require(RandomFields)
require(ResistanceGA)
require(ggplot2)


rm(list = ls())
set.seed(321123)
###################################
# Where are ASCII files and batch files going to be written to?
if("ResistanceGA_Examples"%in%dir("C:/")==FALSE) 
  dir.create(file.path("C:/", "ResistanceGA_Examples")) 

# Create a subdirectory for the first example
dir.create(file.path("C:/ResistanceGA_Examples/","SingleSurface")) 

write.dir <- "C:/ResistanceGA_Examples/SingleSurface/"      # Directory to write .asc files and results

# Random point and raster parameters
r.dim <- 50       # number of cells on a side
cell.size <- 0.025        # raster cell dimension     
min.point <- 0.25*(r.dim*cell.size)       # minimum coordinate for generating random points (multiplied by 0.25 to prevent edge effects)
max.point <- (r.dim*cell.size)-min.point        # maximum coordinate for generating random points

# Number of "Sample locations" to generate. This example will generate points on a square grid, so choose a number that has an even square root
n <- 25 

x <- seq(from=min.point,max.point,by=cell.size)       # set x & y range to draw random samples from
y <- seq(from=min.point,max.point,by=cell.size)
xy <- cbind(x,y)
xy.coords <- SpatialPoints(xy)
Sample.coord <- spsample(Spatial(bbox=bbox(SpatialPoints(xy))), n, type="regular") # Generate regularly spaced points on a grid
coord.id <- cbind((1:n),Sample.coord@coords)       # Combine location ID with coordinates

# Write the table to a file. This is formatted for input into CIRCUITSCAPE
write.table(coord.id,file=paste0(write.dir,"samples.txt"),sep="\t",col.names=F,row.names=F)

#####################################################################
#####################################################################
# Using random fields create two resistance surfaces
model<-RMexp() +
  RMtrend(mean=10)

grid.vars <- GridTopology(cellcentre.offset=c(cell.size/2, cell.size/2),
                          cellsize=c(cell.size, cell.size),
                          cells.dim=rep(r.dim,2))

rf.sim <- RFsimulate(model, x=grid.vars) # Simulate 2 surfaces

cont.rf <- raster(rf.sim) # Define the first as a continuous surface
names(cont.rf)<-"cont"
plot(cont.rf)
plot(Sample.coord, pch=16, col="blue", add=TRUE) # Add randomly generated points


# Write original rasters to file for use with CIRCUITSCAPE
writeRaster(cont.rf,filename=paste0(write.dir,"cont.asc"),overwrite=TRUE)

#########################################
# Run prep functions
# GA.inputs<-GA.prep(ASCII.dir=write.dir,
#                    pop.mult=10,
#                    min.cat=0,
#                    max.cat=500,
#                    max.cont=500,
#                    run=1) # Only two runs selected...THIS WILL NOT OPTIMIZE, done for demostration only

GA.inputs<-GA.prep(ASCII.dir=write.dir,
                   min.cat=0,
                   max.cat=500,
                   max.cont=500,
                   seed = 5)

CS.inputs<-CS.prep(n.POPS=n,
                   CS_Point.File=paste0(write.dir,"samples.txt"),
                   CS.exe=paste('"C:/Program Files/Circuitscape/4.0/cs_run.exe"')) # Note that RESPONSE is omittted because it has not been made yet

# Monomolecular = equation # 3
PARM=c(3,2,275)
Resist<-Resistance.tran(transformation="Monomolecular",shape=2,max=275,r=cont.rf) # Make Combine_Surfaces so that it can take both an R raster object or read a .asc file

# Resist.false<-Resistance.tran(transformation="Monomolecular",shape=2,max=100,r=cont.rf) # Make Combine_Surfaces so that it can take both an R raster object or read a .asc file


plot.t<-PLOT.trans(PARM=c(2,275),Resistance="C:/ResistanceGA_Examples/SingleSurface/cont.asc",transformation="Monomolecular") #print.dir="C:/ResistanceGA_Example/Results/Plots/"

# Run CIRCUITSCAPE to generate pairwise matrix of effective resistance distance
# Only continuous the surface will affect resistance in the first example
CS.Resist<- Run_CS(CS.inputs=CS.inputs,GA.inputs=GA.inputs,r=Resist)

# We add some random noise to the response (i.e. CS resistance output)
NOISE <- rnorm(n=length(CS.Resist),mean=0,(0.005*max(CS.Resist))) # Generate random noise matrix
CS.response<-round((CS.Resist+NOISE),digits=4) # Add random noise to resistance output. Use this as the RESPONSE for single surface testing
plot(CS.response~CS.Resist)
write.table(CS.response,paste0(write.dir,"CS.response.txt"),col.names=F,row.names=F)

# Remake the CS.inputs file, including the newly created CS.Response
CS.inputs<-CS.prep(n.POPS=n,
                   RESPONSE=CS.Resist,
                   CS_Point.File=paste0(write.dir,"samples.txt"),
                   CS.exe=paste('"C:/Program Files/Circuitscape/4.0/cs_run.exe"'))

# Grid search
Grid.Results <- Grid.Search(shape=seq(1,4,by=0.1),max=seq(50,500,by=75),transformation="Monomolecular",Resistance=cont.rf, CS.inputs)

filled.contour(Grid.Results$Plot.data,col=topo.colors(20),xlab="Shape parameter",ylab="Maximum value parameter")

# Best from Grid.Search
Grid.Results$AICc[match(min(Grid.Results$AICc$AICc),Grid.Results$AICc$AICc),]

# Data generating values
Grid.Results$AICc[match(interaction(2,275),interaction(Grid.Results$AICc[,c(1,2)])),]

#############################
# Alternatively
svg("C:/Users/Bill/Dropbox/R_Functions/Git/Packages/ResistanceGA/figure/Grid.Surface.svg",width=6,height=6)
filled.contour(Grid.Results$Plot.data,col=topo.colors(20),xlab="Shape parameter",ylab="Maximum value parameter")
dev.off()

png("C:/Users/Bill/Dropbox/R_Functions/Git/Packages/ResistanceGA/figure/Grid.Surface.png",units="in", width=4,height=4,res=150)
filled.contour(Grid.Results$Plot.data,col=topo.colors(20),xlab="Shape parameter",ylab="Maximum value parameter")
dev.off()
#############################


# Single surface optimization
system.time(SS_RESULTS<-SS_optim(CS.inputs=CS.inputs,
         GA.inputs=GA.inputs))



#############################
set.seed(5)
if("ResistanceGA_Examples"%in%dir("C:/")==FALSE) 
  dir.create(file.path("C:/", "ResistanceGA_Examples")) 

# Create a subdirectory for the second example
dir.create(file.path("C:/ResistanceGA_Examples/","MultipleSurfaces")) 

write.dir <- "C:/ResistanceGA_Examples/MultipleSurfaces/"      # Directory to write .asc files and results


# Simulate another continuous surface. This will be converted into a 3-class categorical surface
rf.sim <- RFsimulate(model, x=grid.vars,n=2) # Create two surfaces

cont.rf <- raster(rf.sim[1]) # Define first as a continuous surface
names(cont.rf)<-"cont"

plot(cont.rf)
plot(Sample.coord, pch=16, col="blue", add=TRUE)

cat.rf <- raster(rf.sim[2]) 
names(cat.rf) <- "cat"
cat.cut <- summary(cat.rf) # Define quartiles, use these to define categories

cat.rf[cat.rf<=cat.cut[2]] <- 0
cat.rf[cat.rf>0 & cat.rf<=cat.cut[4]] <- 1
cat.rf[!cat.rf%in%c(0,1)] <- 2
plot(cat.rf)
plot(Sample.coord, pch=16, col="blue", add=TRUE)

##############################
# Make categorical feature class (like a road)
feature <- matrix(0,r.dim,r.dim)
feature[25,] <- 1
feature[,28:29] <- 1
feature <- raster(feature)
extent(feature)<-extent(cat.rf)
plot(feature)
names(feature)<-"feature"
plot(feature)
plot(Sample.coord, pch=16, col="blue", add=TRUE)

# Write original rasters to file for use with CIRCUITSCAPE
writeRaster(cat.rf,filename=paste0(write.dir,"cat.asc"),overwrite=TRUE)
writeRaster(cont.rf,filename=paste0(write.dir,"cont.asc"),overwrite=TRUE)
writeRaster(feature,filename=paste0(write.dir,"feature.asc"),overwrite=TRUE)

# Write the table to a file. This is formatted for input into CIRCUITSCAPE
write.table(coord.id,file=paste0(write.dir,"samples.txt"),sep="\t",col.names=F,row.names=F)

#############################
# Visualize transformation of continuous surface. The first term of the PARM function refers to the shape parameter, and the second term refers to the maximum value parameter.

plot.t<-PLOT.trans(PARM=c(2,250),Resistance="C:/ResistanceGA_Examples/MultipleSurfaces/cont.asc",transformation="Inverse-Reverse Monomolecular") #print.dir="C:/ResistanceGA_Example/Results/Plots/"

# Combine raster surfaces, apply transformation to continuous surface and change values of categorical and feature surfaces
PARM=c(0,150,50,1,2,250,0,400)


# GA.inputs<-GA.prep(ASCII.dir=write.dir,
#                    pop.mult=5,
#                    min.cat=0,
#                    max.cat=500,
#                    max.cont=500,
#                    run=1) # Only single run selected...THIS WILL NOT OPTIMIZE, done for demostration only

GA.inputs<-GA.prep(ASCII.dir=write.dir,
                   min.cat=0,
                   max.cat=500,
                   max.cont=500,
                   seed = 5)

# Combine resistance surfaces
Resist<-Combine_Surfaces(PARM=PARM,CS.inputs=CS.inputs,GA.inputs=GA.inputs) # Make Combine_Surfaces so that it can take both an R raster object or read a .asc file

# Generate new CS response surface by using Run_CS
CS.response<- Run_CS(CS.inputs=CS.inputs,GA.inputs=GA.inputs,r=Resist)

# # Generate some random noise and add it to the resistance surface
# NOISE <- rnorm(n=length(CS.Resist),mean=0,(0.05*max(CS.Resist)))
# plot(Resist)
# 
# CS.response<-round((CS.Resist+NOISE),digits=4)
# plot(CS.response~CS.Resist)
write.table(CS.response,file=paste0(write.dir,"Combined_response.csv"),sep=",",row.names=F,col.names=F)

# Run prep functions
CS.inputs<-CS.prep(n.POPS=n,
                      RESPONSE=CS.response,
                      CS_Point.File=paste0(write.dir,"samples.txt"),
                      CS.exe=paste('"C:/Program Files/Circuitscape/4.0/cs_run.exe"'))

system.time(Multi.Surface_optim <-MS_optim(CS.inputs=CS.inputs,GA.inputs=GA.inputs))

mantel(CS.Resist~CS.TRUTH)
plot(CS.Resist~CS.TRUTH)


rm(list=".Random.seed", envir=globalenv()) 
