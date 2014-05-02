# Load libraries
require(raster)
require(lme4)
require(plyr)
require(GA)
require(devtools)
require(roxygen2)

install_github("roxygen3","hadley")

GA.inputs<-GA.prep(ASCII.dir="C:/Users/Bill/Dropbox/R_Functions/GA_optim/SimStudy/ASCII/Original/")

# Read in distance matrix
RESPONSE <- lower(matrix=read.csv("C:/Users/Bill/Dropbox/Research/Udzungwa/Analysis/rousset_scale_unique.csv",header=T))

cs.point <- "C:/Users/Bill/Dropbox/Research/Udzungwa/Analysis/FocalData/ASCII/sample_coords.txt"

cs <- paste('"C:/Program Files/Circuitscape/4.0/cs_run.exe"') # 4.0-Beta version


# Combine necessary function inputs
CS.inputs<-CS.prep(n.POPS=36,RESPONSE=RESPONSE,CS_Point.File=cs.point,CS.exe=cs)

#############################
# class       : RasterLayer 
# dimensions  : 40, 40, 1600  (nrow, ncol, ncell)
# resolution  : 0.025, 0.025  (x, y)
# extent      : 0, 1, 0, 1  (xmin, xmax, ymin, ymax)
# coord. ref. : NA 
# data source : C:\Users\Bill\Dropbox\R_Functions\GA_optim\SimStudy\ASCII\Original\cat_original.asc 
# names       : cat_original 

TEST<-GA.inputs$Resistance.stack[[1]]
require(spatstat)
r.spgrd<-as(TEST,"SpatialPointsDataFrame")
r.spgrd = r.spgrd[!is.na(r.spgrd[[1]]),]

selectedPoints = sample(1:length(r.spgrd[[1]]), 1000)
r.sampled = r.spgrd[selectedPoints,]

setwd("C:/Users/Bill/Dropbox/R_Functions/Git/Packages/ResistanceGA")
