r1<-Resistance.tran(transformation="Ricker",shape=2,max=100,r=Udz.rast)

r2<-Resistance.tran(transformation="Reverse Ricker",shape=2,max=100,r=Udz.rast)

(TEST<-matrix(seq(from=0,to=10,length.out=100),nrow=10))
# TEST[c(1:6),c(1,2)]<-NA
# TEST[c(6:8),c(7,8)]<-NA
TEST
(r.TEST<-rev(TEST))
(TEST.r<-matrix(r.TEST,8))

########################
plot(mono<-Resistance.tran(transformation="Monomolecular",shape=2,max=50,r=Udz.rast))

plot(r.mono<-Resistance.tran(transformation="Reverse Monomolecular",shape=2,max=50,r=Udz.rast))

plot(rick<-Resistance.tran(transformation="Ricker",shape=2,max=50,r=Udz.rast))

plot(r.rick<-Resistance.tran(transformation="Reverse Ricker",shape=2,max=50,r=Udz.rast))

plot(r.rick<-Resistance.tran(transformation="Inverse Ricker",shape=2,max=50,r=Udz.rast))

plot(r.rick<-Resistance.tran(transformation="Inverse-Reverse Ricker",shape=2,max=50,r=Udz.rast))

#######################

par(mfrow=c(1,1))
PLOT.trans(PARM=c(2,100),Resistance=Udz.rast,transformation="Monomolecular")
PLOT.trans(PARM=c(2,100),Resistance=Udz.rast,transformation="Reverse Monomolecular")
PLOT.trans(PARM=c(2,100),Resistance=Udz.rast,transformation="Inverse Monomolecular")
PLOT.trans(PARM=c(2,100),Resistance=Udz.rast,transformation="Inverse-Reverse Monomolecular")


plot(test.rast<-raster(TEST))
plot(mono.test<-Resistance.tran(transformation="Monomolecular",shape=2,max=50,r=test.rast))
plot(Resistance.tran(transformation="Reverse Monomolecular",shape=2,max=50,r=test.rast))
plot(Resistance.tran(transformation="Inverse Monomolecular",shape=2,max=50,r=test.rast))
plot(Resistance.tran(transformation="Inverse-Reverse Monomolecular",shape=2,max=50,r=test.rast))
plot(mono.test<-Resistance.tran(transformation="Ricker",shape=2,max=50,r=test.rast))
plot(Resistance.tran(transformation="Reverse Ricker",shape=2,max=50,r=test.rast))
plot(Resistance.tran(transformation="Inverse Ricker",shape=2,max=50,r=test.rast))
plot(Resistance.tran(transformation="Inverse-Reverse Ricker",shape=2,max=50,r=test.rast))

# Manual reversal
test.mat<-rev(mono.test)
plot(setValues(mono.test,values=test.mat))

# Check inverse function
plot(test.rast<-raster(TEST))
plot(mono.test<-Resistance.tran(transformation="Monomolecular",shape=2,max=50,r=test.rast))
plot(inv.mono<-Resistance.tran(transformation="Inverse Monomolecular",shape=2,max=50,r=test.rast))
plot(inv_rev.mono<-Resistance.tran(transformation="Inverse-Reverse Monomolecular",shape=2,max=50,r=test.rast))
####################################
library(ggplot2)
library(gtable)
library(gridExtra)

# Get the widths
scale_x_continuous(limits=c(min(original),max(original)),breaks=x.break) +
  scale_y_continuous(limits=c(min(transformed),max(transformed)),breaks=y.break) 

TEST<-matrix(seq(from=0,to=10,length.out=100),nrow=10)
test.rast<-raster(TEST)

(g1 <- PLOT.trans(PARM=c(2,100),Resistance=test.rast,transformation="Monomolecular"))
(g2 <- PLOT.trans(PARM=c(2,100),Resistance=test.rast,transformation="Reverse Monomolecular"))
(g3 <- PLOT.trans(PARM=c(2,100),Resistance=test.rast,transformation="Inverse Monomolecular"))
(g4 <- PLOT.trans(PARM=c(2,100),Resistance=test.rast,transformation="Inverse-Reverse Monomolecular"))
(g5 <- PLOT.trans(PARM=c(2,100),Resistance=test.rast,transformation="Ricker"))
(g6 <- PLOT.trans(PARM=c(2,100),Resistance=test.rast,transformation="Reverse Ricker"))
(g7 <- PLOT.trans(PARM=c(2,100),Resistance=test.rast,transformation="Inverse Ricker"))
(g8 <- PLOT.trans(PARM=c(2,100),Resistance=test.rast,transformation="Inverse-Reverse Ricker"))


# Arrange the 8 charts
svg("C:/Users/Bill/Dropbox/R_Functions/Git/Packages/ResistanceGA/figure/Transformations.svg",width=10,height=13.33)
grid.arrange(g1, g2, g3, g4, g5,g6,g7,g8, nrow=4)
dev.off()



Resistance.Opt_single(PARM=c(3,2,100),
                      Resistance=CS.Resist,
                      CS.inputs,GA.inputs,Min.Max='max',iter=1)

##################################
Grid.Results<-Grid.Search(shape=seq(1,5,by=0.1),max=seq(50,750,by=50),transformation="Monomolecular",Resistance=cont.rf,CS.inputs)


Grid.Results<-Grid.Search(shape=seq(1,10,by=.25),max=seq(25,500,by=25),transformation="Monomolecular",Resistance=cont.rf,CS.inputs)

# Alternatively
svg("C:/Users/Bill/Dropbox/R_Functions/Git/Packages/ResistanceGA/figure/Grid.Surface.svg",width=4,height=4)
filled.contour(Grid.Results$Plot.data,col=topo.colors(30),xlab="Shape parameter",ylab="Maximum value")
dev.off()


# filled.contour(Grid.Results,col=terrain.colors(30),xlab="Shape parameter",ylab="Maximum value")

Grid.Results$AICc[match(min(Grid.Results$AICc$AICc),Grid.Results$AICc$AICc),]

Grid.Results$AICc[match(interaction(2,250),interaction(Grid.Results$AICc[,c(1,2)])),]


# Good
wireframe(abs(z) ~ y*x, data = Grid.Results,shade=FALSE, aspect=c(1,1),
          xlab = "Max value", ylab = "Shape parameter",zlab=list("abs(AICc)",rot=90),
          main = "Absolute AICc across Max*Shape parameter space",
          drape = TRUE,
          scales=list(arrows=F),
          screen = list(z=-20, x=-85, y=-10)
)
