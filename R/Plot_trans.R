#' Plot continuous surface transformation
#'
#' Plots a transformed continuous resistance surface against the original resistance values
#'
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. A vector of two parameters is required. The first term is the value of shape parameter (c), and the second term is the value of maximum scale parameter (b)
#' @param Resistance Accepts three types of inputs. Provide either the path to the raw, untransformed resistance surface file or specify an R raster object. Alternatively, supply a vector with the minimum and manximum values (e.g., c(1,10))
#' @param transformation Transformation equation to apply. Can be provided as the name of the transformation or its numeric equivalent (see details)
#' @param scale The standard deviation, in number of raster cells, to use when applying Gaussian kernel smoothing. This is the `sigma` parameter in the `kernel2dsmooth` function. (Default = NULL)
#' @param print.dir Specify the directory where a .tiff of the transformation will be written (Default = NULL)
#' @param marginal.plot Logical. Should distribution plots be added to the margins of the response plot (Default = TRUE). Requires that the resistance surface be specified for `Resistance`.
#' @param marg.type Type of marginal plot to add. One of: [density, histogram, boxplot]. See \code{\link[ggExtra]{ggMarginal}}
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
#' Because of the flexibility of the Ricker function to take a monomolecular shape (try \code{Plot.trans(PARM=c(10,100), Resistance=c(1,10), transformation="Ricker")} to see this), whenever a shape parameter >6 is selected in combination with a Ricker family transformation, the transformation reverts to a Distance transformation. In general, it seems that using a combination of intermediate Ricker and Monomolecular transformations provides the best, most flexible coverage of parameter space. This constraint has not been implemented in the \code{Plot.tans} function.
#' @usage Plot.trans(PARM, 
#'                   Resistance, 
#'                   transformation,
#'                   scale, 
#'                   print.dir,
#'                   marginal.plot, 
#'                   marg.type, 
#'                   Name)
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>

Plot.trans <- function(PARM,
                       Resistance,
                       transformation,
                       scale = NULL,
                       print.dir = NULL,
                       marginal.plot = TRUE,
                       marg.type = "histogram",
                       Name = "layer") {
  if (length(Resistance) == 2) {
    marginal.plot <- FALSE
    # stop(
    #   "If you wish to generate marginal plots, please supply the raster surface to the `Resistance` arguement."
    # )
    
    r <- Resistance
    if (!is.null(scale)) {
      stop(
        "You must provide a raster surface if you wish to plot a scaled transformation"
      )
      # r <- k.smooth(raster = r,
      #               sigma = scale,
      #               SCALE = FALSE)
    #   zmat <- as.matrix(r)
    #   
    #   x <- spatstat::as.im(zmat)
    #   
    #   r <- spatstat::blur(x = x,
    #                       sigma = sigma,
    #                       normalise = TRUE,
    #                       bleed = FALSE)
    }
    Mn = min(r)
    Mx = max(r)
    # Mn = cellStats(r, stat = 'min')
    # Mx = cellStats(r, stat = 'max')
    NAME <- Name
  } else if (class(Resistance)[1] != 'RasterLayer') {
    r <- raster(Resistance)
    NAME <- basename(Resistance)
    NAME <- sub("^([^.]*).*", "\\1", NAME)
    names(r) <- NAME
    
    if (!is.null(scale)) {
      r <- k.smooth(raster = r,
                    sigma = scale,
                    SCALE = FALSE)
      
      # zmat <- as.matrix(r)
      # 
      # x <- spatstat::as.im(zmat)
      # 
      # r <- spatstat::blur(x = x,
      #                     sigma = sigma,
      #                     normalise = TRUE,
      #                     bleed = FALSE)
      
      # Mn = min(r)
      # Mx = max(r)
      Mn = cellStats(r, stat = 'min')
      Mx = cellStats(r, stat = 'max')
    } else {
      Mn = cellStats(r, stat = 'min')
      Mx = cellStats(r, stat = 'max')
    }
    
  } else {
    r <- Resistance
    NAME <- names(r)
    # NAME <- sub("^([^.]*).*", "\\1", NAME)
    # names(r) <- NAME
    
    if (!is.null(scale)) {
      r <- k.smooth(raster = r,
                    sigma = scale,
                    SCALE = FALSE)
      
      # zmat <- as.matrix(r)
      # 
      # x <- spatstat::as.im(zmat)
      # 
      # r <- spatstat::blur(x = x,
      #                     sigma = sigma,
      #                     normalise = TRUE,
      #                     bleed = FALSE)
      # Mn = min(r)
      # Mx = max(r)
      Mn = cellStats(r, stat = 'min')
      Mx = cellStats(r, stat = 'max')
    } else {
      Mn = cellStats(r, stat = 'min')
      Mx = cellStats(r, stat = 'max')
    }
  }
  
  if (Name != "layer") {
    NAME <- Name
  }
  
  # Make vector data --------------------------------------------------------
  
  # * No marginal plots -----------------------------------------------------
  if(marginal.plot == FALSE) {
    original <- seq(from = Mn,
                    to = Mx,
                    length.out = 200)
    dat.t <- SCALE.vector(data = original, 0, 10)
    
    SHAPE <- PARM[[1]]
    Max.SCALE <- PARM[[2]]
    
    if (is.numeric(transformation)) {
      equation <- get.EQ(transformation)
      
    } else {
      equation <- transformation
    }
    # Set equation/name combination
    if (equation == "Distance") {
      Trans.vec <- (dat.t * 0) + 1
      TITLE <- "Distance"
      
    } else if (equation == "Inverse-Reverse Monomolecular") {
      SIGN <- -1 # Inverse
      Trans.vec <-
        SIGN * PARM[[2]] * (1 - exp(-1 * dat.t / PARM[[1]])) + SIGN # Monomolecular
      Trans.vec <-
        rev(SCALE.vector(Trans.vec, MIN = abs(max(Trans.vec)), MAX = abs(min(Trans.vec))))
      TITLE <- "Inverse-Reverse Monomolecular"
      
    } else if (equation == "Inverse Monomolecular") {
      SIGN <- -1 # Inverse
      Trans.vec <-
        SIGN * PARM[[2]] * (1 - exp(-1 * dat.t / PARM[[1]])) + SIGN # Monomolecular
      Trans.vec <-
        (SCALE.vector(Trans.vec, MIN = abs(max(Trans.vec)), MAX = abs(min(Trans.vec))))
      TITLE <- "Inverse Monomolecular"
      
    } else if (equation == "Monomolecular") {
      SIGN <- 1
      Trans.vec <-
        SIGN * PARM[[2]] * (1 - exp(-1 * dat.t / PARM[[1]])) + SIGN # Monomolecular
      TITLE <- "Monomolecular"
      
    } else if (equation == "Reverse Monomolecular") {
      SIGN <- 1
      Trans.vec <-
        SIGN * PARM[[2]] * (1 - exp(-1 * dat.t / PARM[[1]])) + SIGN # Monomolecular
      Trans.vec <- rev(Trans.vec) # Reverse
      TITLE <- "Reverse Monomolecular"
      
    } else if (equation == "Inverse Ricker") {
      SIGN <- -1 # Inverse
      Trans.vec <-
        SIGN * (PARM[[2]] * dat.t * exp(-1 * dat.t / PARM[[1]])) + SIGN #  Ricker
      Trans.vec <-
        SCALE.vector(Trans.vec, MIN = abs(max(Trans.vec)), MAX = abs(min(Trans.vec)))
      TITLE <- "Inverse Ricker"
      
    } else if (equation == "Reverse Ricker") {
      SIGN <- 1
      Trans.vec <-
        SIGN * (PARM[[2]] * dat.t * exp(-1 * dat.t / PARM[[1]])) + SIGN #  Ricker
      Trans.vec <- rev(Trans.vec) # Reverse
      TITLE <- "Reverse Ricker"
      
    } else if (equation == "Inverse-Reverse Ricker") {
      SIGN <- -1 # Inverse
      Trans.vec <-
        SIGN * (PARM[[2]] * dat.t * exp(-1 * dat.t / PARM[[1]])) + SIGN #  Ricker
      Trans.vec <-
        rev(SCALE.vector(Trans.vec, MIN = abs(max(Trans.vec)), MAX = abs(min(Trans.vec)))) # Reverse
      TITLE <- "Inverse-Reverse Ricker"
      
    }  else  {
      SIGN <- 1
      Trans.vec <-
        SIGN * (PARM[[2]] * dat.t * exp(-1 * dat.t / PARM[[1]])) + SIGN #  Ricker
      TITLE <- "Ricker"
    }
    transformed <- Trans.vec
    plot.data <- data.frame(original, transformed)
    x.break <- pretty(original * 1.07)
    y.break <- pretty(transformed * 1.07)
    
    (
      p <- ggplot(plot.data, aes(x = original, y = transformed)) +
        # ggtitle(equation) +
        theme_bw() +
        geom_line(size = 1.5) +
        xlab(expression(bold(
          "Original data values"
        ))) +
        ylab(expression(bold(
          "Transformed data values"
        ))) +
        theme(
          # plot.title = element_text(
          #   lineheight = 2,
          #   face = "bold",
          #   size = 20
          # ),
          legend.title = element_blank(),
          legend.key = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)
        ) +
        scale_x_continuous(limits = c(min(original), max(original)), breaks =
                             x.break) +
        scale_y_continuous(limits = c(min(transformed), max(transformed)), breaks =
                             y.break) +
        removeGrid()
    )
    
    # * Marginal plots --------------------------------------------------------
    
  } else {
    original1 <- seq(from = Mn,
                     to = Mx,
                     length.out = 200)
    dat.t1 <- SCALE.vector(data = original1, 0, 10)
    type1 <- rep(2, 200)
    
    original2 <- sort(as.vector(r))
    dat.t2 <- SCALE.vector(data = original2, 0, 10)
    type <- rep(1, length(original2))
    
    Trans.list <- DAT <- list(dat.t2,
                              dat.t1
    )
    
    
    SHAPE <- PARM[[1]]
    Max.SCALE <- PARM[[2]]
    
    if (is.numeric(transformation)) {
      equation <- get.EQ(transformation)
      
    } else {
      equation <- transformation
    }
    # Set equation/name combination
    for(i in 1:2){
      dat.t <- DAT[[i]]
      if (equation == "Distance") {
        Trans.vec <- (dat.t * 0) + 1
        TITLE <- "Distance"
        
      } else if (equation == "Inverse-Reverse Monomolecular") {
        SIGN <- -1 # Inverse
        Trans.vec <-
          SIGN * PARM[[2]] * (1 - exp(-1 * dat.t / PARM[[1]])) + SIGN # Monomolecular
        Trans.vec <-
          rev(SCALE.vector(Trans.vec, MIN = abs(max(Trans.vec)), MAX = abs(min(Trans.vec))))
        TITLE <- "Inverse-Reverse Monomolecular"
        
      } else if (equation == "Inverse Monomolecular") {
        SIGN <- -1 # Inverse
        Trans.vec <-
          SIGN * PARM[[2]] * (1 - exp(-1 * dat.t / PARM[[1]])) + SIGN # Monomolecular
        Trans.vec <-
          (SCALE.vector(Trans.vec, MIN = abs(max(Trans.vec)), MAX = abs(min(Trans.vec))))
        TITLE <- "Inverse Monomolecular"
        
      } else if (equation == "Monomolecular") {
        SIGN <- 1
        Trans.vec <-
          SIGN * PARM[[2]] * (1 - exp(-1 * dat.t / PARM[[1]])) + SIGN # Monomolecular
        TITLE <- "Monomolecular"
        
      } else if (equation == "Reverse Monomolecular") {
        SIGN <- 1
        Trans.vec <-
          SIGN * PARM[[2]] * (1 - exp(-1 * dat.t / PARM[[1]])) + SIGN # Monomolecular
        Trans.vec <- rev(Trans.vec) # Reverse
        TITLE <- "Reverse Monomolecular"
        
      } else if (equation == "Inverse Ricker") {
        SIGN <- -1 # Inverse
        Trans.vec <-
          SIGN * (PARM[[2]] * dat.t * exp(-1 * dat.t / PARM[[1]])) + SIGN #  Ricker
        Trans.vec <-
          SCALE.vector(Trans.vec, MIN = abs(max(Trans.vec)), MAX = abs(min(Trans.vec)))
        TITLE <- "Inverse Ricker"
        
      } else if (equation == "Reverse Ricker") {
        SIGN <- 1
        Trans.vec <-
          SIGN * (PARM[[2]] * dat.t * exp(-1 * dat.t / PARM[[1]])) + SIGN #  Ricker
        Trans.vec <- rev(Trans.vec) # Reverse
        TITLE <- "Reverse Ricker"
        
      } else if (equation == "Inverse-Reverse Ricker") {
        SIGN <- -1 # Inverse
        Trans.vec <-
          SIGN * (PARM[[2]] * dat.t * exp(-1 * dat.t / PARM[[1]])) + SIGN #  Ricker
        Trans.vec <-
          rev(SCALE.vector(Trans.vec, MIN = abs(max(Trans.vec)), MAX = abs(min(Trans.vec)))) # Reverse
        TITLE <- "Inverse-Reverse Ricker"
        
      }  else  {
        SIGN <- 1
        Trans.vec <-
          SIGN * (PARM[[2]] * dat.t * exp(-1 * dat.t / PARM[[1]])) + SIGN #  Ricker
        TITLE <- "Ricker"
      }
      Trans.list[[i]] <- Trans.vec
    } # End for loop
    
    
    transformed <- unlist(Trans.list)
    plot.data <- data.frame(original = c(original2, original1), 
                            transformed = transformed, 
                            type = c(type, type1))
    x.break <- pretty(original1 * 1.07)
    y.break <- pretty(transformed * 1.07)
 
    p <- ggplot(plot.data[type==1,], aes(x = original, y = transformed)) +
      theme_bw() +
      geom_point(size = 1.5, color = "white") +
      # geom_line(size = 1.5, color = "white") +
      xlab(expression(bold(
        "Original data values"
      ))) +
      ylab(expression(bold(
        "Transformed data values"
      ))) +
      theme(
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
      ) +
      scale_x_continuous(
        breaks = x.break) +
      scale_y_continuous(
        breaks = y.break) +
      geom_line(data = plot.data[plot.data$type==2,],
                aes(x = original, y = transformed),
                size = 1.5,
                color = "black") +
      removeGrid()
    
    # Add marginal plot
    (p <- ggMarginal(p, 
                     type = marg.type,
                     size = 8
    ))
    
    
  } # End response plot
  
  if (!is.null(print.dir)) {
    tiff(
      filename = paste0(print.dir, equation, "_Transformation_", NAME , ".tif"),
      width = 160,
      height = 150,
      units = "mm",
      res = 300,
      compression = "lzw"
    )
    print(p)
    dev.off()
    return(p)
  }
  print(p)
  return(p)
  
}

#### ORIGINAL ####
#' #' Plot continuous surface transformation
#' #'
#' #' Plots a transformed continuous resistance surface against the original resistance values
#' #'
#' #' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. A vector of two parameters is required. The first term is the value of shape parameter (c), and the second term is the value of maximum scale parameter (b)
#' #' @param Resistance Accepts three types of inputs. Provide either the path to the raw, untransformed resistance surface file or specify an R raster object. Alternatively, supply a vector with the minimum and manximum values (e.g., c(1,10))
#' #' @param transformation Transformation equation to apply. Can be provided as the name of the transformation or its numeric equivalent (see details)
#' #' @param print.dir Specify the directory where a .tiff of the transformation will be written (Default = NULL)
#' #' @param Name Name of resistance surface being transformed (optional). This will be added to the output file name.
#' #' @return plot of transformed resistance values against original resistance values
#' #' @details This function will create a ggplot object and plot, so it requires \pkg{ggplot2} to be installed.\cr
#' #' Equation names can be:
#' #' \tabular{ll}{
#' #'    \tab 1 = "Inverse-Reverse Monomolecular"\cr
#' #'    \tab 2 = "Inverse-Reverse Ricker"\cr
#' #'    \tab 3 = "Monomolecular"\cr
#' #'    \tab 4 = "Ricker"\cr
#' #'    \tab 5 = "Reverse Monomolecular"\cr
#' #'    \tab 6 = 'Reverse Ricker"\cr
#' #'    \tab 7 = "Inverse Monomolecular"\cr
#' #'    \tab 8 = "Inverse Ricker"\cr
#' #'    \tab 9 = "Distance"\cr
#' #'    }
#' #'
#' #' The "Distance" equation sets all cell values equal to 1.
#' #' Because of the flexibility of the Ricker function to take a monomolecular shape (try \code{Plot.trans(PARM=c(10,100), Resistance=c(1,10), transformation="Ricker")} to see this), whenever a shape parameter >6 is selected in combination with a Ricker family transformation, the transformation reverts to a Distance transformation. In general, it seems that using a combination of intermediate Ricker and Monomolecular transformations provides the best, most flexible coverasge of parameter space. This constraint has not been implemented in the \code{Plot.tans} function.
#' #' @usage Plot.trans(PARM, Resistance, transformation, print.dir, Name)
#' #' @export
#' #' @author Bill Peterman <Bill.Peterman@@gmail.com>
#'
#' Plot.trans <- function(PARM,Resistance,transformation, print.dir=NULL, Name="layer"){
#'   if(length(Resistance)==2) {
#'     r <- Resistance
#'     Mn=min(r)
#'     Mx=max(r)
#'     NAME <- Name
#'   } else if(class(Resistance)[1]!='RasterLayer') {
#'     r<-raster(Resistance)
#'     NAME <- basename(Resistance)
#'     NAME<-sub("^([^.]*).*", "\\1", NAME)
#'     names(r)<-NAME
#'     Mn=cellStats(r,stat='min')
#'     Mx=cellStats(r,stat='max')
#'
#'   } else {
#'     r <- Resistance
#'     Mn=cellStats(r,stat='min')
#'     Mx=cellStats(r,stat='max')
#'     NAME <- basename(r@file@name)
#'     NAME<-sub("^([^.]*).*", "\\1", NAME)
#'     names(r)<-NAME
#'   }
#'
#'   if(Name!="layer"){
#'     NAME<-Name
#'   }
#'
#'   # Make vector of data
#'   original <- seq(from=Mn,to=Mx,length.out=200)
#'   dat.t <- SCALE.vector(data=original,0,10)
#'
#'   SHAPE <- PARM[[1]]
#'   Max.SCALE <- PARM[[2]]
#'
#'   if(is.numeric(transformation)){
#'     equation<-get.EQ(transformation)
#'
#'   } else {
#'     equation<-transformation
#'   }
#'   # Set equation/name combination
#'   if(equation=="Distance") {
#'     Trans.vec <- (dat.t*0)+1
#'     TITLE <- "Distance"
#'
#'   } else if(equation=="Inverse-Reverse Monomolecular"){
#'     SIGN<- -1 # Inverse
#'     Trans.vec <-  SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN # Monomolecular
#'     Trans.vec <- rev(SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec))))
#'     TITLE <- "Inverse-Reverse Monomolecular"
#'
#'   } else if(equation=="Inverse Monomolecular"){
#'     SIGN<- -1 # Inverse
#'     Trans.vec <-  SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN # Monomolecular
#'     Trans.vec <- (SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec))))
#'     TITLE <- "Inverse Monomolecular"
#'
#'   } else if(equation=="Monomolecular") {
#'     SIGN<- 1
#'     Trans.vec <- SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN # Monomolecular
#'     TITLE <- "Monomolecular"
#'
#'   } else if(equation=="Reverse Monomolecular") {
#'     SIGN<- 1
#'     Trans.vec <-  SIGN*PARM[[2]]*(1-exp(-1*dat.t/PARM[[1]]))+SIGN # Monomolecular
#'     Trans.vec <- rev(Trans.vec) # Reverse
#'     TITLE <- "Reverse Monomolecular"
#'
#'   } else if (equation=="Inverse Ricker") {
#'     SIGN <- -1 # Inverse
#'     Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
#'     Trans.vec <- SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec)))
#'     TITLE <- "Inverse Ricker"
#'
#'   } else if (equation=="Reverse Ricker") {
#'     SIGN <- 1
#'     Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
#'     Trans.vec <- rev(Trans.vec) # Reverse
#'     TITLE <- "Reverse Ricker"
#'
#'   } else if (equation=="Inverse-Reverse Ricker") {
#'     SIGN <- -1 # Inverse
#'     Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
#'     Trans.vec <- rev(SCALE.vector(Trans.vec,MIN=abs(max(Trans.vec)),MAX=abs(min(Trans.vec)))) # Reverse
#'     TITLE <- "Inverse-Reverse Ricker"
#'
#'   }  else  {
#'     SIGN <- 1
#'     Trans.vec <- SIGN*(PARM[[2]]*dat.t*exp(-1*dat.t/PARM[[1]]))+SIGN #  Ricker
#'     TITLE <- "Ricker"
#'   }
#'   transformed<-Trans.vec
#'   plot.data<-data.frame(original,transformed)
#'   x.break<-pretty(original*1.07)
#'   y.break<-pretty(transformed*1.07)
#'
#'   ( p<- ggplot(plot.data,aes(x=original,y=transformed)) +
#'       ggtitle(equation) +
#'       theme_bw() +
#'       geom_line(size=1.5) +
#'       xlab(expression(bold("Original data values"))) +
#'       ylab(expression(bold("Transformed data values"))) +
#'       theme(plot.title = element_text(lineheight=2, face="bold",size=20),
#'             legend.title = element_blank(),
#'             legend.key = element_blank(),
#'             axis.text.x = element_text(size=14),
#'             axis.text.y = element_text(size=14),
#'             axis.title.x = element_text(size=16),
#'             axis.title.y = element_text(size=16)
#'             # axis.line = element_line(colour = "black"),
#'             # panel.grid.major = element_blank(),
#'             # panel.grid.minor = element_blank(),
#'             # panel.border = element_blank(),
#'             # panel.background = element_blank()
#'       ) +
#'       scale_x_continuous(limits=c(min(original),max(original)),breaks=x.break) +
#'       scale_y_continuous(limits=c(min(transformed),max(transformed)),breaks=y.break)
#'   )
#'
#'   if(!is.null(print.dir)){
#'     tiff(filename=paste0(print.dir,equation,"_Transformation_",NAME ,".tif"),width=160,height=150,units="mm",res=300,compression="lzw")
#'     print(p)
#'     dev.off()
#'     return(p)
#'   }
#'   print(p)
#'   return(p)
#'
#' }