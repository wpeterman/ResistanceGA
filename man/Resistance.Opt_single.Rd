\name{Resistance.Opt_single}
\alias{Resistance.Opt_single}
\title{Optimize resistance surfaces individually}
\usage{
Resistance.Opt_single(PARM, Resistance, CS.inputs, GA.inputs, Min.Max = "max",
  iter, quiet = FALSE)
}
\arguments{
  \item{PARM}{Parameters to transform conintuous surface or
  resistance values of categorical surface. Should be a
  vector with parameters specified in the order of
  resistance surfaces.These values are selected during
  optimization if called within GA function.}

  \item{Resistance}{Resistance surface to be optimized.
  This should be an R raster object. If not specified, the
  function will attempt to find the a resistance surface
  from \code{GA.inputs}}

  \item{CS.inputs}{Object created from running
  \code{\link[ResistanceGA]{CS.prep}} function}

  \item{GA.inputs}{Object created from running
  \code{\link[ResistanceGA]{GA.prep}} function}

  \item{Min.Max}{Define whether the optimization function
  should minimized ('min') or maximized ('max'). Default in
  'max'}

  \item{iter}{A counter for the number of surfaces that
  will be optimized}

  \item{quiet}{Logical, if TRUE AICc and iteration duration
  will not be printed to the screen at the completion of
  each iteration.}
}
\value{
AIC value from mixed effect model
}
\description{
Optimize all resistance surfaces that are located in the
same directory individually. This optimization function is
designed to be called from GA
}

