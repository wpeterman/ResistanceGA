\name{Resistance.Opt_single}
\alias{Resistance.Opt_single}
\title{Optimize resistance surfaces individually}
\usage{
Resistance.Opt_single(PARM, Resistance, CS.inputs, GA.inputs, Min.Max, iter)
}
\arguments{
  \item{PARM}{Parameters to transform conintuous surface or
  resistance values of categorical surface. Should be a
  vector with parameters specified in the order of
  resistance surfaces.These values are selected during
  optimization if called within GA function.}

  \item{CS.inputs}{Object created from running
  \code{\link{CS.prep}} function}

  \item{GA.inputs}{Object created from running
  \code{\link{GA.prep}} function}

  \item{Min.Max}{Define whether the optimization function
  should minimized ('min') or maximized ('max')}
}
\value{
AIC value from mixed effect model
}
\description{
Optimize all resistance surfaces that are located in the
same directory individually. This optimization function is
designed to be called from GA
}

