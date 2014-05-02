\name{Combine_Surfaces}
\alias{Combine_Surfaces}
\title{Combine multiple resistance surfaces together}
\usage{
Combine_Surfaces(PARM, CS.inputs, GA.inputs)
}
\arguments{
  \item{PARM}{Parameters to transform conintuous surface or
  resistance values of categorical surface. Should be a
  vector with parameters specified in the order of
  resistance surfaces}

  \item{CS.inputs}{Object created from running
  \code{\link{CS.prep}} function}

  \item{GA.inputs}{Object created from running
  \code{\link{GA.prep}} function}
}
\value{
R raster object
}
\description{
Combine multiple resistance surfaces into new composite
surface based on specified parameters
}

