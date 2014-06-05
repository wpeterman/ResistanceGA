\name{Combine_Surfaces}
\alias{Combine_Surfaces}
\title{Combine multiple resistance surfaces together}
\usage{
Combine_Surfaces(PARM, CS.inputs, GA.inputs, out = GA.inputs$Results.dir)
}
\arguments{
  \item{PARM}{Parameters to transform conintuous surface or
  resistance values of categorical surface. Requires a
  vector with parameters specified in the order of
  resistance surfaces}

  \item{CS.inputs}{Object created from running
  \code{\link[ResistanceGA]{CS.prep}} function}

  \item{GA.inputs}{Object created from running
  \code{\link[ResistanceGA]{GA.prep}} function}

  \item{out}{Directory to write combined .asc file. Default
  is GA.inputs$Results.dir}
}
\value{
R raster object that is the sum all transformed and/or
reclassified resistance surfaces provided
}
\description{
Combine multiple resistance surfaces into new composite
surface based on specified parameters
}
\details{
\code{PARM} is designed to accept the output of
\code{MS_optim}. For continuous surfaces, there are three
terms: 1) Transformation, 2) shape, and 3) maximum value.
Transformation must be provided as a numeric value:\cr
\tabular{ll}{ \tab 1 = "Inverse-Reverse Monomolecular"\cr
\tab 2 = "Inverse-Reverse Ricker"\cr \tab 3 =
"Monomolecular"\cr \tab 4 = "Ricker"\cr" \tab 5 = "Reverse
Monomolecular"\cr \tab 6 = "Reverse Ricker"\cr \tab 7 =
"Inverse Monomolecular"\cr \tab 8 = "Inverse Ricker"\cr
\tab 9 = "Distance"\cr }
}

