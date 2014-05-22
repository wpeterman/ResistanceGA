\name{Resistance.tran}
\alias{Resistance.tran}
\title{Apply transformation to continuous resistance surface}
\usage{
Resistance.tran(transformation, shape, max,CS.inputs,r)
}
\arguments{
  \item{transformation}{Transformation equation to apply.
  Can be provided as the name of the transformation or its
  numeric equivalent (see details)}

  \item{shape}{Value of the shape parameter}

  \item{max}{Value of the maximum value parameter}

  \item{r}{Resistance surface to be transformed. Can be
  supplied as full path to .asc file or as a raster object}

  \item{out}{Directory to write transformed .asc file.
  Default is NULL, and will not export .asc file}
}
\value{
R raster object
}
\description{
Apply on the six resistance transformations to a continuous
resistance surface
}
\details{
Valid arguements for \code{transformation} are:\cr
\tabular{ll}{ \tab 1 = "Inverse-Reverse Monomolecular"\cr
\tab 2 = "Inverse-Reverse Ricker"\cr \tab 3 =
"Monomolecular"\cr \tab 4 = "Ricker"\cr \tab 5 = "Reverse
Monomolecular"\cr \tab 6 = 'Reverse Ricker"\cr \tab 7 =
"Inverse Monomolecular"\cr \tab 8 = "Inverse Ricker"\cr
\tab 9 = "Distance"\cr }
}

