\name{Resistance.tran}
\alias{Resistance.tran}
\title{Apply transformation to continuous resistance surface}
\usage{
Resistance.tran(transformation, shape, max,CS.inputs,r)
}
\arguments{
  \item{transformation}{Transformation equation to apply}

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
Valid arguements for \code{Transformation} are:\cr
"Inverse-Reverse Monomolecular"\cr"Inverse
Monomolecular"\cr"Monomolecular"\cr"Reverse
Monomolecular"\cr"Inverse Ricker"\cr"Ricker",\cr"Distance"
}

