\name{PLOT.trans}
\alias{PLOT.trans}
\title{Plot continuous surface transformation}
\usage{
PLOT.trans(PARM, Resistance, transformation, print.dir)
}
\arguments{
  \item{PARM}{Parameters to transform conintuous surface or
  resistance values of categorical surface. A vector of two
  parameters is required. The first term isthe value of
  shape parameter (c), and the second term is the value of
  maximum scale parameter (b)}

  \item{Resistance}{Accepts two types of inputs. Provide
  either the path to the raw, untransformed resistance
  surface file or specify an R raster object}

  \item{transformation}{Transformation equation to apply.
  Can be provided as the name of the transformation or its
  numeric equivalent (see details)}

  \item{print.dir}{Specify the directory where a .tiff of
  the transformation will be written (Default = NULL)}
}
\value{
plot of transformed resistance values against original
resistance values
}
\description{
Plots a transformed continuous resistance surface against
the original resistance values
}
\details{
This function will create a ggplot object and plot, so it
requires \pkg{ggplot2} to be installed.\cr Equation names
can be: #' \tabular{ll}{ \tab 1 = "Inverse-Reverse
Monomolecular"\cr \tab 2 = "Reverse Monomolecular"\cr \tab
3 = "Monomolecular"\cr \tab 4 = "Inverse Monomolecular"\cr
\tab 5 = "Inverse Ricker"\cr \tab 6 = 'Ricker"\cr \tab 7 =
"Reverse Ricker"\cr \tab 8 = "Inverse-Reverse Ricker"\cr
\tab 9 = "Distance"\cr }

The "Distance" equation sets all cell values equal to 1.
}

