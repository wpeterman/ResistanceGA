\name{PLOT.trans}
\alias{PLOT.trans}
\title{Plot continuous surface transformation}
\usage{
PLOT.trans(PARM, Resistance, equation, print.dir)
}
\arguments{
  \item{PARM}{Parameters to transform conintuous surface or
  resistance values of categorical surface. A vector of two
  parameters is required. The first term isthe value of
  shape parameter (c), and the second term is the value of
  maximum scale parameter (b)}

  \item{Resistance}{Path to raw, untransformed resistance
  surface}

  \item{equation}{Name of the transformation equation to
  use:\cr "Inverse-Reverse Monomolecular"\cr"Inverse
  Monomolecular"\cr"Monomolecular"\cr"Reverse
  Monomolecular"\cr"Inverse
  Ricker"\cr"Ricker",\cr"Distance"}

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
can be "Inverse-Reverse Monomolecular", "Inverse
Monomolecular", "Monomolecular", "Reverse Monomolecular",
"Inverse Ricker", "Ricker", or "Distance". The "Distance"
equation sets all cell values equal to 1.
}

