\name{PLOT.trans}
\alias{PLOT.trans}
\title{Plot continuous surface transformation}
\usage{
PLOT.trans(PARM, Resistance, equation, print)
}
\arguments{
  \item{PARM}{Parameters to transform conintuous surface or
  resistance values of categorical surface. A vector with
  the first term identifying the value of shape parameter
  (c), and the second term identifying the maximum scale
  parameter (b)}

  \item{Resistance}{Path to raw, untransformed resistance
  surface}

  \item{equation}{Name of the transformation equation to
  use (see details below)}

  \item{print}{Logical. If TRUE, a PDF file of the plot
  will exported to the current working directory (Default =
  FALSE)}
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
Equation names can be "Inverse-Reverse Monomolecular",
"Inverse Monomolecular", "Monomolecular", "Reverse
Monomolecular", "Inverse Ricker" , "Ricker", or "Distance"
}

