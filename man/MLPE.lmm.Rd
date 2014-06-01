\name{MLPE.lmm}
\alias{MLPE.lmm}
\title{Run maximum likelihood population effects mixed effects model (MLPE)}
\usage{
MLPE.lmm(cs.resistance, pairwise.genetic, REML)
}
\arguments{
  \item{cs.resistance}{Pairwise resistance distance matrix
  (resistances.out) from CS results}

  \item{pairwise.genetic}{Lower half of pairwise genetic
  distance matrix}

  \item{REML}{Logical. If TRUE, mixed effects model will be
  fit using restricted maximum likelihood. Default = FALSE}
}
\value{
A glmer object from the fitted model
}
\description{
Runs MLPE as detailed by Clarke et al. (2002). This
function will run the model and return glmer object
}
\details{
An AIC value will only be returned if \code{REML = FALSE}
}
\references{
Clarke, R. T., P. Rothery, and A. F. Raybould. 2002.
Confidence limits for regression relationships between
distance matrices: Estimating gene flow with distance.
Journal of Agricultural, Biological, and Environmental
Statistics 7:361-372.
}

