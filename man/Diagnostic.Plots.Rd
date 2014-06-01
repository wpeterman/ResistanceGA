\name{Diagnostic.Plots}
\alias{Diagnostic.Plots}
\title{Create diagnostic plots}
\usage{
Diagnostic.Plots(cs.resistance.mat, genetic.dist, XLAB,YLAB, plot.dir)
}
\arguments{
  \item{cs.resistance.mat}{Path to CIRCUITSCAPE
  "_resistances.out" file}

  \item{genetic.dist}{Vector of pairwise genetic distances
  (lower half of pairwise matrix). Can be input as
  CS.inputs$RESPONSE}

  \item{XLAB}{Label for x-axis (Defaults to "Estimated
  resistance")}

  \item{YLAB}{Label for y-axis (Defaults to "Genetic
  distance")}

  \item{plot.dir}{Directory to output PDF of diagnostic
  plots}
}
\value{
A four panel PDF including residual scatterplot, historgram
of residuals, qqplot, and
}
\description{
This function will generate mixed effect model diagnostic
plots following optimization
}

