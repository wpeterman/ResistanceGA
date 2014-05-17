\name{SS_optim}
\alias{SS_optim}
\title{Single surface optimization}
\usage{
SS_optim(CS.inputs, GA.inputs)
}
\arguments{
  \item{CS.inputs}{Object created from running
  \code{\link{CS.prep}} function}

  \item{GA.inputs}{Object created from running
  \code{\link{GA.prep}} function}

  \item{nlm}{Logical, if TRUE, the final step of
  optimization will use nlm to fine-tune parameter
  estimates. This may lead to overfitting in some cases.
  Default = FALSE.}
}
\value{
This function optimizes multiple resistance surfaces.
Following optimization, several summary objects are
created.\cr \enumerate{ \item Diagnostic plots of model fit
are output to the "Results/Plots" directory that is
automatically generated within the folder containing the
optimized ASCII files. \item A .csv file with the Maximum
Likelihood Population Effects mixed effects model
coefficient estimates (MLPE_coeff_Table.csv) \item Three
summary .csv files are generated: CategoricalResults.csv,
ContinuousResults.csv, & All_Results_AICc.csv. These tables
contain AICc values and optimization summaries for each
surface. } All results tables are also summarized in a
named list ($ContinuousResults, $CategoricalResults, $AICc,
$MLPE)
}
\description{
Optimize all surfaces contained in a directory using a
genetic algorithm executed with the \code{\link[GA]{ga}}
function in the Genetic Algorithms package \pkg{GA}
}

