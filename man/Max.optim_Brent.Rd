\name{Max.optim_Brent}
\alias{Max.optim_Brent}
\title{Determine maximum value of resistance surfaces using Brent optimization algorithm}
\usage{
Max.optim_Brent(PARM, CS.inputs, GA.inputs, Min.Max = "min", quiet = FALSE,
  GA.opt)
}
\arguments{
  \item{PARM}{A single multiplier value to pass to Brent
  optimization algorithm. Can range from 0--25}

  \item{CS.inputs}{Object created from running
  \code{\link[ResistanceGA]{CS.prep}} function}

  \item{GA.inputs}{Object created from running
  \code{\link[ResistanceGA]{GA.prep}} function}

  \item{Min.Max}{Define whether the optimization function
  should minimized ('min') or maximized ('max')}

  \item{quiet}{Logical, if TRUE, AICc and iteration time
  will not be printed to the screen at the completion of
  each iteration. Default = FALSE}

  \item{GA.opt}{Optimied parameters from GA optimization}
}
\value{
AIC value from mixed effect model
}
\description{
Following GA optimization, determine if the maximum
resistance values are correct. This optimization function
is designed to take results from GA optimization
}

