\name{GA.prep}
\alias{GA.prep}
\title{Create R object with genetic algorithm optimization settings}
\usage{
GA.prep(ASCII.dir,
Min.Max="max",
min.cat=0,
max.cat=2500,
max.cont=2500,
cont.shape=NULL,
pop.mult = 15,
percent.elite = 0.05,
type= "real-valued",
pcrossover=0.85,
pmutation=0.1,
maxiter=1000,
run=25,
keepBest=TRUE,
population = gaControl(type)$population,
selection = gaControl(type)$selection,
crossover="gareal_blxCrossover",
mutation = gaControl(type)$mutation)
}
\arguments{
  \item{ASCII.dir}{Directory containing all raster objects
  to optimized.}

  \item{min.cat}{The minimum value to be assessed during
  optimization of of categorical resistance surfaces
  (Default = 0)}

  \item{max.cat}{The maximum value to be assessed during
  optimization of of categorical resistance surfaces
  (Default = 2500)}

  \item{max.cont}{The maximum value to be assessed during
  optimization of of continuous resistance surfaces
  (Default = 2500)}

  \item{cont.shape}{A vector of hypothesized relationships
  that each continuous resistance surface will have in
  relation to the genetic distance reposnse (Default =
  NULL; see details)}

  \item{Neighbor.Connect}{Select 4 or 8 to designate the
  connection scheme to use in CIRCUITSCAPE (Default = 8)}

  \item{pop.mult}{Value will be multiplied with number of
  parameters in surface to determine 'popSize' in GA. By
  default this is set to 15.}

  \item{percent.elite}{Percent used to determine the number
  of best fitness individuals to survive at each generation
  ('elitism' in GA). By default the top 5\% individuals
  will survive at each iteration.}

  \item{type}{Default is "real-valued"}

  \item{population}{Default is gareal_Population from GA}

  \item{selection}{Default is gareal_lsSelection from GA}

  \item{mutation}{Default is gareal_raMutation from GA}

  \item{pcrossover}{Probability of crossover. Default =
  0.85}

  \item{pmutation}{Probability of mutation. Default = 0.1}

  \item{crossover}{Default = "gareal_blxCrossover". This
  crossover method greatly improved optimization during
  preliminary testing}

  \item{maxiter}{Maximum number of iterations to run before
  the GA search is halted (Default = 1000)}

  \item{run}{Number of consecutive generations without any
  improvement in AICc before the GA is stopped (Default =
  25)}

  \item{keepBest}{A logical argument specifying if best
  solutions at each iteration should be saved (Default =
  TRUE)}

  \item{Min.Max}{Define whether the optimization function
  should minimized ('min') or maximized ('max' = Default).
  Optimization with \code{ga} maximizes the objective
  criteria}
}
\value{
An R object that is a required input into optimization
functions
}
\description{
This function prepares and compiles objects and commands
for optimization with the GA package
}
\details{
Only files that you wish to optimize, either in isolation
or simultaneously, should be included in the specified
\code{ASCII.dir}. If you wish to optimize different
combinations of surfaces, different directories contaiing
these surfaces must be created.

\code{cont.shape} can take values of "Increase",
"Decrease", or "Peaked". If you believe a resistance
surface is related to your reposnse in a particular way,
specifying this here may decrease the time to optimization.
\code{cont.shape} is used to generate an initial set of
parameter values to test during optimization. If specified,
a greater proportion of the starting values will include
your believed relatiosnship. If unspecified (the Default),
a completely random set of starting values will be
generated.

It is recommended to first run GA optimization with the
default settings
}

